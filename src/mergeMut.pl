use strict;
use File::Glob ':glob';
use Data::Dumper;
use Getopt::Long;
use FindBin qw($RealBin);
use File::Basename;
use List::Util qw(max min sum);


my $list;   #filename of all vcfs
my $type;
my $normal;
my $task;
my $qualfilter;
my $prefix;
my $dbsnp = "no";
my $tmpdir = "./";
my $bin = $RealBin;
my $tolerance = 0;
my $withChr;
my $strandBiasTh = 0.005;
my $tailDisBiasTh = 0.005;
my $nonsegdup;
my $exonic;
my $qualTitan = 40;
my $vafTitan = 0.15;

GetOptions (
            "list|l=s"       => \$list,             #filename of all vcfs
            "type|t=s"       => \$type,             #snv or indel
            "normal|n=s"     => \$normal,           #comma seperated id of normal samples
            "task|k=s"       => \$task,             #task type
            "prefix|p=s"     => \$prefix,
            "qualfilter"     => \$qualfilter,
            "dbsnp|d=s"      => \$dbsnp,
            "tmpdir|y=s"     => \$tmpdir,
            "tolerance=i"    => \$tolerance,
            "withChr"        => \$withChr,
            "strandBiasTh=f" => \$strandBiasTh,
            "tailDisBiasTh=f"=> \$tailDisBiasTh,
            "nonsegdup"      => \$nonsegdup,
            "exonic"         => \$exonic,
            "qualTitan=i"    => \$qualTitan,
            "vafTitan=f"     => \$vafTitan,
            "help|h"         => sub{
                               print "usage: $0 get all somatic and rare variants from a bunch of vcf files\n\nOptions:\n\t--list\t\tthe filename of all vcfs\n";
                               print "\t--type\t\tthe type of variants, snv or indel\n";
                               print "\t--normal\tcomma seperated id of normal samples\n";
                               print "\t--prefix\tthe prefix of samples' names\n";
                               print "\t--qualfilter\tset the number when you want to skip low qual calls\n";
                               print "\t--task\t\tthe task, such as tcga or rnaediting\n";
                               print "\t--dbsnp\t\tyes or no, whether to keep dbsnp variants into the table\n";
                               print "\t--tolerance\tthe tolerance for comparing indels with certain distance shift\n";
                               print "\t--withChr\tconvert all chromosome names to be started with chr\n";
                               print "\t--tmpdir\tthe temporary dir to write tmp files\n";
                               print "\t--strandBiasTh\tthe p value threshold for filtering out vars with strand bias <0.005>\n";
                               print "\t--tailDisBiasTh\tthe p value threshold for filtering out vars with Tail distance bias <0.005>\n";
                               print "\t--nonsegdup\tdo not collect segdup ones\n";
                               print "\t--exonic\tonly collect exonic ones (including UTR and splicing)\n";
                               print "\t--qualTitan\tmin qual for titan selection of heteroSNPs\n";
                               print "\t--vafTitan\tmin VAF for titan selection of heteroSNPs\n";
                               print "\t--help\t\tprint this help message\n";
                               print "\n";
                               exit 0;
                             },
           );


#define needed chrs
my %chrs;
for (1..22){
  $chrs{$_} = '';
  $chrs{'chr'.$_} = '';
}
$chrs{'chrX'} = '';
$chrs{'chrY'} = '';
$chrs{'chrM'} = '';
$chrs{'chrMT'} = '';
$chrs{'X'} = '';
$chrs{'Y'} = '';
$chrs{'M'} = '';
$chrs{'MT'} = '';
print STDERR Dumper(\%chrs);

my @list;
open IN, "$list";
while ( <IN> ) {
  chomp;
  next if /^#/;
  push(@list, $_);
}
close IN;

my @prefix = split(',', $prefix);
my $prefixReg = join('|', @prefix);
print STDERR "prefixReg is $prefixReg\n";

my %normals;
foreach my $normalControl (split(',', $normal)){
  $normals{$normalControl} = '';
}
print STDERR Dumper(\%normals);

my $taskOri = $task;
my %somatic;
my %samples;
foreach my $file (@list) {
  $task = $taskOri;
  my $name;
  my $filebase = basename($file);
  if ($prefixReg ne '' and $filebase =~ /(($prefixReg)([a-zA-Z0-9\-\_]+)?)[^a-zA-Z0-9\-\_]/) {     #changed to adapt to lowGI data
    $name = $1;
  }
  elsif ($task =~ /tcga/i) {
    $filebase =~ /(TCGA\-[A-Z0-9]+\-[A-Z0-9]+)/;
    $name = $1;
  }
  print STDERR "$name\t$file\t$filebase\n";

  $samples{$name} = '';

  my $openway = $file;
  if ($file =~ /\.gz$/) {
    $openway = "gzip -dc $file |";
  }
  print STDERR "$openway\n";

  open IN, "$openway";
  my $revertornot = "no";
  my $printerror = 0;
  my $singlecalling = "no";
  while ( <IN> ) {
     chomp;

     if ($task =~ /mutectTemp/) {  #gathering information for mutect calls "nopara.txt"
       next if $_ =~ /^chr\tpos/;
       my ($chr, $pos, $ref, $alt, $normal_reads, $tumor_reads, $normal_varfrac, $tumor_varfrac, $tumor_nonadj_varfrac, $purity) = split /\t/;
       my $coor = $chr.':'.$pos;
       my $id = '.';
       $somatic{$coor}{$name} = $tumor_varfrac.'|'.$tumor_reads;
       my $idrefalt = join("\t", ($id,$ref,$alt));
       $somatic{$coor}{'info'}{$idrefalt} = 'unknown';  #not necessarily just one!
       $somatic{$coor}{'somatic'} .= $name.',';
       next;
     }

     #judge whether it is single calling or paired calling
     if ($_ =~ /^#/) {
       if ($_ =~ /^#CHROM\tPOS\tID/) {   #the common three column header in vcf file
         my @cols = split /\t/;
         my $minusI = 1;
         if (exists($normals{$cols[$#cols - $minusI]}) or ($cols[$#cols - $minusI] =~ /NORMAL/i)) {
           $revertornot = "yes";
         } elsif ( $cols[$#cols - $minusI] eq 'FORMAT' ) {
           $singlecalling = "yes";
           if ( ! exists($normals{$cols[$#cols]}) ) {  #it is not normal, then start db hetero germline guessing or collecting rare variants only
             if ($task =~ /titan/) {
               $task .= ",guessNormal";
             } elsif ($task =~ /muTect/) {
               $task = "rare";
             }
           }
         }
         print STDERR "revert or not: $revertornot\n";
         print STDERR "singlecalling: $singlecalling\n";
         print STDERR "taskNow: $task\n";
         next;
       } else {
         next;
       }
     }
     #judging whether it is single calling or paired calling

     my ($chr, $pos, $id, $ref, $alt, $qual, $pass, $info, $format, $sample, $blood) = split /\t/;
     $chr =~ s/^chr//;
     if ($withChr) {
       $chr = 'chr'.$chr if ($chr !~ /^chr/);
       $chr = 'chrM' if ($chr eq 'chrMT');
     }

     #chr filter
     if (! exists($chrs{$chr})) {
       next;
     }

     if ($qualfilter) {             #if do qual filter
       next if ($qual ne '.' and $qual < 30 and $pass ne 'PASS');
     }

     if ($task =~ /titan/ or $task =~ /errorEst/) {     #for titan, pick good ones
       next if ($qual ne '.' and $qual < $qualTitan);
     }

     if ( $task =~ /vardict/ ) {
       next if ( $pass ne 'PASS' );
     }

     if ($revertornot eq 'yes') {   #revert sample and blood
        my $tmp = $sample;
        $sample = $blood;
        $blood = $tmp;
     }

     my @formats = split(':', $format);
     my %formindex;
     for(my $f = 0; $f <= $#formats; $f++) {
       $formindex{$formats[$f]} = $f;
     }
     if ($printerror == 0) {
       print STDERR Dumper(\%formindex);
       $printerror ++;
     }

     my @sample = split(/\:/,$sample);


     ###########################################################################decide somatic
     my $somatic = 0;
     if ($singlecalling eq 'no') {      #if it is paired calling
       if ($type eq 'snv') {            #for snp
         if ($task =~ /muTect/) {       #if mutect
           if ($pass eq 'QSI_ref' or $pass eq 'PASS') {
             $somatic = 1;
           }
         } else {
           my @blood = split(/\:/,$blood);
           if ($blood[$formindex{'GT'}] !~ /1/) {
             $somatic = 1;
           }
         }
       } elsif ($type eq 'indel') { #for indel
         if ($task =~ /strelka/) {   #if strelka
           if ($pass =~ /PASS/) {
             $somatic = 1;
           }
         } else {
           my @blood = split(/\:/,$blood);
           if ($blood[$formindex{'GT'}] !~ /1/) {
             $somatic = 1;
           }
         }
       }
     }
     ############################################################################# Decide Somatic, Only for SNV

     if ( $type eq 'indel' and $task =~ /strelka/ ) {                    #if clinvar, skip normal filter
       goto PRODUCE;
     }
     if ( $task =~ /errorEst/ ) {
       if ($somatic == 1 or $id eq '.') {                                #using known vars for estimation of errors
         next;
       }
       goto PRODUCE;
     }


     my $popFreqMax = -1;
     if ($id ne '.' or $info =~ /dbSNP/ or $info =~ /1KG\=/ or $info =~ /1000g[0-9a-z\_\-]+\=/ or $info =~ /ESP\d+\=/ or $info =~ /ExAC\_ALL\=/ or $info =~ /PopFreqMax\=/) {  #snp in population, is it a somatic one?

       my $freq = -1;
       my $freq_eac = -1;
       if ($info =~ /(1KG=(.+?));/) {
         my $kid = $1;                                 #re define $id to 1KG when absent
         $freq = $2;
         if ($id !~ /^[rR][sS]/) {
           $id = $kid;
         }
       }
       if ($info =~ /(1000g[0-9a-z\_\-]+\=(.+?));/) {
         my $kid = $1;                                 #re define $id to 1KG when absent
         $freq = $2;
         if ($id !~ /^[rR][sS]/) {
           $id = $kid;
         }
       }
       if ($info =~ /((ESP\d+)=(.+?));/) {
           my $kid = $1;                               #re define $id to ESP when absent
           $freq = $3;
           if ($id !~ /^[rR][sS]/) {
              $id = $kid;
           }
       }
       if ($info =~ /(ExAC\_ALL=(.+?));/) {
         my $kid = $1;                                 #re define $id to ExAC when absent
         $freq = $2;
         $freq_eac = $2;
         if ($id !~ /^[rR][sS]/) {
           $id = $kid;
         }
       }
       if ($info =~ /(PopFreqMax=(.+?));/) {
         my $kid = $1;                                 #re define $id to PopFreqMax when TOTALLY absent
         my $popmax = $2;
         if ( $popmax > 0.005 and $popmax <= 0.022 and $freq_eac != -1 and $freq_eac < 0.005 ) {                         #save some with rare EAC
           $info =~ s/PopFreqMax\=$popmax\;/PopFreqMax\=$freq_eac\;/;
         } else {
           $freq = $popmax;
         }
         if ($id eq '.') {
           $id = $kid;
         }
       }

       $popFreqMax = $freq;

       if ( $info =~ /clinvar/ ) {                    #if clinvar, skip normal filter
         goto PRODUCE;
       }

       if ( $somatic == 0 ) {                         #if it is not somatic, then only rare ones should be kept
         if ($freq == -1) {
           if ($dbsnp eq "no" or $task =~ /rare/i) {
             next;
           }
         } else {   #freq is defined
           if ($id =~ /^[rR][sS]/) {                  #dbSNP ones with reported MAF (with dbSNP id)
             next if ($dbsnp eq "no");
           }
           if ($task =~ /rare/i) {
             next if $freq > 0.005;                    #rare kept
           }
         }
       }

     } #known variation

   PRODUCE:

     unless ($qual ne '.' and $qual > 90) { #unless very high quality

       if ($info =~ /MQ0F=(.+?);/) {    #if with MQ0F
         next if ($1 > 0.1);
       }

       if ($info =~ /\;PV4\=(.+?)\,(.+?)\,(.+?)\,(.+?);/) {    #if with PV4 info
         my $strandb = $1;
         my $baseqb = $2;
         my $mapqb = $3;
         my $tailb = $4;
         if ($strandb =~ /e/) {               #strand bias
           next;
         } elsif ($strandb < $strandBiasTh) { #strand bias make it very stringent! for net data
           next;
         #} elsif ($baseqb =~ /e/) {           #basequality bias
         #  next;
         #} elsif ($baseqb < 0.0005) {         #basequality bias
         #  next;
         #} elsif ($mapqb =~ /e/) {         #mapquality bias   #mask it for sid's data
         #  next;
         #} elsif ($mapqb < 0.0001) {       #mapquality bias   #mask it for sid's data
         #  next;/MQ0
         } elsif ($tailb =~ /e/) {            #tailbias
           next;
         } elsif ($tailb < $tailDisBiasTh) {  #tailbias
           next;
         } else {
           #pass
         }
       }

     }  #unless very high quality

     my $coor = $chr.':'.$pos;

     my $maf = -1;    #get maf
     my $tdp = -1;    #get total depth
     my $altd = -1;
     if ( $type eq 'indel' and $task =~ /strelka/ ) {

       #GT:DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50

       my @sampleInfo = split(":", $sample);
       if ( $sampleInfo[$formindex{'TIR'}] =~ /^(\d+)\,(\d+)$/ ) {
         $altd = $2;
       } else {
         die ("strelka vcf record without TIR, $_\n");
       }

       if ( $sampleInfo[$formindex{'DP2'}] =~ /(\d+)/ ) {
         $tdp = $1;
       } else {
         die ("strelka vcf record without DP2, $_\n");
       }

       if ( $info =~ /\;QSI\=(\d+)\;/ ) {
         $qual = $1;
       } else {
         die ("strelka vcf record without QSI, $_\n");
       }

       $maf = ($tdp > 0)? sprintf("%.3f", $altd/$tdp) : 0.;

     } else {  #SNV

       if ($info =~ /DP4\=(\d+)\,(\d+)\,(\d+)\,(\d+)\;/) { #DP4 information
         $maf = sprintf("%.3f", ($3+$4)/($1+$2+$3+$4));
         $tdp = $1+$2+$3+$4;
         $altd = $3+$4;
       } else {                 #no DP4 information
         if ($info =~ /DP\=(\d+)\;/) {
           $tdp = $1;
         }
         my @sampleinfo = split(":", $sample);
         if ($sampleinfo[$formindex{'AD'}] =~ /^(\d+)\,(\d+)$/) {
           $maf = (($1+$2) > 0)? sprintf("%.3f", $2/($1+$2)) : 0;
           if ($tdp == -1) {
             $tdp = $1+$2;
           }
           if ($altd == -1) {
             $altd = $2;
           }
         } else {               #multiple depths
           my @ads = split(',', $sampleinfo[$formindex{'AD'}]);
           my $adsum = sum(@ads);
           shift(@ads);
           my $altadsum = sum(@ads);
           $maf = sprintf("%.3f", $altadsum/$adsum);
           if ($tdp == -1) {
             $tdp = $adsum;
           }
           if ($altd == -1) {
             $altd = $altadsum;
           }
         }
       }

     }

     #depth filter
     unless ($qualfilter) {
       unless ( $somatic == 1 ) {
         next if ( $tdp < 5 );
         next if ( $altd < 2 );
         #next if ( $maf < 0.001 );
       }
     }
     #depth filter

     if ( $task =~ /errorEst/ ) {
       next if ( $maf < 0.2 or $maf > 0.8 );
     }
     if ( $task =~ /vardict/) {
       next if ( $altd < 3 );
     }

     # generate reference panel for her2 brca
     if ($task =~ /refpanel/) {
       next if ($pass ne 'PASS');
       next if ($info =~ /SOMATIC/);
       next if (($info !~ /exonic/ and $info !~ /function\=UTR/ and $info !~ /splicing/) and $qual < 80);   #skip non exonic ones
       unless ($task =~ /huzheng/) {
         next if ($info =~ /dbSNP/ or $info =~ /1KG\=/ or $info =~ /ESP\d+\=/);  #skip known ones
       }
     }
     # generate reference panel for her2 brca

     #for titan purpose
     my $GTblood;
     my $DPblood;
     if ($task =~ /titan/) {
       if ($blood ne '') {
         my @bloodInfo = split(/\:/, $blood);             #for titan, get blood info
         $GTblood = $bloodInfo[$formindex{'GT'}];
         $DPblood = $bloodInfo[$formindex{'DP'}];
         if ( $GTblood !~ /1/ or $DPblood == 0 ) {     #skip vars that are not found in blood
           next;
         }
       } else {
         my @sampleInfo = split(":", $sample);
         $GTblood = $sampleInfo[$formindex{'GT'}];
         $DPblood = $sampleInfo[$formindex{'DP'}];
         if ($task =~ /guessNormal/) {                     # I made some changes for no normal data from allesio AFX/UPS
           if ($popFreqMax != -1 and $popFreqMax < 0.5) {
             $GTblood = '0/1';
           } else {
             if ($maf > 0.05){
               $GTblood = '0/1';
             } else {
               next;
             }
           }
         } else {
           next if ( $maf < $vafTitan or $DPblood == 0 );           #retain only with high frequency
         }
       }
     }
     #for titan purpose

     #functional filter
     my $function;
     $info =~ /(function=.+?$)/;
     $function = $1;
     if ($nonsegdup) {
       unless ($somatic == 1) {
         next if $function =~ /segdup\.score/;
       }
     }
     if ($exonic) {
       next if ($function !~ /exonic/ and $function !~ /UTR[35]/ and $function !~ /splicing/);
     }

     my $muTectLod = '';
     if ($task =~ /muTect/ and $info =~ /t\_lod\_fstar\=(.+?)\;init\_n\_lod\=(.+?)\;/){
       $muTectLod = $1.','.$2;
     }

     #now start to store data
     $somatic{$coor}{$name} = ($task =~ /muTect/)? $maf.'|'.$tdp.'|'.$muTectLod : $maf.'|'.$qual;
     if ($task =~ /strelka/) {
       $somatic{$coor}{$name} .= '|'.$altd.','.$tdp.'|'.$somatic;
     }
     if ($task =~ /titan/) {
       $somatic{$coor}{$name} .= '|'.$GTblood;
     }

     my $idrefalt = join("\t", ($id,$ref,$alt));
     $somatic{$coor}{'info'}{$idrefalt} = $function;  #not necessarily just one!

     if ($somatic == 1) {
        $somatic{$coor}{'somatic'} .= $name.',';
     } else {
        $somatic{$coor}{'germline'} .= $name.',';
     }
  }
  close IN;
  my $tmpc = scalar(keys %somatic);
  print STDERR "table tmp: $tmpc\n";

  my $cmd = "rm $tmpdir/tmp -f";
  RunCommand($cmd, 0 ,0);
}

#######------------------print header-----------------------#####
if ($task =~ /refpanel/){
  print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
  goto CONTENT;
}

print "#chr\tpos\tid\tref\talt";
if ($prefixReg ne ''){
#  foreach my $name (sort {$a =~ /($prefixReg)(\d+)?([a-zA-Z0-9\-\_]+)?/; my $pa = $1; my $ia = $2; my $ias = $3; $b =~ /($prefixReg)(\d+)?([a-zA-Z0-9\-\_]+)?/; my $pb = $1; my $ib = $2; my $ibs = $3; $pa cmp $pb or $ia <=> $ib or $ias cmp $ibs} keys %samples) {
  foreach my $name ( sort {$a cmp $b} keys %samples) {      #update rank
    print "\t$name";
  }
} elsif ($task =~ /tcga/i) {
  foreach my $name (sort {$a =~ /TCGA\-([A-Z0-9]+)\-([A-Z0-9]+)/; my $tsa = $1; my $inda = $2; $b =~ /TCGA\-([^\-]+)\-([^\-]+)/; my $tsb = $1; my $indb = $2; $tsa cmp $tsb or $inda cmp $indb} keys %samples) {
    print "\t$name";
  }
}
print "\tfunction\ttrace";
print "\n";
#################################################################

CONTENT:

foreach my $coor (sort {$a =~ /^([A-Za-z0-9\-\_\.]+)\:(\d+)$/; my $ca = $1; my $pa = $2; $b =~ /^([A-Za-z0-9\-\_\.]+)\:(\d+)$/; my $cb = $1; my $pb = $2; $ca cmp $cb or $pa <=> $pb} keys %somatic) {
  $coor =~ /^(\w+):(\d+)$/;
  my $chrom = $1;
  my $position = $2;
  my $pmaf = 0;
  my $pfreq = 0;
  my $aqual = 0;
  foreach my $info (keys (%{$somatic{$coor}{'info'}})) {
    print "$chrom\t$position\t$info";
    if ($prefixReg ne '') {
    #  foreach my $name (sort {$a =~ /($prefixReg)(\d+)?([a-zA-Z0-9\-\_]+)?/; my $pa = $1; my $ia = $2; my $ias = $3; $b =~ /($prefixReg)(\d+)?([a-zA-Z0-9\-\_]+)?/; my $pb = $1; my $ib = $2; my $ibs = $3; $pa cmp $pb or $ia <=> $ib or $ias cmp $ibs} keys %samples) {
      foreach my $name ( sort {$a cmp $b} keys %samples) {
        if ($somatic{$coor}{$name} ne '') {
          print "\t$somatic{$coor}{$name}";
        } else {
          print "\t0";
        }
      }
    } elsif ($task =~ /tcga/i) {
      my $totalSamps = 0;
      my $totalQualSamps = 0;
      my $totalMafs = 0;
      my $totalQuals = 0;
      foreach my $name (sort {$a =~ /TCGA\-([A-Z0-9]+)\-([A-Z0-9]+)/; my $tsa = $1; my $inda = $2; $b =~ /TCGA\-([^\-]+)\-([^\-]+)/; my $tsb = $1; my $indb = $2; $tsa cmp $tsb or $inda cmp $indb} keys %samples) {
        $totalSamps += 1;
        if ($task =~ /refpanel/) {
           my ($maf, $qual) = split(/\|/, $somatic{$coor}{$name});
           $totalMafs += $maf if ($maf > 0);
           $totalQuals += $qual if ($qual > 0);
           $totalQualSamps += 1 if ($qual > 0);
           next;
        }
        if ($somatic{$coor}{$name} ne '') {
          print "\t$somatic{$coor}{$name}";
        } else {
          print "\t0";
        }
      }
      $pmaf = sprintf("%g", $totalMafs/$totalSamps);
      $pfreq = $totalQualSamps.'/'.$totalSamps;
      $aqual = sprintf("%.1f", $totalQuals/$totalQualSamps);
    } #tcga
    my $function = $somatic{$coor}{'info'}{$info};
    my $somatic = ($somatic{$coor}{'somatic'} eq '')? 0 : $somatic{$coor}{'somatic'};       #temporarily silence somatic and germline info
    my $germline = ($somatic{$coor}{'germline'} eq '')? 0 : $somatic{$coor}{'germline'};    #temporarily silence somatic and germline info
    my $trace = 'somatic='.$somatic.';germline='.$germline;

    if ($task =~ /refpanel/) {
      $function .= 'pmaf='.$pmaf.';pfreq='.$pfreq;
      print "\t$aqual\tPASS\t$function\n";
      next;
    }
    print "\t$function\t$trace";
    print "\n";
  } #for each different variant in this coordinate
}


sub RunCommand {
  my ($command,$noexecute,$quiet) = @_ ;
  unless ($quiet){
    print STDERR "$command\n\n";
  }
  unless ($noexecute) {
    system($command);
  }
}
