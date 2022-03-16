use strict;
use Data::Dumper;
use File::Glob ':glob';
use List::Util qw[min max];
use Math::CDF qw[:all];
use Data::Dumper;
use Getopt::Long;
use Text::NSP::Measures::2D::Fisher::right;  #fisher test right tail

my $file;
my $task;
my $type;
my $somaticInfo;        # if for somatic judgement
my $bloodCall;          # whether the blood is also single called
my $rna;                # if for rna found
my $noStrandBias = "no";
my $loosefound = "no";
my $skipArFilter = "no";

my $Th_tumorLOD = 4.0;
my $Th_normalLOD = 2.3;
my $Th_maf = 0.03;
my $Th_endsratio = 0.9;
my $Th_vard = 2;
my $Th_badQualFrac = 0.6;    #0.6
my $Th_cmeancmedian = 5.5;   #5.5
my $Th_cmedian = 2;          #2
my $Th_bloodZero = 0.1;
if ($type eq 'indel') {
  $Th_endsratio = 0.95;
  $Th_cmeancmedian = 6.5;
  $Th_cmedian = 3;
}


GetOptions (
           "file|f=s"          => \$file,             #filename
           "type|t=s"          => \$type,             #snv or indel
           "task|k=s"          => \$task,             #task type
           "somaticInfo|s=s"   => \$somaticInfo,      #info for somatic sample pairs
           "bloodCall|b=s"     => \$bloodCall,        #bloodCall
           "cmeme=f"           => \$Th_cmeancmedian,  #cmeancmedian
           "cmedian=f"         => \$Th_cmedian,       #cmedian
           "rna|r=s"           => \$rna,              #rna info
           "noStrandBias=s"    => \$noStrandBias,     #treatment of strandbias
           "loosefound=s"      => \$loosefound,       #loose found samples
           "skipArFilter=s"      => \$skipArFilter,   #skip arbitrary filter
           "help|h"         => sub {
                               print "usage: $0 final preparation of somatic/germline variant calls \n\nOptions:\n\t--file\t\tthe filename of the res file\n";
                               print "\t--type\t\tthe type of variants, snv or indel\n";
                               print "\t--task\t\tthe aim of the analysis, maf, somatic, filter, (sam)founds, (rna)trace, depth, freq, etc\n";
                               print "\t--somaticInfo\tthe paired sample information\n";
                               print "\t--cmeme\tthe sum of concensus mismatch count (median and mean), default 5.5\n";
                               print "\t--cmedian\tthe concensus mismatch count (median only), default 2\n";
                               print "\t--bloodCall\twhether the blood has called for variants (yes or no), for a sanity check\n";
                               print "\t--rna\t\tthe dna-rna matching table, when task is rnatrace\n";
                               print "\t--help\t\tprint this help message\n";
                               print "\n";
                               exit 0;
                             },
           );


print STDERR "Th_maf: $Th_maf\n";
print STDERR "Th_endsratio: $Th_endsratio\n";
print STDERR "Th_vard: $Th_vard\n";
print STDERR "Th_badQualFrac: $Th_badQualFrac\n";
print STDERR "Th_cemancmedian: $Th_cmeancmedian\n";
print STDERR "Th_cmedian: $Th_cmedian\n";

my %VarClassMap;
$VarClassMap{'downstream'} = 'downstream';
$VarClassMap{'frameshift deletion'} = 'Frame_Shift_Del';
$VarClassMap{'frameshift insertion'} = 'Frame_Shift_Ins';
$VarClassMap{'intergenic'} = 'IGR';
$VarClassMap{'intronic'} = 'Intron';
$VarClassMap{'ncRNA_exonic'} = 'Non-coding_Transcript';
$VarClassMap{'ncRNA_intronic'} = 'Non-coding_Transcript';
$VarClassMap{'ncRNA_splicing'} = 'Non-coding_Transcript';
$VarClassMap{'nonframeshift deletion'} = 'In_frame_Del';
$VarClassMap{'nonframeshift insertion'} = 'In_frame_Ins';
$VarClassMap{'nonsynonymous SNV'} = 'Missense';
$VarClassMap{'splicing'} = 'Splice_site';
$VarClassMap{'stopgain'} = 'Nonsense';
$VarClassMap{'stoploss'} = 'Nonstop_Mutation';
$VarClassMap{'synonymous SNV'} = 'Synonymous';
$VarClassMap{'unknown'} = 'Silent';
$VarClassMap{'upstream'} = 'upstream';
$VarClassMap{'UTR3'} = "3\'UTR";
$VarClassMap{'UTR5'} = "5\'UTR";


my @all;
my %all;
my %somatic;
my %germline;  #may have multiple tumors
if ($somaticInfo and -s "$somaticInfo") {

  open IN, "$somaticInfo";
  while ( <IN> ) {
    chomp;
    s/[\s\n]$//;
    my @columns = split /\t/;
    my $tumor = $columns[0];
    my $normal = $columns[1];

    $somatic{$tumor} = $normal;
    push(@{$germline{$normal}}, $tumor) if $normal ne 'undef';
    if (!exists($all{$tumor})){
      push(@all, $tumor);
      $all{$tumor} = '';
    }
    if (!exists($all{$normal}) and $normal ne 'undef') {
      push(@all, $normal);
      $all{$normal} = '';
    }
  }
  close IN;
  #print STDERR Dumper (\%somatic);
  #print STDERR Dumper (\%germline);

}


my %rnasamps;
if ($rna) {
  if ( (! -e $rna) and $rna =~ /,/ ) {
    my @rnas = split(',', $rna);
    foreach my $rnas (@rnas) {
      $rnasamps{$rnas} = '';
    }
  } elsif (-e $rna) {
    my %rcol;
    open IN, "$rna";
    while ( <IN> ) {
      chomp;
      my @cols = split /\t/;
      if ($_ =~ /^([\#])?dna\t/) { #header
        for (my $i=0; $i<=$#cols; $i++) {
          $rcol{$cols[$i]} = $i;
        }
      } else {
        next if $cols[$rcol{'rna'}] eq 'NA';
        $rnasamps{$cols[$rcol{'rna'}]} = $cols[$rcol{'dna'}];
      }
    }
    close IN;
  }
  print STDERR Dumper(\%rnasamps);
}


open IN, "$file";
my %colnames;
my %colindex;
my $header;
while ( <IN> ) {
  chomp;
  if (/^[\#]?chr\t/) {
    #it is header
    $_ =~ s/^#//;
    my @cols = split /\t/;
    for(my $i = 0; $i <= $#cols; $i++) {
      $colnames{$cols[$i]} = $i;
      $colindex{$i} = $cols[$i];
    }
    if ($task eq 'maf') {
      print "$_\tmaf\n";
    } elsif ($task eq 'split') {
      #print "$_\n";
      $header = $_;
    } elsif ($task eq 'freq') {
      print "$_\tfreq\n";
    } elsif ($task eq 'trace') {
      print "$_\ttrace\n";
    } elsif ($task eq 'rnatrace') {
      print "$_\trnatrace\n";
    } elsif ($task eq 'founds') {
      print "$_\tfounds\n";
    } elsif ($task eq 'rnafounds') {
      print "$_\trnafounds\n";
    } elsif ($task =~ /filter/) {
      print "$_\tfilter\n";
    } elsif ($task =~ /somatic/) {
      print "$_\tsomatic\tgermline\n";
    } elsif ($task eq 'depth'){
      print "$_\tdepthav\n";
    } elsif ($task eq 'samfounds'){
      print "$_\tfounds\trsam\n";
    } elsif ($task eq 'mut2maf') {
      printf("%s\n", join("\t", "Hugo_Symbol","Chromosome","Start_Position","End_position","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele2","Tumor_Sample_Barcode"));
    } else {
      print STDERR "task is wierd\n";
      exit 22;
    }
  } else {
    my @cols = split /\t/;
    if ($task eq 'maf') {    #get real maf
      my $smaf = 0;
      my $sampleCounts = scalar(@all);
      foreach my $sample (@all) {
        my $sampmaf = $sample.'maf';
        my $sampd = $sample.'d';
        my $maf = $cols[$colnames{$sampmaf}];
        my $endsratio = 0;
        my $cmean = 0;
        my $cmedian = 0;
        my $strandRatio = 0;
        my $strandRatioRef = 0;
        my $strandFisherP = 0;
        my $badQualFrac = 0;
        my $depth = $cols[$colnames{$sampd}];
        if ($cols[$colnames{$sampmaf}] =~ /\|/) {
          my @infos = split(/\|/, $cols[$colnames{$sampmaf}]);
          $maf = $infos[0];
          $endsratio = $infos[1];
          ($cmean, $cmedian) = split(',', $infos[2]);
          ($strandRatio, $strandRatioRef, $strandFisherP) = split(',', $infos[3]);
          $badQualFrac = $infos[4];
        }
        my $vard = sprintf("%.1f", $maf*$depth);
        my $refd = $depth-$vard;
        unless (($endsratio < 0.9 or ((1-$endsratio)*$vard >= 2)) and $badQualFrac <= 0.7 and (($strandRatio > 0 and $strandRatio < 1) or ($noStrandBias eq 'no' and $strandFisherP > 0.7 and $refd >= 10 and $vard >= 5 and $maf >= 0.1)) and (($cmean+$cmedian) < $Th_cmeancmedian or $cmedian <= $Th_cmedian)) {
          $maf = 0;
        }
        if ($maf >= 0.1 and $vard >= 2) {
            $smaf += $maf;
        }
      }
      $smaf = sprintf("%.6f",$smaf/$sampleCounts);
      print "$_\t$smaf\n";
    } elsif ($task eq 'split' or $task eq 'mut2maf') {

      (my $Hugo_Symbol = $cols[$colnames{'geneName'}]) =~ s/\(dist=\d+\)//g;
      my $Chromosome = $cols[$colnames{'chr'}];
      my $Start_Position = $cols[$colnames{'pos'}];
      my $End_Position = $cols[$colnames{'pos'}];
      my $geneLoc = $cols[$colnames{'geneLoc'}];
      my $Variant_Classification = ($cols[$colnames{'functionalClass'}] eq 'NA')? $cols[$colnames{'geneLoc'}]:$cols[$colnames{'functionalClass'}];
      $Variant_Classification = $VarClassMap{$Variant_Classification};
      my $Variant_Type = "SNP";
      my $Reference_Allele = $cols[$colnames{'ref'}];
      my $Tumor_Seq_Allele2 = $cols[$colnames{'alt'}];
      my $Tumor_Sample_Barcode = "";

      my @spliceIndexes;
      foreach my $sample (@all) {
        push(@spliceIndexes,$colnames{$sample}) if ( exists($colnames{$sample}) );
        my $sampmaf = $sample.'maf';
        my $sampd = $sample.'d';
        my $maf = $cols[$colnames{$sampmaf}];
        if (!exists($colnames{$sampmaf})) {   #sample realmaf not present, skip
          next;
        }
        my $endsratio = 0;
        my $cmean = 0;
        my $cmedian = 0;
        my $strandRatio = 0;
        my $strandRatioRef = 0;
        my $strandFisherP = 0;
        my $badQualFrac = 0;
        my $depth = $cols[$colnames{$sampd}];
        if ($cols[$colnames{$sampmaf}] =~ /\|/) {
          my @infos = split(/\|/, $cols[$colnames{$sampmaf}]);
          $maf = $infos[0];
          $endsratio = $infos[1];
          ($cmean, $cmedian) = split(',', $infos[2]);
          ($strandRatio, $strandRatioRef, $strandFisherP) = split(',', $infos[3]);
          $badQualFrac = $infos[4];
        }
        my $vard = sprintf("%.1f", $maf*$depth);
        my $refd = $depth-$vard;
        if (exists($somatic{$sample})) { #do this only for tumorsample
          unless (($endsratio < 0.9 or ((1-$endsratio)*$vard >= 2)) and $badQualFrac <= 0.7 and (($strandRatio > 0 and $strandRatio < 1) or ($noStrandBias eq 'no' and $strandFisherP > 0.7 and $refd >= 10 and $vard >= 5 and $maf >= 0.1)) and (($cmean+$cmedian) < $Th_cmeancmedian or $cmedian <= $Th_cmedian)) {
            $maf = 0;
          }
        }
        $cols[$colnames{$sampmaf}] = $maf;

        if ($task eq 'mut2maf') {  #turning mut to maf format
          $Tumor_Sample_Barcode = $sample;
          next if exists($germline{$sample});
          if ($maf > 0 and $cols[$colnames{'germline'}] eq 'NA' and (($cols[$colnames{'rep'}] == 0 and $cols[$colnames{'rep'}] == 0) or $cols[$colnames{'cmedianav'}] < 2) ) {
            printf("%s\n",join("\t",$Hugo_Symbol,$Chromosome,$Start_Position,$End_Position,$Variant_Classification,$Variant_Type,$Reference_Allele,$Tumor_Seq_Allele2,$Tumor_Sample_Barcode));
          }
        }
      }
      if ($task eq 'split') {
        @cols = &splice_entry(\@cols, \@spliceIndexes) if ($#spliceIndexes > 0);
        if ($header ne '') {    #print header cut
          my @headercols = split("\t", $header);
          @headercols = &splice_entry(\@headercols, \@spliceIndexes) if ($#spliceIndexes > 0);
          printf("%s\n", join("\t", @headercols));
          $header = '';
        }
        printf("%s\n", join("\t", @cols));
      }
    } elsif ($task eq 'freq') {
      my $freq = 'NA';
      #my $freq1KG;
      #my $freqESP;
      my $function = (exists($colnames{'function'}))? $cols[$colnames{'function'}] : die("no function column.\n");
      #my $id = $cols[$colnames{'id'}];
      if ($function =~ /PopFreqMax=(.+?)\;/){
        $freq = $1;
      } elsif ($function =~ /ExAC_ALL=(.+?)\;/){
        $freq = $1;
      } elsif ($function =~ /ESP\d+=(.+?)\;/){
        $freq = $1;
      } elsif ($function =~ /1000g2014oct_all=(.+?)\;/){
        $freq = $1;
      } elsif ($function =~ /1KG=(.+?)\;/){
        $freq = $1;
      }
      print "$_\t$freq\n";
    } elsif ($task =~ /founds/ or $task =~ /trace/) {    #trace all samples and check whether it is originally called
      my $founds = 0;
      my $samfounds = 0;
      my $trace = '';
      for (my $i = 0; $i <= $#cols; $i++) {
        if ($colindex{$i} =~ /^(.+?)maf$/) {
          my $samp = $1;
          if ($task =~ /rna/) {      #for rna
            next if (!exists($rnasamps{$samp}));
          }
          my $maf = $cols[$i];
          my $endsratio = 0;
          my $cmean = 0;
          my $cmedian = 0;
          my $strandRatio = 0;
          my $strandRatioRef = 0;
          my $strandFisherP = 0;
          my $badQualFrac = 0;

          if ($cols[$i] =~ /\|/) { #split the var surrounding information
            my @infos = split(/\|/, $cols[$i]);
            $maf = $infos[0];
            $endsratio = $infos[1];
            ($cmean, $cmedian) = split(',', $infos[2]);
            ($strandRatio, $strandRatioRef, $strandFisherP) = split(',', $infos[3]);
            $badQualFrac = $infos[4];
          }

          my $depth = $cols[$i+1];
          if ($depth =~ /\,/){
            my @depths = split(',', $cols[$i+1]);   #spandepth, jumpdepth
            $depth = $depths[0];
          }
          my $vard = sprintf("%.1f", $maf*$depth);
          my $refd = $depth-$vard;

          my $SatisfyCondition = 0;
          if ($loosefound eq 'no') {
            $SatisfyCondition = (($endsratio <= $Th_endsratio or ((1-$endsratio)*$vard >= $Th_vard)) and $badQualFrac <= $Th_badQualFrac and (($strandRatio > 0 and $strandRatio < 1) or ($noStrandBias eq 'no' and $strandFisherP > 0.7 and $refd >= 10 and $vard >= 5 and $maf >= 0.1)) and (($cmean+$cmedian) < $Th_cmeancmedian or $cmedian <= $Th_cmedian))? 1 : 0;
          } else {    #loosen criteria
            $SatisfyCondition = ($vard > 5 and ($cmean+$cmedian) <= 7 and $cmedian <= 3)? 1 : 0;
          }
          if ($SatisfyCondition == 1) {  #it looks good
            if ($task =~ /rna/) {
              if ($maf >= ($Th_maf - 0.02) and $vard >= $Th_vard) {
                $founds++;
                $trace .= "$samp,";
              }
            } else {
              if ($maf >= $Th_maf and $vard >= ($Th_vard+1)) {
                if ($somaticInfo ne '') { #count only tumor
                  if ( exists($somatic{$samp}) ) {
                    $founds++;
                    $trace .= "$samp,";
                    if ($colindex{$i-1} =~ /^$samp$/) { #calling information for the same sample
                      if ( $cols[$i-1] =~ /\|/ ) { #it is called originally
                        $samfounds++;
                      }
                    }
                  }
                } else {
                  $founds++;                       #count all samples
                  $trace .= "$samp,";
                  if ($colindex{$i-1} =~ /^$samp$/) { #calling information for the same sample
                    if ( $cols[$i-1] =~ /\|/ ) { #it is called originally
                      $samfounds++;
                    }
                  }
                }
              }                 #maf and vard requirements
            } #for dna
          }  # it looks good
        } #maf
      } #each column
      if ($task =~ /founds/) {
        print "$_\t$founds";
        if ($task =~ /sam/) {    #original founds
          my $rsam = ($founds > 0)? sprintf("%.2f", $samfounds/$founds) : 0;
          print "\t$rsam";
        }
        print "\n";
      } elsif ($task =~ /trace/) {
        if ($trace ne ''){
          print "$_\t$trace\n";
        } else {
          print "$_\tNA\n";
        }
      }
    } elsif ($task eq 'depth') {   #find av depth
      my $dep = 0;
      my $Ndep = 0;
      for (my $i = 0; $i <= $#cols; $i++) {
        if ($colindex{$i} =~ /maf$/) {
          if ($cols[$i] >= 0.1){    #found clonal
            my $vard = sprintf("%.1f", $cols[$i]*$cols[$i+1]);
            if ($vard >= $Th_vard) {
              $Ndep++;
              $dep += $cols[$i+1];
            }
          }
        } #maf
      } #each column
      my $depav = ($Ndep > 0)? sprintf("%.1f", $dep/$Ndep):0;
      print "$_\t$depav\n";
    } elsif ($task =~ /filter/) {     #filter
      my @detectedSample;
      my %detectedSample;
      my $somaticCalled = 0;
      my %mafs;
      my $chr;
      my $pos;
      my $endsratio = 0;
      my @strandRatio;
      my $strandRatio = 0;
      my $strandRatioRef = 0;
      my $strandFisherP = 0;
      my $badQualFrac = 0;
      my $cmean = 0;
      my $cmedian = 0;
      my $cmeanav = 0;
      my $cmedianav = 0;
      my $mmaf = 0;
      my $mlod = 0;
      my $rep = 0;
      my $sc = 0;
      for ( my $i = 0; $i <= $#cols; $i++ ) {
        if ($colindex{$i} eq 'chr') {
          $chr = $cols[$i];
        } elsif ($colindex{$i} eq 'pos') {
          $pos = $cols[$i];
        } elsif ($colindex{$i} eq 'trace') {
          my ($traceSomatic, $traceGermline) = split(';', $cols[$i]);
          $traceSomatic =~ /somatic\=(.+?)$/;
          $traceSomatic = $1;
          $traceGermline =~ /germline\=(.+?)$/;
          $traceGermline = $1;
          if ($traceSomatic =~ /\,/) {
            my @traceSomatic = split(/\,/, $traceSomatic);
            push(@detectedSample, @traceSomatic);
            $somaticCalled = 1;
          }
          if ($traceGermline =~ /\,/) {
            my @traceGermline = split(/\,/, $traceGermline);
            push(@detectedSample, @traceGermline);
          }
          foreach my $detectedSamp (@detectedSample) {
            $detectedSample{$detectedSamp} = '';
          }
        } elsif ($colindex{$i} eq 'rep') {
          $rep = $cols[$i];
        } elsif ($colindex{$i} eq 'sc') {
          $sc = $cols[$i];
        } elsif ($colindex{$i} =~ /^(.+?)maf$/) {   #store maf information
          my $samp = $1;
          $mafs{$samp} = $cols[$i].'|'.$cols[$i+1];  #maf + depth
        } elsif ($colindex{$i} eq 'cmeanav') {
          $cmeanav = $cols[$i];
          $cmedianav = $cols[$i+1];
        }
      } #each column

      foreach my $samp (keys %mafs) {  #get endsratio and cmean cmedian info
        next if (! exists($detectedSample{$samp}));   #only look at the sample where it is called
        my $sampmaf = $mafs{$samp};
        if ($sampmaf =~ /\|/) {        #split the var surrounding information
          my @infos = split(/\|/, $sampmaf);
          $endsratio = ($infos[5] > $mlod)? $infos[1]:$endsratio;
          ($cmean, $cmedian) = ($infos[5] > $mlod)? split(',', $infos[2]):($cmean, $cmedian);
          @strandRatio = ($infos[5] > $mlod)? split(',', $infos[3]):@strandRatio;
          $strandRatio = $strandRatio[0];
          $strandRatioRef = $strandRatio[1];
          $strandFisherP = $strandRatio[2];
          $badQualFrac = ($infos[5] > $mlod)? $infos[4]:$badQualFrac;
          #$mmaf = ($infos[5] > $mlod)? $infos[0]:$mmaf;
          $mlod = ($infos[5] > $mlod)? $infos[5]:$mlod;
        }
      }

      my $status;
      if ($rep == 1 and $sc == 1) {
        $status = ($endsratio < $Th_endsratio and $badQualFrac <= ($Th_badQualFrac-0.1) and (($cmean+$cmedian) < ($Th_cmeancmedian-1) or $cmedian < $Th_cmedian) and ($cmeanav + $cmedianav) <= ($Th_cmeancmedian-0.3))? 'PASS':'FOUT';   #conservative for rep and sc
      } elsif ($rep == 1 or $sc == 1) {
        $status = ($endsratio < $Th_endsratio and $badQualFrac <= $Th_badQualFrac and (($cmean+$cmedian) <= $Th_cmeancmedian or $cmedian <= $Th_cmedian) and ($cmeanav + $cmedianav) <= $Th_cmeancmedian)? 'PASS':'FOUT';
      } else {
        $status = ($endsratio < $Th_endsratio and $badQualFrac <= $Th_badQualFrac and (($cmean+$cmedian) <= $Th_cmeancmedian or $cmedian <= $Th_cmedian) and ($cmeanav + $cmedianav) <= $Th_cmeancmedian)? 'PASS':'FOUT';
      }

      if ($type eq 'indel' and $somaticCalled == 1){
        $status = 'PASS';
      }
      if ($task !~ /keep/) {
        print "$_\t$status\n" if ($status eq 'PASS');
        #print STDERR "$_\t$status\n" if ($status eq 'FOUT');
      } else {
        print "$_\t$status\n";
      }
    } elsif ($task =~ /somatic/) {  #find somatic ones
      my %tumor;
      my %blood;
      my %bloodcalled;   #store the called information for blood
      my %nonblood;
      my %unknown;
      for ( my $i = 0; $i <= $#cols; $i++ ) {

        if ($colindex{$i} =~ /^(.+?)maf$/) {
          my $samp = $1;
          my $maf = $cols[$i];
          my $endsratio = 0;
          my @strandRatio;
          my $strandRatio = 0;
          my $strandRatioRef = 0;
          my $strandFisherP = 0;
          my $badQualFrac = 0;
          my $cmean = 0;
          my $cmedian = 0;
          my $lod = 0;

          if ($cols[$i] =~ /\|/) {  #split the var surrounding information
            my @infos = split(/\|/, $cols[$i]);
            $maf = $infos[0];
            $endsratio = $infos[1];
            ($cmean, $cmedian) = split(',', $infos[2]);
            @strandRatio = split(',', $infos[3]);
            $strandRatio = $strandRatio[0];
            $strandRatioRef = $strandRatio[1];
            $strandFisherP = $strandRatio[2];
            $badQualFrac = $infos[4];
            $lod = $infos[5];
          }

          my $depth = $cols[$i+1];
          #my $vard = sprintf("%.1f", $maf*$depth);
          my $vard = round($maf*$depth);
          my $refd = $depth-$vard;

          #DECIDE MAF
          if (exists $somatic{$samp} and $skipArFilter eq 'no') {     #for tumor samples require some additional thing, if the arbitrary filter is turned off, the original MAF will be retained.
            if (($endsratio <= $Th_endsratio or ((1-$endsratio)*$vard >= $Th_vard)) and $badQualFrac <= $Th_badQualFrac and (($strandRatio > 0 and $strandRatio < 1) or ($noStrandBias eq 'no' and $strandFisherP > 0.7 and $refd >= 10 and $vard >= 5 and $maf >= 0.1)) and (($cmean+$cmedian) < ($Th_cmeancmedian-0.3) or $cmedian <= $Th_cmedian)) {    #simple filter to decide true event
              if ( $maf >= 0.1 ) {  #clonal ones
                $maf = $maf;
              } else {              #subclonal ones, subject to additional constrains
                if ($badQualFrac < ($Th_badQualFrac-0.2) and $cmedian <= ($Th_cmedian-0.3)) {
                  $maf = $maf;
                } else {
                  $maf = 0;         #not reliable somatic
                }
              }
            } else {
              $maf = 0;         #not reliable somatic
            }
          }

          #save info
          if (exists $somatic{$samp}) {       #it is tumor sample name
            my $tumorLOD = $lod;
            if ($vard >= ($Th_vard+1) and $maf >= ($Th_maf+0.01) and $depth >= 8) {
              if ($vard >= 50) {              #must be sig anyway
                $tumor{$samp} = join(',', $maf, $vard, $depth);               #maf - > maf, vard, depth
              } else {                        #cal for less vard
                #print STDERR "$samp\t$maf\t$vard\t$depth\t$tumorLOD\n";
                if ($tumorLOD >= $Th_tumorLOD) {
                  $tumor{$samp} = join(',', $maf, $vard, $depth);             #maf - > maf, vard, depth
                }
              }
            }
          } elsif (exists $germline{$samp}) {    #it is blood
            my $normalLOD = $lod;
            foreach my $ct (@{$germline{$samp}}) {
              if ($bloodCall eq 'yes' and $cols[$i-1] =~ /\|/) {                         #it is originally called
                $bloodcalled{$ct} = '';
              }
              #print STDERR "$samp\t$maf\t$vard\t$depth\t$normalLOD\n";
              if ( $maf == 0 ) {   #no blood alt found
                if ( ($normalLOD > $Th_normalLOD or $cols[$i] == 0) and $depth >= 8 ) {  #good normalLOD or no variant reads
                  $nonblood{$ct} = '';
                } else {
                  $unknown{$ct} = '';
                }
              } else {          #maf != 0, should keep normalLOD in mind
                $blood{$ct} = join(',', $maf, $vard, $depth, $normalLOD);              #maf,normalLOD - > maf, vard, depth,normalLOD
              }
            }
          }                     #it is blood
        }                       #maf
      }                         #each column

      my $soma = 'NA';
      my $germ = 'NA';
      foreach my $tumorSamp (keys %tumor) {   #check blood to confirm that it is somatic

        my ($tumorMaf, $tumorVard, $tumorDepth) = split(',', $tumor{$tumorSamp});
        my $stype = 'NA';
        if (exists $nonblood{$tumorSamp}) {
          $stype = 'good';
          $stype .= ($tumorMaf < 0.1)? 'Sub' : '';  #add subclonal info
          $soma = ($soma eq 'NA')? $tumorSamp."\[$stype\]".',':$soma.$tumorSamp."\[$stype\]".',';

        } elsif (exists($blood{$tumorSamp})) {

          my ($bloodMaf, $bloodVard, $bloodDepth, $bloodLOD) = split(',',$blood{$tumorSamp});
          my $pFisher = calculateStatistic(n11=>$tumorVard, n1p=>$tumorDepth, np1=>($tumorVard+$bloodVard), npp=>($tumorDepth+$bloodDepth));    #fisher exact test to see whether tumor enrich for variants
          #print STDERR "$tumorVard\t$tumorDepth\t$bloodVard\t$bloodDepth\t$pFisher\n";

          if ( $bloodMaf < $Th_bloodZero and $tumorMaf/$bloodMaf >= 4 and $bloodLOD > $Th_normalLOD and $pFisher < 0.05 ) {
            $stype = 'doubt';
            $stype .= ($tumorMaf < 0.1)? 'Sub' : '';  #add subclonal info
            $soma = ($soma eq 'NA')? $tumorSamp."\[$stype\]".',':$soma.$tumorSamp."\[$stype\]".',';
          } elsif (exists($unknown{$tumorSamp})) {
            $stype = 'undef';
            $stype .= ($tumorMaf < 0.1)? 'Sub' : '';  #add subclonal info
            $soma = ($soma eq 'NA')? $tumorSamp."\[$stype\]".',':$soma.$tumorSamp."\[$stype\]".',';
          } else {
            $germ = ($germ eq 'NA')? $tumorSamp.',':$germ.$tumorSamp.',';
          }

        } else { #unknown ones, all germ
          $germ = ($germ eq 'NA')? $tumorSamp.',':$germ.$tumorSamp.',';
        }
      }

      if ($soma eq 'NA' and $germ eq 'NA') {     #check whether it is germline
        foreach my $ct (sort keys %bloodcalled){
          $germ .= $ct.',' if ($ct ne '');
        }
      }
      print "$_\t$soma\t$germ\n" if ($soma ne 'NA' or $germ ne 'NA');
    } #somatic
  }
}
close IN;

exit 0;

sub round {
    my $number = shift;
    my $tmp = int($number);
    if ($number >= ($tmp+0.5)){
      $tmp++;
    }
    return $tmp;
}

sub splice_entry {
  my $cols = shift;
  my @newcols = @{$cols};
  my $spliceIndexes = shift;
  my $decr = 0;
  for (sort {$a <=> $b} @{$spliceIndexes}) {
    splice(@newcols, $_+$decr, 1);
    --$decr;
  }
  return(@newcols);
}
