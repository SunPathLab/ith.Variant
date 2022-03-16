use strict;
use File::Glob ':glob';
use List::Util qw[min max];
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use FindBin qw($RealBin);
use lib "$RealBin/lib";
use snvCalling;
use Text::NSP::Measures::2D::Fisher::right;


my $file;      #filename of all rechecked files
my $type;
my $original;
my $task;      #only for tcga
my $prefix;
my $blood;
my $readlen = 76;

GetOptions (
           "file|f=s"       => \$file,             #filename
           "type|t=s"       => \$type,             #snv or indel
           "original|o=s"   => \$original,         #original big table
           "task|k=s"       => \$task,             #task type
           "prefix|p=s"     => \$prefix,
           "readlen=i"      => \$readlen,
           "blood|b=s"      => \$blood,            #blood
           "help|h"         => sub{
                               print "usage: $0 get all minor allele frequency for samples under recheck\n\nOptions:\n\t--file\t\tthe filename of all rechecked files\n";
                               print "\t--type\t\tthe type of variants, snv or indel\n";
                               print "\t--original\tthe original mutation big table\n";
                               print "\t--prefix\tthe prefix of samples' names, comma separated\n";
                               print "\t--blood\t\tthe sample names of blood samples\n";
                               print "\t--task\t\tthe task, such as tcga, or errorEst (estimate global sequencing error rate)\n";
                               print "\t--help\t\tprint this help message\n";
                               print "\n";
                               exit 0;
                             },
           );



my %blood;
if ($blood ne ''){
  foreach my $bl (split(/\,/, $blood)) {
    $blood{$bl} = '';
  }
}

my @prefix = split(',', $prefix);
my $prefixReg = join('|', @prefix);
print STDERR "prefixReg is $prefixReg\n";

my @list;
open IN, "$file";
while ( <IN> ) {
  chomp;
  next if /^#/;
  push(@list, $_) unless ($_ eq '');
}
close IN;


####################################################get chr pos
my %chrJumper;
$chrJumper{'original'} = getchrpos($original);

my %samples;
foreach my $file (@list) {
  my $filebase = basename($file);
  my $name;
  if ($prefix ne '' and $filebase =~ /(($prefixReg)([A-Za-z0-9\-\_]+)?)/) {
    $name = $1;
  } elsif ($task eq 'tcga' and $filebase =~ /((TCGA\-([^\-]+\-[^\-]+))\-[^\-]+\-[^\-]+\-[^\-]+\-\d+)$/) {
    $name = $1;
  }
  if ($name ne '') {
    $samples{$name} = '';
    $chrJumper{$file} = getchrpos($file);   #remember the chr start for each file
  }
}

print STDERR Dumper(\%chrJumper);
####################################################get chr pos


print "#chr\tpos\tid\tref\talt";
#foreach my $name (sort {$a =~ /($prefixReg)(\d+)?([A-Za-z0-9\-\_]+)?/; my $pa = $1; my $ia = $2; my $ias = $3; $b =~ /($prefixReg)(\d+)?([A-Za-z0-9\-\_]+)?/; my $pb = $1; my $ib = $2; my $ibs = $3; $pa cmp $pb or $ia <=> $ib or $ias cmp $ibs} keys %samples) {
foreach my $name (sort {$a cmp $b} keys %samples) {     #avoid problem when names overlap, strict rule for names
   print "\t$name\t$name".'d';
}
print "\tcmeanav\tcmedianav\n";


my %errorEst;                    #for the estimation of error rate
foreach my $chrc (sort keys %{$chrJumper{'original'}}) {

  my $jumperO = $chrJumper{'original'}->{$chrc};
  open OR, "$original";
  seek(OR, $jumperO, 0);
  my %OR;
  my $diji = 1;
  my $oldCoor = 'SRP';
  while ( <OR> ) {
    next if /^[@#]/;
    chomp;
    my @cols = split /\t/;
    my $chr = $cols[0];

    last if $chr ne $chrc;

    my $pos = $cols[1];
    my $id  = $cols[2];
    my $coor = $chr.':'.$pos;
    if ($coor eq $oldCoor) {
      $diji += 1;
    } else {
      $diji = 1;
    }
    my $ref = $cols[3];
    my $alt = $cols[4];
    $OR{$coor}{$diji} = $id.','.$ref.','.$alt;

    #redefine
    $oldCoor = $coor;
  }
  close OR;
  print STDERR "$original $chrc loaded\n";


  my %somatic;
  foreach my $file (@list) {       #foreach file
    my $filebase = basename($file);
    my $name;
    if ($prefix ne '' and $filebase =~ /(($prefixReg)([A-Za-z0-9\-\_]+)?)/) {
      $name = $1;
    } elsif ($task eq 'tcga' and $filebase =~ /((TCGA\-([^\-]+\-[^\-]+))\-[^\-]+\-[^\-]+\-[^\-]+\-\d+)$/) {
      $name = $1;
    }

    print STDERR "$name\n";     #now the name is sample Name

    my $jumperI = $chrJumper{$file}->{$chrc};    #rechecked file jumper
    my $djindex = 1;
    my $prevCoor = 'SRP';
    open IN, "$file";
    seek(IN, $jumperI, 0);
    while ( <IN> ) {
      chomp;
      next if /^[#@]/;
      next if /^chr\t/;
      if ($type =~ /indel/) {       #indel
        my ($chr, $pos, $ref, $alt, $indelType, $depth, $vardp, $vardn, $vends, $junction, $badqual, $cmean, $cmedian) = split /\t/;
        my $vard = $vardp + $vardn;

        if ($cmean =~ /e/) {
          $cmean = 0;
        }
        last if $chr ne $chrc;

        my $coor = $chr.':'.$pos;
        if ($coor eq $prevCoor) {
          $djindex += 1;
        } else {
          $djindex = 1;
        }

        if ($vard > 0 and $depth > 0) {
          my $endratio = sprintf("%.4f", $vends/$vard);
          my $strandRatio = sprintf("%.4f", $vardp/$vard);
          $somatic{$coor}{$djindex}{$name} = sprintf("%.4f", $vard/$depth);
          $somatic{$coor}{$djindex}{$name} .= '|'.$endratio.'|'.$cmean.','.$cmedian.'|'.$strandRatio.'|'.$badqual;
        } else {
          $somatic{$coor}{$djindex}{$name} = 0;
        }
        $somatic{$coor}{$djindex}{$name} .= "\t$depth";
        if ($junction != 0) {  #there are some junction reads
          $somatic{$coor}{$djindex}{$name} .= ",$junction";
        }
        $somatic{$coor}{$djindex}{'info'} = join("\t", ($ref,$alt));
        $somatic{$coor}{$djindex}{'consecutive'} .= $cmean.','.$cmedian.','.$vard.';';

        #redefine
        $prevCoor = $coor;

      } elsif ($type =~ /snv/) {    #snv

        my ($chr, $pos, $depth, $pstrand, $nstrand, $F1R2all, $F2R1all, $F1R2alt, $F2R1alt, $vard, $A, $An, $C, $Cn, $G, $Gn, $T, $Tn, $vends, $junction, $badqual, $cmean, $cmedian, $indmean, $indmedian, $vrlen, $localEr, $phred) = split /\t/;

        if ($cmean =~ /e/) {
          $cmean = 0;
        }
        last if $chr ne $chrc;

        my $coor = $chr.':'.$pos;
        if ($coor eq $prevCoor) {   #for the same coordinate but different var
          $djindex += 1;
        } else {
          $djindex = 1;
        }

        if ( $somatic{$coor}{$djindex}{$name} ne '' and $somatic{$coor}{$djindex}{$name} !~ /^0\t/ ) {   #already detected with certain reads in a previous file, for merging purposes
          my @alreadyHave = split(/\t/, $somatic{$coor}{$djindex}{$name});
          my @alreadyDepths = split(',', $alreadyHave[1]);
          my $alreadyDepth = $alreadyDepths[0];
          if ($alreadyDepths[0] > $depth) {
            goto SEENB;
          }
        }

        ############################ here adding data ###########################################
        if ( $vard > 0 and $depth > 0 ) {
          if ( $OR{$coor}{$djindex} ne '' ) {
            my @information = split(",", $OR{$coor}{$djindex});
            my $ref = $information[1];
            my $alt = $information[2];
            my $altd;
            my $strandRatio;
            my $strandRatioRef;        #strand bias ratio from reference allele
            my $stranFisherP;          #strand bias fisher p value, calculate from $pstrand and $nstrad
            if ($alt eq 'A') {
              $altd = $A + $An;
              $strandRatio =($altd > 0)? sprintf("%.4f", $A/$altd) : 0;
              $strandRatioRef = (($depth-$altd) > 0)? sprintf("%.4f", ($pstrand-$A)/($depth-$altd)) : 0;
              $stranFisherP = ($altd > 0)? calculateStatistic(n11=>max($A, $An), n1p=>$altd, np1=>max($A, $An)+max(($pstrand-$A), ($nstrand-$An)), npp=>$depth) : 1;
            } elsif ($alt eq 'C') {
              $altd = $C + $Cn;
              $strandRatio =($altd > 0)? sprintf("%.4f", $C/$altd) : 0;
              $strandRatioRef = (($depth-$altd) > 0)? sprintf("%.4f", ($pstrand-$C)/($depth-$altd)) : 0;
              $stranFisherP = ($altd > 0)? calculateStatistic(n11=>max($C, $Cn), n1p=>$altd, np1=>max($C, $Cn)+max(($pstrand-$C), ($nstrand-$Cn)), npp=>$depth) : 1;
            } elsif ($alt eq 'G') {
              $altd = $G + $Gn;
              $strandRatio =($altd > 0)? sprintf("%.4f", $G/$altd) : 0;
              $strandRatioRef = (($depth-$altd) > 0)? sprintf("%.4f", ($pstrand-$G)/($depth-$altd)) : 0;
              $stranFisherP = ($altd > 0)? calculateStatistic(n11=>max($G, $Gn), n1p=>$altd, np1=>max($G, $Gn)+max(($pstrand-$G), ($nstrand-$Gn)), npp=>$depth) : 1;
            } elsif ($alt eq 'T') {
              $altd = $T + $Tn;
              $strandRatio =($altd > 0)? sprintf("%.4f", $T/$altd) : 0;
              $strandRatioRef = (($depth-$altd) > 0)? sprintf("%.4f", ($pstrand-$T)/($depth-$altd)) : 0;
              $stranFisherP = ($altd > 0)? calculateStatistic(n11=>max($T, $Tn), n1p=>$altd, np1=>max($T, $Tn)+max(($pstrand-$T), ($nstrand-$Tn)), npp=>$depth) : 1;
            } else {
              print STDERR "$coor\t$alt\talt is not ACGT\n";
              exit 22;
            }
            $stranFisherP = sprintf("%.5f", $stranFisherP);   #keep precision 5
            if ($altd > 0) {

              #prepare ToxoG
              my $ToxoG = join(',', $F1R2all, $F2R1all, $F1R2alt, $F2R1alt);

              #prepare LODs
              if ( $depth < 5000/$readlen ) {
                $localEr = 0.001;
              }

              #for some wierd bug when doing RNA-seq ###need to check latter!!!!
              if ($altd > $depth){
                $altd = $depth;
              }
              #for some wierd bug when doing RNA-seq ###need to check latter!!!!


              if (exists($blood{$name}) or $blood eq 'yes') { #it is blood
                my $endratio = sprintf("%.4f", $vends/$vard);
                $somatic{$coor}{$djindex}{$name} = sprintf("%.4f", $altd/$depth);
                ########################### NLOD ###########################
                my $nlod = snvCalling->calNormalLOD($phred, $localEr, 0.001, $somatic{$coor}{$djindex}{$name}, $altd, $depth);
                ############################################################
                $somatic{$coor}{$djindex}{$name} .= '|'.$endratio.'|'.$cmean.','.$cmedian.'|'.$strandRatio.','.$strandRatioRef.','.$stranFisherP.'|'.$badqual.'|'.$nlod.'|'.$ToxoG.'|'.$indmean.','.$indmedian.'|'.$vrlen;
              } else {  #it is tumor
                my $endratio = sprintf("%.4f", $vends/$vard);
                if (($endratio <= 0.8 or ($altd - $vends) >= 2) and (($cmean+$cmedian) < 6 or $cmedian <= 2)) {  #limiting endsratio and mismatch stuff
                  $somatic{$coor}{$djindex}{$name} = sprintf("%.4f", $altd/$depth);

                  #print STDERR "$coor\t$phred\t$localEr\t0.001\t$somatic{$coor}{$djindex}{$name}\t$altd\t$depth\n";    #for debugging

                  ########################### TLOD ###########################
                  my $tlod = snvCalling->calTumorLOD($phred, $localEr, 0.001, $somatic{$coor}{$djindex}{$name}, $altd, $depth);
                  ############################################################
                  $somatic{$coor}{$djindex}{$name} .= '|'.$endratio.'|'.$cmean.','.$cmedian.'|'.$strandRatio.','.$strandRatioRef.','.$stranFisherP.'|'.$badqual.'|'.$tlod.'|'.$ToxoG.'|'.$indmean.','.$indmedian.'|'.$vrlen;
                } else {  #looks like artifact
                  $somatic{$coor}{$djindex}{$name} = sprintf("%.4f", $altd/$depth);                   #now accept everything for further filtration!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ########################### TLOD ###########################
                  my $tlod = snvCalling->calTumorLOD($phred, $localEr, 0.001, $somatic{$coor}{$djindex}{$name}, $altd, $depth);
                  ############################################################
                  $somatic{$coor}{$djindex}{$name} .= '|'.$endratio.'|'.$cmean.','.$cmedian.'|'.$strandRatio.','.$strandRatioRef.','.$stranFisherP.'|'.$badqual.'|'.$tlod.'|'.$ToxoG.'|'.$indmean.','.$indmedian.'|'.$vrlen;
                  #$cmean = 0;   #reset for artifact like stuff
                  #$cmedian = 0; #reset
                }
              }
            } else {   #vard !=0 but it is other variant bases
              $somatic{$coor}{$djindex}{$name} = 0;
              $cmean = 0;   #reset for artifact like stuff
              $cmedian = 0; #reset
            }

            ####### here for errorEst #############################
            if ($task =~ /errorEst/) {
              if ($altd > 0) {
                my $hsnpendratio = sprintf("%.4f", $vends/$vard);
                my $hsnpmaf = sprintf("%.3f", $altd/$depth);
                if ($hsnpmaf > 0.2 and $hsnpmaf < 0.8 and $hsnpendratio < 0.7 and $depth > 20 and ($cmean < 2 and $cmedian < 1.5)) { #only clean snp
                  my %baseCounts;
                  $baseCounts{'A'} = $A + $An;
                  $baseCounts{'C'} = $C + $Cn;
                  $baseCounts{'G'} = $G + $Gn;
                  $baseCounts{'T'} = $T + $Tn;
                  delete($baseCounts{$ref});
                  delete($baseCounts{$alt});
                  foreach my $errorBase (keys %baseCounts){
                    if ($baseCounts{$errorBase} < 4){
                      $errorEst{$name}{'errorBase'} += $baseCounts{$errorBase};
                    }
                  } #error bases
                  $errorEst{$name}{'totalBase'} += $depth;
                }
              }
            }
            ####### here for errorEst #############################

          } else {  #coor is not found in the original table
            print STDERR "error: coor is not found in the original file, must be found!\n";
            exit 22;
          }
        } else {   #either vard or depth == 0
          $somatic{$coor}{$djindex}{$name} = 0;
        }

        $somatic{$coor}{$djindex}{$name} .= "\t$depth";
        if ($junction != 0) {    #there are some junction reads
          $somatic{$coor}{$djindex}{$name} .= ",$junction";
        }
        $somatic{$coor}{$djindex}{'info'} = join("\t", ("ref","alt"));
        $somatic{$coor}{$djindex}{'consecutive'} .= $cmean.','.$cmedian.','.$vard.';';
        ############################ above adding data ###########################################

      SEENB:

        #redefine
        $prevCoor = $coor;
      }                         #snv
    }
    close IN;
  }

  foreach my $coor (sort {$a =~ /^(\w+):(\d+)$/; my $ca = $1; my $pa = $2; $b =~ /^(\w+):(\d+)$/; my $cb = $1; my $pb = $2; $ca cmp $cb || $pa <=> $pb} keys %somatic) {
    $coor =~ /^(\w+):(\d+)$/;
    my $chrom = $1;
    my $pos = $2;

    foreach my $djindex (sort {$a <=> $b} keys %{$somatic{$coor}}) { #each identical coordinates
      my $info = $somatic{$coor}{$djindex}{'info'};
      my @consecutive = split (';', $somatic{$coor}{$djindex}{'consecutive'});
      my $n = 0;
      my $sumCmean = 0;
      my $sumCmedian = 0;

      #my $teststring = "";
      foreach my $consecutive (@consecutive) {
        next if $consecutive eq '';
        my @tmp = split (',', $consecutive);
        next if $tmp[0] == 0;   #cmean
        next if $tmp[1] == 0;   #cmedian
        next if $tmp[2] == 0;   #vard
        $sumCmean += $tmp[2]*$tmp[0];    #vard tmp[2] used as a weight
        $sumCmedian += $tmp[2]*$tmp[1];  #vard tmp[2] used as a weight
        $n += $tmp[2];                   #total vdepth
        #$teststring .= "$tmp[2]"."\*"."$tmp[0]"." +";
      }

      if ($n > 0) {         #if you have cmean and cmedian information
        $sumCmean = sprintf("%.2f", ($sumCmean/$n));
        $sumCmedian = sprintf("%.2f", ($sumCmedian/$n));
      }

      my @information = split(",", $OR{$coor}{$djindex});
      my $id = $information[0];
      if ($type =~ 'snv') {
        $info = $information[1]."\t".$information[2];
      }
      next if ($id eq '');

      print "$chrom\t$pos\t$id\t$info";
      #foreach my $name (sort {$a =~ /($prefixReg)(\d+)?([A-Za-z0-9\-\_]+)?/; my $pa = $1; my $ia = $2; my $ias = $3; $b =~ /($prefixReg)(\d+)?([A-Za-z0-9\-\_]+)?/; my $pb = $1; my $ib = $2; my $ibs = $3; $pa cmp $pb or $ia <=> $ib or $ias cmp $ibs} keys %samples) {
      foreach my $name (sort {$a cmp $b} keys %samples) {        #strict name rule, avoid problems when names overlap
        if ($somatic{$coor}{$djindex}{$name} ne '') {
          print "\t$somatic{$coor}{$djindex}{$name}";
        } elsif ($blood eq 'yes') {
          print "\t0\t0";
        }
      }
      print "\t$sumCmean\t$sumCmedian";
      print "\n";
    } #each identical coordinates
  } #each coor
} #each chr


if ($task =~ /errorEst/) {
  open OUT, ">errorEst.tsv";
  foreach my $sample (keys %errorEst) {
    my $allbaseCount = $errorEst{$sample}{'totalBase'};
    my $errbaseCount = $errorEst{$sample}{'errorBase'};
    my $errorRate = sprintf("%.5f", $errbaseCount/$allbaseCount);
    print OUT "$sample\t$allbaseCount\t$errbaseCount\t$errorRate\n";
  }
  close OUT;
}

exit 0;



sub getchrpos {

  my $dbfile = shift;

  my $chr_old = "UNDEF";

  my %chr_start;

  my $jumper = 0;

  open DBFILE, "$dbfile" or die "The db file read error!";

  while ( <DBFILE> ) {

    if ($_ =~ /^[\#\@]/){
        $jumper = tell DBFILE;
        next;
    }
    chomp;

    my @cols = split /\t/;
    my $chr  = $cols[0];

    #$chr = 'chr'.$chr unless ($chr =~ /^chr/);
    #$chr = 'chrM' if ($chr eq 'chrMT');

    if ($chr ne $chr_old) {
        $chr_start{$chr} = $jumper;
    }

    $chr_old = $chr;
    $jumper = tell DBFILE;
  }

  close DBFILE;
  print STDERR "$dbfile chr_start loaded\n";

  return \%chr_start;
}
