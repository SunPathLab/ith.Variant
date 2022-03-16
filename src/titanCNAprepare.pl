use strict;
use Data::Dumper;
use File::Glob ':glob';
use File::Basename;

my $data = shift;
my $somaticInfo = shift;
my $pairedCall = shift;
my $homoThred = shift;
my $ndepthThred = shift;
my $lohRegion = shift;
my $split = 1;
my $noNormal = 0;

if ($homoThred eq ''){
  $homoThred = 0.85;
}
if ($ndepthThred eq ''){
  $ndepthThred = 8;
}

my %somatic;
my %germline;  #may have multiple tumors
if ($somaticInfo ne '' and -s "$somaticInfo") {

  open IN, "$somaticInfo";
  while ( <IN> ){
    chomp;
    s/[\s\n]$//;
    my @columns = split /\t/;
    my $tumor = $columns[0];
    my $normal = $columns[1];

    $somatic{$tumor} = $normal;
    push(@{$germline{$normal}}, $tumor) if $normal ne 'undef';
  }
  close IN;
  print STDERR Dumper (\%somatic);
  print STDERR Dumper (\%germline);
}


my $outdir = dirname($data)."/titan/";
print STDERR "$outdir\n";
unless (-e "$outdir") {
  system("mkdir -p $outdir");
}


my %lohr;
if ($lohRegion ne '') {
  open LR, "$lohRegion";
  while ( <LR> ) {
    chomp;
    next if /^ID/;
    my ($sample, $chr, $start, $end, $nm, $segmean) = split /\t/;
    push(@{$lohr{$chr}}, $start.','.$end.','.$sample);
  }
  close LR;
}
print STDERR "predetermined germline LOH region:\n";
print STDERR Dumper(\%lohr);


if ($split == 1) {

  my %colnames;
  my %colindex;
  my %fhs;
  open IN, "$data";
  while (<IN>) {
    chomp;
    my @cols = split /\t/;
    if ($_ =~ /^[\#]?chr\t/) {
      $_ =~ s/^\#//;
      for (my $i = 0; $i <= $#cols; $i++) {
        $colnames{$i} = $cols[$i];
        $colindex{$cols[$i]} = $i;
      }
      next;
    } else {
      my $chr = $cols[$colindex{'chr'}];
      my $pos = $cols[$colindex{'pos'}];
      my $ref = $cols[$colindex{'ref'}];
      my $alt = $cols[$colindex{'alt'}];

      next if ($chr =~ /M(T)?$/);    #skip mitochon

      #  start loh Region sample list #
      my %lohSample;    #contain sample names of lohr for this coordinate
      if ( exists( $lohr{$chr} ) ) {
        foreach my $lregion ( @{$lohr{$chr}} ) {
          my ($lstart, $lend, $lsample) = split(',', $lregion);
          if ($pos >= $lstart and $pos <= $lend) {   #overlaps
            $lohSample{$lsample} = 1;
          }
        } #each lregion
      } #if chr is in the loh region list
      #  end loh Region sample list   #

      for (my $i = 0; $i <= $#cols; $i++) {
        if ($colnames{$i} =~ /^(.+?)maf$/) {  #now it is sample maf

          my $sample = $1;
          if ($noNormal == 1) {  #if no normal, for testing purpose only
            if ( $cols[$i] =~ /\|/ ) { #split the var surrounding information
                      my @tsinfo = split(/\|/, $cols[$i]);
                      my $tsmaf = $tsinfo[0];
                      my $tsendsratio = $tsinfo[1];
                      my ($tscmean, $tscmedian) = split(',', $tsinfo[2]);
                      my $tsd = $cols[$i+1];
                      if (($cols[$i] =~ /\|/ and $tsendsratio <= 0.9 and (($tscmean+$tscmedian) < 5.5 or $tscmedian <= 2)) or ($cols[$i] == 0 and $homoThred >= 0.85)) {  #print
                        my $fh = $sample;
                        unless (-e "$outdir/$sample\_titan") {
                          open ( my $fh, ">>", "$outdir/$sample\_titan" )  || die $!;
                          $fhs{$sample} = $fh;
                          print {$fhs{$sample}} "chr\tpos\tref\trefCount\talt\taltCount\n";
                        }
                        my $NrefCount = 0;
                        my $refCount = 0;
                        $NrefCount = round($tsmaf*$tsd);
                        $refCount = $tsd - $NrefCount;
                        if (($refCount + $NrefCount) >= 5) {
                          print {$fhs{$sample}} "$chr\t$pos\t$ref\t$refCount\t$alt\t$NrefCount\n";
                        }
                      } #true event print
print STDERR "$chr\t$pos\t$ref\t$alt\t$sample\n";
            } #split tumor info
next;
          }



          if (exists($germline{$sample})) {                                             #it is a blood/normal control, then process each tumor sample

            #  start determine if is loh Region #
            my $lohSamplePos = 'no';
            if ($lohRegion ne '') {
              if (exists($lohSample{$sample})) {   #overlaps
                $lohSamplePos = 'yes';
              }
            } #check germline loh
            if ($lohSamplePos eq 'yes') {
              print STDERR "GLOH: $sample\t$chr\t$pos\n";
            }
            #  end determine if is loh Region #

            my $calledBlood = $cols[$i-1];
            if ( $pairedCall == 1 ) {
              $calledBlood = $cols[$colindex{${$germline{$sample}}[0]}];                #paired-T-original-column
            }
            if ($calledBlood =~ /\|/) {                                                 #originally called
              my @calledBloodInfo = split(/\|/, $calledBlood);

              next if ($calledBloodInfo[2] ne '0/1' and $lohSamplePos eq 'no');         #only focus on originally hetero ones unless germline loh

              my @calledBloodRecheck = split(/\|/, $cols[$i]);                          #it is the N column rechecked
              my $calledBloodDepth = $cols[$i+1];                                       #it is the N depth column rechecked
              unless ($lohRegion ne '' or exists($somatic{$sample})) {                  #either loh region or it is both a normal and tumor (so ignore the subsequent filter)
                next if ($calledBloodRecheck[0] > $homoThred);                          #if blood has greater than 0.85 VAF, indicating wrong genotyping
              }
              next if $calledBloodDepth < $ndepthThred;                                 #if blood has too low dept


              if ($cols[$i] =~ /\|/) { #split the var surrounding information
                my @infos = split(/\|/, $cols[$i]);
                my $bmaf = $infos[0];
                my $bendsratio = $infos[1];
                my ($bcmean, $bcmedian) = split(',', $infos[2]);
                my ($strandRatio, $strandRatioRef, $strandFisherP) = split(',', $infos[3]);
                my $badQualFrac = $infos[4];
                if ($bendsratio <= 0.9 and ($strandRatio != 0 and $strandRatio != 1) and $badQualFrac < 0.6 and (($bcmean+$bcmedian) < 5.5 or $bcmedian <= 2)) { #make sure it looks real in normal

                  foreach my $tumorSamp (@{$germline{$sample}}) {   ##now should start checking for each tumor samples

                    my $indexts = $colindex{$tumorSamp.'maf'};
                    if ($cols[$indexts] =~ /\|/ or ($cols[$indexts] == 0 and $homoThred >= 0.85)) { #split the var surrounding information
                      my @tsinfo = split(/\|/, $cols[$indexts]);
                      my $tsmaf = $tsinfo[0];
                      my $tsendsratio = $tsinfo[1];
                      my ($tscmean, $tscmedian) = split(',', $tsinfo[2]);
                      my $tsd = $cols[$indexts+1];
                      if (($cols[$indexts] =~ /\|/ and $tsendsratio <= 0.9 and (($tscmean+$tscmedian) < 5.5 or $tscmedian <= 2)) or ($cols[$indexts] == 0 and $homoThred >= 0.85)) {  #print

                        my $fh = $tumorSamp;
                        unless (-e "$outdir/$tumorSamp\_titan") {
                          open ( my $fh, ">>", "$outdir/$tumorSamp\_titan" )  || die $!;
                          $fhs{$tumorSamp} = $fh;
                          print {$fhs{$tumorSamp}} "chr\tpos\tref\trefCount\talt\taltCount\n";
                        }
                        my $NrefCount = 0;
                        my $refCount = 0;
                        $NrefCount = round($tsmaf*$tsd);
                        $refCount = $tsd - $NrefCount;
                        if (($refCount + $NrefCount) >= 5) {
                          print {$fhs{$tumorSamp}} "$chr\t$pos\t$ref\t$refCount\t$alt\t$NrefCount\n";
                        }
                      } #true event print
                    } #split tumor info
                  }  ##now should start checking for each tumor samples

                } #true blood event
              } #split blood recheck info
            } #originally called
          } #blood

        } #maf
      } #each col
    } #each non header
  } #each line
  close IN;

} #split samples



sub round {
  my $number = shift;
  my $tmp = int($number);
  if ($number >= ($tmp+0.5)){
    $tmp++;
  }
  return $tmp;
}
