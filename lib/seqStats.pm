package seqStats;

use strict;

#
# sequencing Stats
#

sub mappingStats {

  my ($class, $bamStatsBin, $BAM, $readlen, $mappingStatsOut) = @_;

  #my $cmd = "$bamStatsBin --mapping $BAM --readlength $readlen --maxIntron 23000 --type multiMis >$mappingStatsOut";
  my $cmd = "$bamStatsBin stats $BAM >$mappingStatsOut";

  return $cmd;

}

sub xenoStats {

  my ($class, $bamStatsBin, $BAM, $outBam, $readlen, $xenoStatsOut) = @_;

  my $cmd = "$bamStatsBin --mapping $BAM --writer $outBam --readlength $readlen --maxIntron 23000 --type xeno >$xenoStatsOut";

  return $cmd;

}


sub grepStarts {

  my  ($class, $grepStartsBin, $targetRegion, $BAM, $bedCover, $chrInBam) = @_;

  my $cmd = "$grepStartsBin --region $targetRegion --mapping $BAM >$bedCover";
  if ($chrInBam ne 'SRP') {
    $cmd = "$grepStartsBin --region $targetRegion --mapping $BAM --chr $chrInBam >$bedCover";
  }

  return $cmd;

}


sub getLorenz {

  my  ($class, $lorenzCurveBin, $bedCover, $lorenzCover, $scaleFactor) = @_;

  my $cmd = "perl $lorenzCurveBin $bedCover $scaleFactor >$lorenzCover";

  return $cmd;

}

sub bed2wig {

  my ($class, $bed2wigBin, $bedCount, $wigOut) = @_;

  my $cmd = "perl $bed2wigBin $bedCount >$wigOut";

  return $cmd;

}


sub insertSize {

  my ($class, $samtoolsBin, $BAM, $outINS, $maxLine) = @_;

  my $cmd = "$samtoolsBin view -f 0x2 -F 0x400 $BAM | cut -f 9 | awk \'\$1\>0 \&\& \$1\<1000\'";

  if ($maxLine == 0){
    $cmd .= " >$outINS";
  } else {
    $cmd .= " | head -$maxLine >$outINS";
  }
  return $cmd;

}


sub plotInsertSize {

  my ($class, $RscriptBin, $insertSizeRbin, $path, $sampleName, $insGZ, $outPDF) = @_;

  my $cmd = "$RscriptBin $insertSizeRbin $path $sampleName $insGZ $outPDF";

  return $cmd;

}


1;


