package cnaCalling;

use strict;

#
# cna Calling
#

sub runTitan {

  my ($class, $RscriptBin, $titanRBin, $PATH, $pbin, $sampleName, $gender, $alleleCount, $tumorWig, $normalWig, $gcWig, $mapWig, $plp, $plpe, $nc, $ncm, $sym, $transtate, $tranclone, $exons) = @_;

  my $cmd = "$RscriptBin $titanRBin $PATH $pbin $sampleName $gender $alleleCount $tumorWig $normalWig $gcWig $mapWig $plp $plpe $nc $ncm $sym $transtate $tranclone $exons";

  return $cmd;

}


1;


