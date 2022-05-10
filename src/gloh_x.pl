use strict;

my $somaticInfo = shift;
my $length = 156040895;

open IN, "$somaticInfo";
printf("%s\n", join("\t", "ID","chrom","loc-start","loc.end","num.mark","seg.mean"));

my %samples;
while ( <IN> ) {
  chomp;
  next if /^#/;
  my ($sample, $normal) = split /\t/;
  if (! exists($samples{$sample})) {
    printf("%s\n", join("\t", $sample, "X", 1, $length, 10000, 1));
    $samples{$sample} = '';
  }
  if (! exists($samples{$normal})) {
    printf("%s\n", join("\t", $normal, "X", 1, $length, 10000, 1));
    $samples{$normal} = '';
  }
}
close IN;
