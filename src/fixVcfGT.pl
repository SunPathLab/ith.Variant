use strict;

#GT GQ DP NR AD

my $vcf = shift;
my $vcfout = $vcf;
$vcfout =~ s/\.vcf/\.gtfix\.vcf/;

open IN, "$vcf";
open OUT, ">$vcfout";
while ( <IN> ){
  chomp;
  if (/^#/) {
    print OUT "$_\n";
  } else {
    my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$NORMAL,$TUMOR) = split /\t/;
    $FORMAT = 'GT:'.$FORMAT;
    $NORMAL = '0:'.$NORMAL;
    $TUMOR = '0/1:'.$TUMOR;
    printf OUT ("%s\n", join("\t", $CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,$NORMAL,$TUMOR));
  }
}
close IN;
close OUT;
