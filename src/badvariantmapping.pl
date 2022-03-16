use strict;

open IN, shift;
my %tmp;
while ( <IN> ) {
  chomp;
  $_ =~ s/^(\w+)\:(\d+\t)/\1\t\2/;
  my ($chr, $pos, $chrom, $mstart, $mend, $strand, $edit, $better, $mismatches) = split(/\t/, $_);
  $chr = 'chr'.$chr;
  $chr =~ s/chrMT/chrM/;
  my $out = 0;
  if ($better > 1){
    $out = 1;
  }
  elsif ($chrom ne $chr){
    $out = 1;
  }
  elsif ( abs($pos - $mstart) > 70 or abs($pos - $mstart) < 30 ) {
    $out = 1;
  }
  if ($out == 1) {
    $tmp{$chr}{$pos} = $_;
  }
}
close IN;

foreach my $chr (sort keys %tmp){
  foreach my $pos (sort {$a <=> $b} keys %{$tmp{$chr}}){
      print "$tmp{$chr}{$pos}\n";
  }
}
