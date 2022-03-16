use strict;


my $bedcov = shift;
my $scaleFactor = shift;
if ($scaleFactor eq '') {
  $scaleFactor = 1;
}


my $totalb = 0;
my $totalr = 0;


my %lorenz;
open IN, "$bedcov";
while ( <IN> ) {
  chomp;
  my ($chr, $start, $end, $depth, $read) = split /\t/;
  $depth = round($depth*$scaleFactor);
  $read = round($read*$scaleFactor);
  $lorenz{$depth}{'cov'} += 1;
  $lorenz{$depth}{'read'} += $read;
  $totalr += $read;
  $totalb += 1;
}
close IN;
print STDERR "total bases (w) is $totalb\n";
print STDERR "total reads is $totalr\n";


foreach my $depth (sort {$a <=> $b} keys %lorenz){
  my $depc = $lorenz{$depth}{'cov'};
  my $depr = $lorenz{$depth}{'read'};
  my $cumc = 0;
  my $cumr = 0;
  my $cumc2 = 0;
  my $cumr2 = 0;
  foreach my $depth2 (sort {$a <=> $b} keys %lorenz) {
    if ($depth2 > $depth) {
       $cumc += $lorenz{$depth2}{'cov'};
       $cumr += $lorenz{$depth2}{'read'};
    } elsif ($depth2 < $depth and $depth2 > 0) {
       $cumc2 += $lorenz{$depth2}{'cov'};
       $cumr2 += $lorenz{$depth2}{'read'};
    }
  }
  $cumc = sprintf("%.5f", $cumc/$totalb);
  $cumr = sprintf("%.5f", $cumr/$totalr);
  $cumc2 = sprintf("%.5f", $cumc2/$totalb);
  $cumr2 = sprintf("%.5f", $cumr2/$totalr);
  print "$depth\t$depc\t$depr\t$cumc\t$cumr\t$cumc2\t$cumr2\n";
}


sub round {
    my $number = shift;
    my $tmp = int($number);
    if ($number >= ($tmp+0.5)){
      $tmp++;
    }
    return $tmp;
}
