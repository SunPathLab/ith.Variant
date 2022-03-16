use strict;

my $genome = shift;
my $mut = shift;
my $type = shift;

open HS, "$genome";
my %genome = ();
my $chr = undef;
while ( <HS> ) {
  if (/^>(\w+).*?\n$/) {$chr = $1; $chr =~ s/^chr//;}
  else { s/\n//g; s/\s//g; $genome{$chr}.=$_;}
}
close HS;
print STDERR "genome loaded\n";

open IN, "$mut";

while ( <IN> ) {

  chomp;
  next if ($_ =~ /^[#@]/ or $_ =~ /^chr\t/);
  my @cols = split /\t/;
  my $chr = $cols[0];
  $chr =~ s/^chr//;
  my $pos = $cols[1];
  my $ref = $cols[3];
  next if ($ref !~ /[ACGT\-]/);
  my $alt = $cols[4];
  $ref =~ s/\-//;
  $alt =~ s/\-//;

  my $lengthRef = length($ref);
  my $add = ($lengthRef == 0)? 0:1;
  #print STDERR "$chr\t$pos\t$ref\t$alt\tadd:$add\tlref:$lengthRef\n";

  my $seq = substr($genome{$chr}, ($pos-(50+$add)), (100+$lengthRef));
  $seq =~ tr/a-z/A-Z/;
  my $realref = substr($seq, 50, 1);
  if ($realref ne $ref and $type eq 'snv') {
     print STDERR "shit_not_same_ref\t$chr\t$pos\t$ref\t$realref\n";
     #exit 22;
  }

  my $oldseq = $seq;
  my $newseq;

  $newseq = substr($seq, 50, $lengthRef, $alt);
  $newseq = $seq;

  my $identifier = $chr.":".$pos;
  print ">$identifier\n";
  print "$newseq\n";

}
close IN;

exit;
