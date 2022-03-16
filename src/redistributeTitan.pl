use strict;
use File::Glob ':glob';
use Data::Dumper;
use File::Basename;

my $dir = shift;
my $root = shift;

my @files = bsd_glob("$dir/*_titan");

foreach my $titan (@files){
  my  $titanbase = basename($titan);
  my $sn = "SRP";
  if ($titanbase =~ /^([A-Za-z0-9\-\_]+)\_titan$/){
    $sn = $1;
    my $cmd = "cp $titan $root/$sn/04_SNV/";
    print STDERR "$cmd\n";
    system($cmd);
  } else {
    die "strange file name: $titan\n";
  }
}

exit 0;
