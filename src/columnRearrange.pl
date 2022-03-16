use strict;
use Data::Dumper;
use Getopt::Long;

my $file;
my $order;
my $prefix;
my $sep;

GetOptions (
           "file|f=s"       => \$file,             #filename
           "order|o=s"      => \$order,            #comma seperated indexes
           "prefix|p=s"     => \$prefix,
           "sep|s=s"        => \$sep,
           "help|h"         => sub {
                               print "usage: $0 rearrange columns according to your order or for maf rearrange\n\nOptions:\n\t--file\t\tthe filename to be reordered\n";
                               print "\t--order\t\tcomma seperated indexes\n";
                               print "\t--prefix\tthe prefix of samples' names\n";
                               print "\t--sep\tthe seperator to rearrange columns\n";
                               print "\t--help\t\tprint this help message\n";
                               print "\n";
                               exit 0;
                             },
           );


my @order = split(/\,/, $order);
my @prefix = split(',', $prefix);
my $prefixReg = join('|', @prefix);
print STDERR "prefixReg is $prefixReg\n";


open IN, "$file";
while ( <IN> ) {
  chomp;
  my @cols = split /\t/;
  if ($order ne ''){ #order
    my $print;
    foreach my $index (@order) {
      $print .= $cols[$index]."\t";
    }
    $print =~ s/\t$//;
    print "$print\n";
  } else { #order not specified, for variants maf table
    if (/^#/) {
      my %samples;
      my $secondstart = 0;

      for (my $i = 0; $i <= $#cols; $i++){       # header of the realmaf file
        if ($cols[$i] eq $sep) {
          $secondstart = $i+1;
          push(@order, $i);
          last;
        } else {
          push(@order, $i);
          if ($cols[$i] =~ /^TCGA-.+?$/ or $cols[$i] =~ /^($prefixReg)([A-Za-z0-9\-\_]+)?$/) {
            $samples{$cols[$i]} = $i;                  #store column name and array idx
          }
        } #else
      } #for

      my %inserted;
      for (my $i = $secondstart; $i <= $#cols; $i++) {
        if ( exists($samples{$cols[$i]}) ) {
          my $insertpos = $samples{$cols[$i]};
          my $offset;
          #foreach my $sample ( sort {$a =~ /($prefixReg)(\d+)?([A-Za-z0-9\-\_]+)?/; my $pa = $1; my $ia = $2; my $ias = $3; $b =~ /($prefixReg)(\d+)?([A-Za-z0-9\-\_]+)?/; my $pb = $1; my $ib = $2; my $ibs = $3; $pa cmp $pb or $ia <=> $ib or $ias cmp $ibs} keys %samples ) {
          foreach my $sample ( sort {$a cmp $b} keys %samples) {      #update rank
             my $rank = $samples{$sample};
             last if ($insertpos == $rank);              #stop at the current sample that is being added
             if ( exists($inserted{$sample}) ) {
                $offset += 2;
             }
          }
          $insertpos += $offset+1;
          splice(@order, $insertpos, 0, $i, $i+1);      #insert the two columns of realmaf
          $inserted{$cols[$i]} = '';
          $cols[$i] .= 'maf';                           #add maf suffix
        } elsif ( $cols[$i] !~ /^($prefixReg)([A-Za-z0-9\-\_]+)?d$/ ) {                #it is not depth for a sample, for residual columns
          if ( $cols[$i] =~ /^($prefixReg)([A-Za-z0-9\-\_]+)?/ ) {
            $cols[$i] .= 'maf';
            push(@order, $i);
            if ($cols[$i] =~ /^TCGA/ or $cols[$i] =~ /^($prefixReg)([A-Za-z0-9\-\_]+)?/) {
              push(@order, $i+1) if $cols[$i+1] =~ /^($prefixReg)([A-Za-z0-9\-\_]+)?d$/;      #depth info pushed
            }
          } else {
            push(@order, $i);
          }
        }
      }

    } #if it is the header

    my $print;
    foreach my $index (@order) {
      $print .= $cols[$index]."\t";
    }
    $print =~ s/\t$//;
    print "$print\n";

  } #column not specified
}
close IN;
