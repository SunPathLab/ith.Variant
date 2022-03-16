#!/usr/bin/perl
#this is a script for finding the common variations among multiple samples using a sliding window method

use strict;
use Getopt::Long;
use File::Glob ':glob';
use File::Basename;
use Data::Dumper;

my $cmdline=join(" ",@ARGV);
#print "@ $0 $cmdline\n";


my %opt = (
	   'nonrepeat'  => undef,
           'nonselfchain' => undef,
           'file' => undef,
           'bpfile' => undef,
          );

GetOptions (
	   "nonrepeat=s"    => \$opt{nonrepeat},
           "nonselfchain=s" => \$opt{nonselfchain},
           "file=s"         => \$opt{file},
           "bpfile=s"       => \$opt{bpfile},
	   "help"           => sub{
                                 print "\t--nonrepeat\toptionally choose whether allow variations in repetitive regions to be searched, repeatmask file\n";
                                 print "\t--nonselfchain\toptionally choose whether allow variations in selfchain regions to be searched, selfchain file\n";
			         print "\t--file\twhat is the file to be checked?\n";
                                 print "\t--bpfile\tbreakpointer file\n";
                                 print "\t--help\t\tprint this help message\n";
                                 print "\n";
       	                       },
           );



#############################repeat region loaded in 2D hash %repeatmask############################################
my %repeatmask;
my @rs_rm;                      #start array for each chr
my $old_chr_rm;                 #checking the chr
my $ptr_rm;                     #pointer for repeatmask sub
my $old_run_rm;

if ($opt{nonrepeat}) {
  open REP, "$opt{nonrepeat}";
  while (<REP>) {
    next if /^#/;
    chomp;
    my @tmp = split /\t/;
    (my $chr = $tmp[0]) =~ s/^chr(\w+)$/$1/;
    my $repeat_s = $tmp[3];
    my $repeat_e = $tmp[4];
    $repeatmask{$chr}{$repeat_s} = $repeat_e;
  }
  close REP;
  print STDERR "#repeatmasker file loaded\n";
}


#############################Ucsc Self Chain %selfChain#############################################
my %selfChain;
my @rs_selfChain;
my $old_chr_selfChain;
my $ptr_selfChain;
my $old_run_selfChain;

if ($opt{nonselfchain}) {
  open SELFCHAIN, "$opt{nonselfchain}";
  while ( <SELFCHAIN> ) {
    next if /^#/;
    chomp;
    my ($bin, $score, $tName, $tSize, $tStart, $tEnd, $qName, $qSize, $qStrand, $qStart, $qEnd, $id, $normScore) = split /\t/;

    $tName =~ s/^chr(\w+)$/$1/;
    $qName =~ s/^chr(\w+)$/$1/;

    next if (($normScore > 0 and $normScore < 20) or $score < 2000);

    next if ($selfChain{$tName}{$tStart} ne '' and $selfChain{$tName}{$tStart} >= $tEnd);
    $selfChain{$tName}{$tStart} = $tEnd;

    next if ($selfChain{$qName}{$qStart} ne '' and $selfChain{$qName}{$qStart} >= $qEnd);
    $selfChain{$qName}{$qStart} = $qEnd;
  }
  close SELFCHAIN;
  print STDERR "#selfchain annotation loaded\n";
}


#############################repeat region loaded in 2D hash %breakpointer############################################
my %breakpointer;
my @rs_bp;                      #start array for each chr
my $old_chr_bp;                 #checking the chr
my $ptr_bp;                     #pointer for breakpointer sub
my $old_run_bp;

if ($opt{bpfile}) {
  open BP, "gzip -dc $opt{bpfile} |";
  while (<BP>) {
    next if /^#/;
    chomp;
    my @tmp = split /\t/;
    (my $chr = $tmp[0]) =~ s/^chr(\w+)$/$1/;
    my $bp_s = $tmp[3];
    my $bp_e = $tmp[4];
    $breakpointer{$chr}{$bp_s} = $bp_e;
  }
  close BP;
  print STDERR "#breakpointer file loaded\n";
}
#################################################################################################################


open IN, "$opt{'file'}";
while ( <IN> ){
   chomp;
   if ($_ =~ /^@/){
      print "$_\n";
      next;
   }
   if ($_ =~ /^#/){
      print "$_\trep\tsc\tbp\n";
      next;
   }
   my @cols = split /\t/;
   my $chr = $cols[0];
   $chr =~ s/^chr(\w+)$/$1/;
   my $coor = $cols[1];
   my $repeatflag = $opt{'nonrepeat'}? repeatmask("run1", $chr, $coor, $coor):0;
   my $selfchainflag = $opt{'nonselfchain'}? selfChainMask("run1", $chr, $coor):0;
   my $bpflag = $opt{'bpfile'}? bpmask("run1", $chr, $coor, $coor): 0;
   print "$_\t$repeatflag\t$selfchainflag\t$bpflag\n";
}
close IN;


exit 0;

##########################################################################################################
#                                        Subroutine Region                                               #
##########################################################################################################
sub repeatmask {
    my ($run, $chr, $start, $end) = @_;
    my $flag = 0;
    if (($chr ne $old_chr_rm) or ($run ne $old_run_rm)){
       @rs_rm = sort {$a <=> $b} keys %{$repeatmask{$chr}};
       $ptr_rm = 0;
    }
    while (($ptr_rm<=$#rs_rm) and ($repeatmask{$chr}{$rs_rm[$ptr_rm]} < $start)){
      $ptr_rm++;
    }
    if ($rs_rm[$ptr_rm] <= $end){
      $flag = 1;
    }
    $old_chr_rm = $chr;
    $old_run_rm = $run;
    return $flag;
}

sub bpmask {
    my ($run, $chr, $start, $end) = @_;
    my $flag = 0;
    if (($chr ne $old_chr_bp) or ($run ne $old_run_bp)) {
       @rs_bp = sort {$a <=> $b} keys %{$breakpointer{$chr}};
       $ptr_bp = 0;
    }
    while (($ptr_bp<=$#rs_bp) and ($breakpointer{$chr}{$rs_bp[$ptr_bp]} < $start)){
      $ptr_bp++;
    }
    if ($rs_bp[$ptr_bp] <= $end){
      $flag = 1;
    }
    $old_chr_bp = $chr;
    $old_run_bp = $run;
    return $flag;
}


sub selfChainMask {
    my ($run, $chr, $coor) = @_;
    my $flag = 0;
    if (($chr ne $old_chr_selfChain) or ($run ne $old_run_selfChain)){
        @rs_selfChain = sort {$a <=> $b} keys %{$selfChain{$chr}};
        $ptr_selfChain = 0;
    }
    while (($ptr_selfChain <= $#rs_selfChain) and ($selfChain{$chr}{$rs_selfChain[$ptr_selfChain]} < $coor)){
        $ptr_selfChain++;
    }
    if ($rs_selfChain[$ptr_selfChain] <= $coor){
        $flag = 1;
    }
    $old_chr_selfChain = $chr;
    $old_run_selfChain = $run;
    return $flag;
}
