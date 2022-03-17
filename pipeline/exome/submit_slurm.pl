use strict;
use File::Glob ':glob';
use File::Basename;
use Data::Dumper;
use Getopt::Long;


my $somaticInfo;
my $workhard;
my $readpool;
my $config;
my $root;
my $samples;
my $prembamDir;
my $xloh;


GetOptions (
           "root|t=s"            => \$root,             #analysis dir
           "workhard|w=s"        => \$workhard,         #slurm script
           "readpool|r=s"        => \$readpool,         #where fastq is
           "somaticInfo|s=s"     => \$somaticInfo,      #somaticInfo table
           "config|c=s"          => \$config,
           "samples|p=s"         => \$samples,
           "prembamDir|b=s"      => \$prembamDir,       #old bam locate
           "xloh|x=s"            => \$xloh,             #xloh file
           "help|h"         => sub{
                               print "usage: $0 run seqare pipeline with recheck\n\nOptions:\n\t--root\t\tthe analysis dir\n";
                               print "\t--workhard\tthe slurm script\n";
                               print "\t--readpool\twhere fastq files locate\n";
                               print "\t--somaticInfo\tthe sample name table\n";
                               print "\t--config\tconfiguration file\n";
                               print "\t--samples\tthe samples needed to be RE-run, comma sperated\n";
                               print "\t--prembamDir\twhere the pre-mapping bams (for remapping, e.g., from hg19 to hg38) are located\n";
                               print "\t--xloh\t\txloh file\n";
                               #print "\t--help\t\tprint this help message\n";
                               print "\n";
                               exit 0;
                             },
           );




my @chrs;
for (1..22) {
  $_ = 'chr'.$_;
  push(@chrs, $_);
}
push(@chrs, 'chrX');
push(@chrs, 'chrY');


################ which sample to run ? #################
#my @rerun = qw();
my @rerun = split(',', $samples);
my %rerun;
foreach my $rerun (@rerun) {
  $rerun{$rerun} = '';
}
print STDERR "samples (specified) to run\n";
print STDERR Dumper(\%rerun);
########################################################


################ load all sample info ##################
my %all;
my %somaticInfo;
if ($somaticInfo and -s "$somaticInfo") {

  open IN, "$somaticInfo";
  while ( <IN> ) {
    chomp;
    s/[\s\n]$//;
    next if /^#/;
    my @columns = split /\t/;
    my $tumor = $columns[0];
    my $normal = $columns[1];

    $all{$tumor} = '';
    $all{$normal} = '' unless ($normal eq 'undef');
    $somaticInfo{$tumor} = $normal;
  }
  close IN;

}
print STDERR Dumper(\%all);
########################################################


foreach my $need (sort keys %all) {

   my $sample = $need;
   my $patient = '';
   if ($need =~ /^(D\d+)/){
     $patient = $1;
   }

   if ($samples) {
     next if (!exists($rerun{$need}));     #work for only specified samples
   }

   next if ( !exists($somaticInfo{$need}) and !exists($rerun{$need}) );     #skip germline sample
   my $normal = $somaticInfo{$need};

   my @fastqs;
   my $fastq1s = '';
   my $fastq2s = '';
   my %fastq1s;
   my %fastq2s;
   #revise to allow multiple fastq files
   if ($readpool) {
     @fastqs = bsd_glob("$readpool/$need*");
     foreach my $fastq (@fastqs) {
       if ($fastq =~ /_1\.f(ast)?q/){
         my $fastq1 = basename($fastq);
         $fastq1 =~ /^(.+?)_1\.f(ast)?q/;
         $fastq1s{$1} = $fastq1;
       } elsif ($fastq =~ /_2\.f(ast)?q/){
         my $fastq2 = basename($fastq);
         $fastq2 =~ /^(.+?)_2\.f(ast)?q/;
         $fastq2s{$1} = $fastq2;
       }
     }
   }
   foreach my $fqprefix (sort keys %fastq1s) {
     if ($fastq1s eq '') {
       $fastq1s = $fastq1s{$fqprefix};
       $fastq2s = $fastq2s{$fqprefix};
     } else {
       $fastq1s .= ",".$fastq1s{$fqprefix};
       $fastq2s .= ",".$fastq2s{$fqprefix};
     }
   }

   my @prembams;
   my $prembam;
   if ($prembamDir) {
     @prembams = bsd_glob("$prembamDir/$need*.bam");
     $prembam = $prembams[0];
   }

   foreach my $chrom ( @chrs ) {
     next if ($chrom ne 'chr1');    # donot run for each chromosome
     my $cmd;
     if ($readpool and $workhard =~ /mapping/) {   #readpool
       $cmd = "sbatch $workhard $config $root $somaticInfo $readpool $sample $fastq1s $fastq2s";    #for first step mapping
     } elsif ($prembamDir and $workhard =~ /reMapping/) {
       $cmd = "sbatch $workhard $config $root $somaticInfo $sample $prembam"; #for remapping
     } elsif ($workhard =~ /muTect|samtools/) {
       $cmd = "sbatch $workhard $config $root $somaticInfo $sample $chrom";    #for muTect calling/samtools calling
     } elsif ($workhard =~ /extractFeatures|stats|CNA/) {
       $cmd = "sbatch $workhard $config $root $somaticInfo $sample";
     } elsif ($workhard =~ /mergeCalls/) {
       $cmd = "sbatch $workhard $config $root $somaticInfo"
     } elsif ($workhard =~ /varSummary/) {
       $cmd = "sbatch $workhard $config $root $somaticInfo $xloh";     #variant classification
     }
     print "$cmd\n";
     system($cmd);
   }

}

exit 0;
