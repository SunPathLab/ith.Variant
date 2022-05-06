#!/usr/bin/perl -w

use Getopt::Long;
use Data::Dumper;
use strict;
use File::Glob ':glob';
use File::Basename;
use FindBin qw($RealBin);
use Parallel::ForkManager;
use Statistics::Basic qw(:all);
use lib "$RealBin/../lib";
use bwaMapping;
use snvCalling;
use cnaCalling;
use seqStats;

my %options;
my %runlevel;
my %runTask;

$options{'noexecute'}   = 0;
$options{'quiet'}       = 0;
$options{'runlevels'}   = undef;
$options{'runTask'}     = undef;
$options{'readlen'}     = 0;
$options{'mapper'}      = "bwa";
$options{'sampleName'}  = 'SRP';
$options{'samplePairNames'} = "SRP";
$options{'gender'} = "female";
$options{'FASTQ1'}      = 'SRP';
$options{'FASTQ2'}      = 'SRP';
$options{'fastqFiles1'} = 'SRP';
$options{'fastqFiles2'} = 'SRP';
$options{'platform'}    = "ILLUMINA";
$options{'skipTask'}    = 'SRP';
$options{'bams'}        = 'SRP';
$options{'bamID'}       = 1;
$options{'mutectCall'}  = 'SRP';
$options{'lanepath'}    = 'SRP';
$options{'threads'}     = 1;
$options{'splitChr'}    = undef;
$options{'help'}        = undef;
$options{'qcOFF'}       = undef;
$options{'keepBedCount'}= undef;
$options{'root'}        = "$RealBin/../PIPELINE";
$options{'readpool'}    = 'SRP';
$options{'gf'}          = "png";                 #the format used in html report
$options{'bzip'}        = undef;                 #to allow bzip copressed fastq files
$options{'Rbinary'}     = 'R';
$options{'seqType'}     = 'WXS,paired-end';      #experimental types
$options{'tmpDir'}      = '';
$options{'bin'}         = "$RealBin/";
$options{'configure'}   = "SRP";
$options{'skipPileup'}  = "yes";
$options{'samSens'} = 0.005;
$options{'snvPattern'} = "SRP";                  #for grepping samtools vcf files (and other)
$options{'indelPattern'} = "SRP";                #for grepping samtools vcf files (and other)
$options{'lorenzScaleFactor'} = 1.0;

$options{'chrPrefInBam'} = "SRP";
$options{'somaticInfo'} = "SRP";
$options{'germline'}    = "SRP";
$options{'samCallmaxDepth'} = 400;
$options{'indel'}       = "SRP";
$options{'recheck'}     = "SRP";
$options{'recheckBams'} = 'SRP';
$options{'plpTitan'}    = 2.0;
$options{'plpeTitan'}   = "TRUE";
$options{'ncTitan'}     = 0.5;
$options{'ncmTitan'}    = "map";
$options{'symmetric'}   = "TRUE";
$options{'transtate'}   = 1e12;
$options{'tranclone'}   = 1e9;
$options{'maxMem'} = '4g';
$options{'javaTmp'} = '/local/scratch';
$options{'nwigString'}  = 'SRP';
$options{'homoThred'}   = 0.85;
$options{'ndepthThred'} = 8;

$options{'mergeNonsegdup'} = 1;
$options{'mergeRare'}      = 1;
$options{'qualTitan'}   = 50;
$options{'vafTitan'}  = 0.15;
$options{'rareVariants'} = undef;
$options{'germlineLOH'} = '';                  #this file contains ID for both tumors and normals
$options{'maxInsLine'} = 0;
$options{'ignoreRG'} = 0;
$options{'CfiftyOff'} = 0;

$options{'chrProcess'} = 'SRP';
$options{'chrProcessRegion'} = 'SRP';

$options{'noStrandBias'} = 'no';
$options{'customCalling'} = 'SRP';

$options{'cmeme'} = 5.5;
$options{'cmedian'} = 2;


if (@ARGV == 0) {
  helpm();
} else {
  printf STDERR "\n# $0 %s\n",join(" ",@ARGV);
}


GetOptions(
           "sampleName=s" => \$options{'sampleName'},
           "samplePairNames=s" => \$options{'samplePairNames'},
           "gender=s" => \$options{'gender'},
           "FASTQ1=s"     => \$options{'FASTQ1'},
           "FASTQ2=s"     => \$options{'FASTQ2'},
           "fastqFiles1=s"=> \$options{'fastqFiles1'},
           "fastqFiles2=s"=> \$options{'fastqFiles2'},
           "platform=s"   => \$options{'platform'},
           "bams=s"       => \$options{'bams'},            #already mapped -> halfway enter the pipe
           "bamID=s"      => \$options{'bamID'},
           "chrPrefInBam=s" => \$options{'chrPrefInBam'},
           "mutectCall=s" => \$options{'mutectCall'},      #already called -> halfway enter the pipe
           "qcOFF"        => \$options{'qcOFF'},
           "keepBedCount" => \$options{'keepBedCount'},
           "runID=s"      => \$options{'runID'},
           "runlevel=s"   => \$options{'runlevels'},
           "runTask=s"    => \$options{'runTask'},
           "skipTask=s"   => \$options{'skipTask'},
           "seqType=s"    => \$options{'seqType'},
           "noexecute"    => \$options{'noexecute'},
           "quiet"        => \$options{'quiet'},
           "maxMem=s"     => \$options{'maxMem'},
           "javaTmp=s"    => \$options{'javaTmp'},
           "splitChr"     => \$options{'splitChr'},
           "readlen=i"    => \$options{'readlen'},
           "mapper=s"     => \$options{'mapper'},
           "threads=i"    => \$options{'threads'},
           "gf=s"         => \$options{'gf'},
           "root=s"       => \$options{'root'},
           "readpool=s"   => \$options{'readpool'},
           "bzip"         => \$options{'bzip'},
           "Rbinary=s"    => \$options{'Rbinary'},
           "help|h"       => \$options{'help'},
           "configure=s"  => \$options{'configure'},
           "somaticInfo=s"=> \$options{'somaticInfo'},
           "germline=s"   => \$options{'germline'},
           "samCallmaxDepth=i" => \$options{'samCallmaxDepth'},
           "indel=s"      => \$options{'indel'},
           "lorenzScaleFactor=f" => \$options{'lorenzScaleFactor'},
           "recheck=s"    => \$options{'recheck'},
           "recheckBams=s" => \$options{'recheckBams'},
           "tmpDir=s"     => \$options{'tmpDir'},
           "qualTitan=i"  => \$options{'qualTitan'},
           "vafTitan=f" => \$options{'vafTitan'},
           "rareVariants" => \$options{'rareVariants'},
           "plpTitan=f"   => \$options{'plpTitan'},
           "plpeTitan=s"  => \$options{'plpeTitan'},
           "ncTitan=f"    => \$options{'ncTitan'},
           "ncmTitan=s"    => \$options{'ncmTitan'},
           "symmetric=s"  => \$options{'symmetric'},
           "transtate=f"  => \$options{'transtate'},
           "tranclone=f"  => \$options{'tranclone'},
           "nwigString=s" => \$options{'nwigString'},
           "homoThred=f"  => \$options{'homoThred'},
           "ndepthThred=i" => \$options{'ndepthThred'},
           "mergeNonsegdup=i" => \$options{'mergeNonsegdup'},
           "mergeRare=i"  => \$options{'mergeRare'},
           "germlineLOH=s"=> \$options{'germlineLOH'},
           "maxInsLine=i" => \$options{'maxInsLine'},
           "ignoreRG=i"   => \$options{'ignoreRG'},
           "CfiftyOff=i"  => \$options{'CfiftyOff'},
           "chrProcess=s" => \$options{'chrProcess'},
           "skipPileup=s" => \$options{'skipPileup'},
           "samSens=f"    => \$options{'samSens'},
           "snvPattern=s" => \$options{'snvPattern'},
           "indelPattern=s" => \$options{'indelPattern'},
           "noStrandBias=s" => \$options{'noStrandBias'},     #treatment of strandbias
           "customCalling=s" => \$options{'customCalling'},   #customized calling file
           "cmeme=f"      => \$options{'cmeme'},              #cmeancmedian
           "cmedian=f"    => \$options{'cmedian'}             #cmedian
          );

#print help
helpm() if ($options{'help'});


### Read configuration and set all paths----------------------------------
my %confs;
open IN, "$options{'configure'}";
while ( <IN> ) {
  chomp;
  next if /^#/;
  my @cols = split /\t/;
  $confs{$cols[0]} = $cols[1];
}
close IN;

#translate environment variable
foreach my $confele (keys %confs){
  while ($confs{$confele} =~ /\$([A-Za-z0-9]+)/g) {
    my $eleName = $1;
    my $eleTranslate;
    if (exists ($confs{$eleName})) {
      $eleTranslate = $confs{$eleName};
      $confs{$confele} =~ s/\$$eleName/$eleTranslate/;
    } else {
      die("can't translate eleName: $eleName\n");
    }
  }
}
print STDERR Dumper (\%confs);
#-------------------------------------------------------------------------

### Frequently used names-------------------------------------------------
my @chrs = split(/\,/, $confs{'chrs'});
#-------------------------------------------------------------------------

#decompression option-----------------------------------------------------
$options{'decompress'} = "gzip -d -c";
$options{'compress'}   = "gzip";
$options{'zipSuffix'}  = "gz";
if ( $options{'bzip'} ) {
  $options{'decompress'} = "bzip2 -d -c";
  $options{'compress'}   = "bzip2";
  $options{'zipSuffix'}  = "bz2";
}
#-------------------------------------------------------------------------

### Already specified full path fastq files-------------------------------
if ($options{'fastqFiles1'} ne 'SRP'){
  $options{'fastqFiles1'} =~ s/\,/ /g;
}
if ($options{'fastqFiles1'} ne 'SRP'){
  $options{'fastqFiles2'} =~ s/\,/ /g;
}
#-------------------------------------------------------------------------

### Runlevel/Task check up------------------------------------------------
if ($options{'runlevels'}) { #true runlevels
  foreach my $r (split /\,/,$options{'runlevels'}) {
    my $from=1;
    my $to=20;
    if ($r=~/^(\d+)/) {
      $from=$1;
    }
    if ($r=~/\-(\d+)$/) {
      $to=$1;
    } elsif ($r!~/\-/) {
      $to=$from;
    }
    for (my $i=$from;$i<=$to;$i++) {
      $runlevel{$i}=1;
    }
  }
}
if ($options{'runTask'}) {
  foreach my $task (split(/\,/, $options{'runTask'})) {
    $runTask{$task} = '';
  }
}
if (! $options{'runlevels'} and ! $options{'runTask'}) {
  print STDERR "no runlevel or runTask has been set, exit.\n";
  helpm();
}
#-------------------------------------------------------------------------

if ($options{'root'} eq "$RealBin/../PIPELINE") {
  if (-e "$RealBin/../PIPELINE") {
    print STDERR "no root dir given, analysis will be run under $options{'root'}.\n";
  }
  else {
    print STDERR "no root dir given, $options{'root'} does not exist, please do -h or --help to check how to set root dir.\n";
    helpm();
  }
} else {
  $options{'readpool'} = $options{'root'} if $options{'readpool'} eq 'SRP';
}


#store somatic information------------------------------------------------
my %somatic;
my %germline;                   #may have multiple tumors
if (-s "$options{'somaticInfo'}") {

  open IN, "$options{'somaticInfo'}";
  while ( <IN> ) {
    chomp;
    s/[\s\n]$//;
    my @columns = split /\t/;
    my $tumor = $columns[0];
    my $normal = $columns[1];

    $somatic{$tumor} = $normal;
    push(@{$germline{$normal}}, $tumor) if $normal ne 'undef';
  }
  close IN;
  #print STDERR Dumper (\%somatic);
  #print STDERR Dumper (\%germline);
}
#-------------------------------------------------------------------------


#-----get the chr that need to be processed-------------------------------
if ($options{'chrProcess'} ne 'SRP') {
  open IN, "$confs{'chromosomeSize'}";
  while ( <IN> ) {
    chomp;
    my ($chr, $size) = split /\t/;
    if ( $chr !~ /^chr/ ) {
      $chr = 'chr'.$chr;
    }
    my $chrWithchr = $options{'chrProcess'};
    if ( $chrWithchr !~ /^chr/ ) {
      $chrWithchr = 'chr'.$chrWithchr;
    }
    if ($chr eq $chrWithchr) {
      $options{'chrProcessRegion'} = $options{'chrProcess'}.':1-'.$size;
    }
  }
  close IN;
}
print STDERR "chr region that need to be taken care is $options{'chrProcess'}\n";
#-------------------------------------------------------------------------


###
###preparation the lane and read path enviroment
###

if ($options{'lanepath'} eq 'SRP' and $options{'sampleName'} ne 'SRP') {
  printtime();
  $options{'lanepath'} = "$options{'root'}/$options{'sampleName'}";   #define lane path
  print STDERR "####### lane name is set to $options{'sampleName'} #######\n\n";
  unless (-e "$options{'lanepath'}") {
    my $cmd = "mkdir -p $options{'lanepath'}";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
}

if ($options{'bams'} ne 'SRP') {                              #bam -> halfway enter pipe
  unless (-e "$options{'lanepath'}/02_MAPPING") {
    my $cmd = "mkdir -p $options{'lanepath'}/02_MAPPING";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  my $linkBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.br\.rmDup\.md\.bam";
  if ( exists($runTask{'remap'}) ) {                          #need to do remapping
    $linkBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.premap\.bam";
  }
  if ( exists($runTask{'indelRealignment'}) ) {
    $linkBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.bam";
  }
  if ( exists($runTask{'BaseRecalibration'}) ) {
    $linkBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.bam";
  }
  if ( exists($runTask{'MarkDuplicates'}) ) {
    $linkBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.br\.bam";
  }
  if ( exists($runTask{'recalMD'}) ) {
    $linkBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.br\.rmDup\.bam";
  }

  unless (-s "$linkBam") {
    my $cmd = "ln -s $options{'bams'} $linkBam";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  unless (-s "$linkBam\.bai" or !-s "$options{'bams'}\.bai") {  #if bam bai available and not linked
    my $cmd = "ln -s $options{'bams'}\.bai $linkBam\.bai";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  if (-s "$linkBam" and !-s "$linkBam\.bai" and ! exists($runTask{'remap'}) ) {                         # not indexed before, fix
    my $cmd = bwaMapping->bamIndex($confs{'samtoolsBin'}, $options{'threads'}, $linkBam);  # index it
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  goto REALSTEPS;
}

if ($options{'mutectCall'} ne 'SRP') {                           #mutectCall -> halfway enter pipe
  unless (-e "$options{'lanepath'}/04_SNV") {
    my $cmd = "mkdir -p $options{'lanepath'}/04_SNV";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  my $linkMutect = "$options{'lanepath'}/04_SNV/$options{'sampleName'}\.mutect";
  if ($options{'mutectCall'} =~ /\.gz$/) { #gzipped
    $linkMutect .= '.gz';
  }

  unless (-s "$linkMutect") {
    my $cmd = "ln -s $options{'mutectCall'} $linkMutect";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  goto REALSTEPS;
}

if ($options{'readpool'} ne 'SRP' and $options{'FASTQ1'} ne 'SRP' and $options{'fastqFiles1'} eq 'SRP') {

  printtime();
  print STDERR "####### preparing directories #######\n\n";

  unless (-e "$options{'lanepath'}/01_READS/") {
    my $cmd = "mkdir -p $options{'lanepath'}/01_READS/";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  my @fastqFile1 = split(/\,/, $options{'FASTQ1'});
  $options{'fastqFiles1'} =~ s/^SRP//;
  foreach my $fastqFile1 (@fastqFile1) {
    if ($fastqFile1 !~ /\.[bg]z2?$/){
      die "\[error\]: $fastqFile1 must be gzip or bzipped!\n";
    }
    $fastqFile1 = $options{'readpool'}.'/'.$fastqFile1;
    $options{'fastqFiles1'} .= $fastqFile1." ";
  }
  $options{'fastqFiles1'} =~ s/\s$//;

  if ($options{'FASTQ2'} ne 'SRP') {
    my @fastqFile2 = split(/\,/, $options{'FASTQ2'});
    $options{'fastqFiles2'} =~ s/^SRP//;
    foreach my $fastqFile2 (@fastqFile2) {
      $fastqFile2 = $options{'readpool'}.'/'.$fastqFile2;
      $options{'fastqFiles2'} .= $fastqFile2." ";
    }
    $options{'fastqFiles2'} =~ s/\s$//;
  }

  print STDERR "lanefile1:\t$options{'fastqFiles1'}\n";
  print STDERR "lanefile2:\t$options{'fastqFiles2'}\n"; #if paired end

}

if ($options{'fastqFiles1'} ne 'SRP') {
  printtime();
  print STDERR "####### preparing directories #######\n\n";

  unless (-e "$options{'lanepath'}/01_READS/") {
    my $cmd = "mkdir -p $options{'lanepath'}/01_READS/";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  foreach my $fastqFile1 (split(" ", $options{'fastqFiles1'})) {
    my $cmd = "ln -b -s $fastqFile1 $options{'lanepath'}/01_READS/";
    my $fastqFile1Basename = basename($fastqFile1);
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'}) unless (-s "$options{'lanepath'}/01_READS/$fastqFile1Basename");
  }
}

if ($options{'fastqFiles2'} ne 'SRP' and $options{'fastqFiles2'} ne 'interleaved') {
  foreach my $fastqFile2 (split(" ", $options{'fastqFiles2'})){
    my $cmd = "ln -b -s $fastqFile2 $options{'lanepath'}/01_READS/";
    my $fastqFile2Basename = basename($fastqFile2);
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'}) unless (-s "$options{'lanepath'}/01_READS/$fastqFile2Basename");
  }
}



REALSTEPS:
############################ determine read length if not provided ###################################
if ( $options{'readlen'} == 0 and $options{'sampleName'} ne 'SRP') { #read length not set
  if ($options{'fastqFiles1'} ne 'SRP') {
    my @fastqFiles1Temp = split(/\s/, $options{'fastqFiles1'});
    if ( -s "$fastqFiles1Temp[0]" ) {
      my $first_second_line = `$options{'decompress'} "$fastqFiles1Temp[0]" | head -2 | grep -v "^@"`;
      $options{'readlen'} = length($first_second_line) - 1;
    }
  }
  if ( $options{'readlen'} == 0 ) {  #still zero, check bams
    my @bamTmp = bsd_glob("$options{'lanepath'}/02_MAPPING/*.bam");
    foreach my $bamTmp (@bamTmp){
      if (-s "$bamTmp"){
        my $samSix = `samtools view $bamTmp \| awk \-F\"\t\" \'\{print \$6\}\' \| awk \'\$1 \!\~ \/\[IDNHS\\\*\]\/\' \| head \-1000 \| tr \"\\n\" \"\,\"`;
        chomp($samSix);
        my @matchLen;
        foreach my $matchLen (split(/\,/, $samSix)){
          $matchLen =~ /^(\d+)M$/;
          $matchLen = $1;
          push(@matchLen, $matchLen);
        }
        $options{'readlen'} = median(\@matchLen);   #set length to the median of the first 1000 reads with complete matching
      }
      if ($options{'readlen'} > 0){
        last;
      }
    } #loop each bams to find existing bam
  } #check bam
  print STDERR "read length is not set, will take the original read length ($options{'readlen'} bp)\n";
}

#######################################################################################################

###
###runlevel1: QC
###

my $runlevels = 1;
if (exists($runlevel{$runlevels}) or exists($runTask{'QC'})) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  my $qc_out1 = "$options{'lanepath'}/01_READS/mate1.qc";    ######
  my $qc_out2;

  if ($options{'fastqFiles2'} ne 'SRP' and $options{'fastqFiles2'} ne 'interleaved') {
    $qc_out2 = "$options{'lanepath'}/01_READS/mate2.qc";
  }

  unless ((-e "$qc_out1")) {
    my $cmd = "$options{'decompress'} $options{'fastqFiles1'} | $options{'bin'}/fastx_quality_stats -Q33 -o $qc_out1";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'}) unless ($options{'qcOFF'});
  }
  unless (($options{'fastqFiles2'} eq 'SRP' or $options{'fastqFiles2'} eq 'interleaved') or ($qc_out2 ne '' and -e "$qc_out2")) {
    my $cmd = "$options{'decompress'} $options{'fastqFiles2'} | $options{'bin'}/fastx_quality_stats -Q33 -o $qc_out2";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'}) unless ($options{'qcOFF'});
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";

}


###
###runlevel2: do the mapping and generate the statistics
###


$runlevels = 2;
if (exists($runlevel{$runlevels}) or exists($runTask{'mapping'}) or exists($runTask{'indelRealignment'}) or exists($runTask{'BaseRecalibration'}) or exists($runTask{'MarkDuplicates'}) or exists($runTask{'recalMD'})) {

  my $premapBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.premap\.bam";
  my $rawBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.bam";
  my $sortedBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.bam";
  #my $irBam = ($options{'chrProcess'} eq 'SRP')? "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.bam" : "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.$options{'chrProcess'}\.bam";
  my $irBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.bam";
  my $brBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.br\.bam";
  my $rmDupBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.br\.rmDup\.bam";
  my $finalBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.br\.rmDup\.md\.bam";

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  unless (-e "$options{'lanepath'}/02_MAPPING") {
    my $cmd = "mkdir -p $options{'lanepath'}/02_MAPPING";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  my @allBams = bsd_glob("$options{'lanepath'}/02_MAPPING/$options{'sampleName'}*\.bam");
  if ($#allBams == -1 or exists($runTask{'mapping'}) or exists($runTask{'remap'})) {
    my $ReadGroup = '@RG'."\tID:".$options{'bamID'}."\tSM\:".$options{'sampleName'}."\tPL\:".$options{'platform'};
    my $cmd;
    if (exists($runTask{'remap'}) and -e "$premapBam") {
      $cmd = bwaMapping->bwaRemapping($confs{'bwaBin'}, $confs{'samtoolsBin'}, $ReadGroup, $options{'threads'}, $confs{'BWAINDEX'}, $rawBam, $premapBam);
    } elsif ($options{'fastqFiles2'} eq 'interleaved') {  #need smart paring
      $cmd = bwaMapping->bwaSmartMapping($confs{'bwaBin'}, $confs{'samtoolsBin'}, $ReadGroup, $options{'threads'}, $confs{'BWAINDEX'}, $rawBam, $options{'fastqFiles1'});
    } elsif ($options{'fastqFiles2'} eq 'SRP') {     #single end
      $cmd = bwaMapping->bwaSingleMapping($confs{'bwaBin'}, $confs{'samtoolsBin'}, $ReadGroup, $options{'threads'}, $confs{'BWAINDEX'}, $rawBam, $options{'fastqFiles1'});
    } else {                                         #paired-end
      $cmd = bwaMapping->bwaPairMapping($confs{'bwaBin'}, $confs{'samtoolsBin'}, $ReadGroup, $options{'threads'}, $confs{'BWAINDEX'}, $rawBam, $options{'fastqFiles1'}, $options{'fastqFiles2'});
    }
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  if ( $options{'seqType'} =~ /xenograft/) {  #need to filter out mouse reads

    my $xenoStats = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.xenoStats";
    my $noMouseBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.noMouse\.bam";
    my $cmd = seqStats->xenoStats("$options{'bin'}/Rseq_bam_stats", $rawBam,  $noMouseBam, $options{'readlen'}, $xenoStats);
    unless (-e "$xenoStats") {
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    if (-s "$noMouseBam" and -s "$rawBam") {
      (my $originalBam = $rawBam) =~ s/\.bam/\.original\.bam/;
      $cmd = "mv $rawBam $originalBam -f";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      $cmd = "mv $noMouseBam $rawBam -f";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
  }

  if ($options{'seqType'} =~ /WGS/ and $options{'seqType'} !~ /ignore/) {
    print STDERR "stop for disk space checking\n";
    exit 0;
  }

  unless (-s "$finalBam" and !(exists($runTask{'indelRealignment'}) or exists($runTask{'MarkDuplicates'}) or exists($runTask{'recalMD'}) or exists($runTask{'BaseRecalibration'})) ) {    #processing bam
    if (-s "$rawBam" and !(-s "$sortedBam")) {     #must sort
      my $cmd = bwaMapping->bamSort($confs{'samtoolsBin'}, $options{'threads'}, $rawBam."\.tmp", $sortedBam, $rawBam);
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      $cmd = bwaMapping->bamIndex($confs{'samtoolsBin'}, $options{'threads'}, $sortedBam);     #index it
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    if (-s "$rawBam" and -s "$sortedBam") {
      my $cmd = "rm $rawBam -f";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    if ((-s "$sortedBam" and !(-s "$irBam")) or exists($runTask{'indelRealignment'})) { #indel realignment
      if ($options{'skipTask'} !~ /indelRealignment/) {
        my $indelTargetList = $sortedBam."\.target_intervals.list";
        my $CHR = 'ALL';
        #if ( $options{'chrProcess'} ne 'SRP' ) {
        #  $CHR = $options{'chrProcess'};
        #  $indelTargetList = $sortedBam."\.$CHR\.target_intervals.list";
        #}
        if ($options{'splitChr'}) { #if split chr, folk it up, not for hpc clusters, only for workstations, need to be merged later
          my $chrBatches = partitionArray(\@chrs, $options{'threads'});
          foreach my $chrBatch (@{$chrBatches}) {
            my $manager = Parallel::ForkManager->new($options{'threads'});
            my $processedChroms = "chromosome ";
            foreach my $chrom (@{$chrBatch}) {
              $manager->start and next;
              $CHR = $chrom;
              $indelTargetList =~ s/\.target_intervals.list/\.$CHR\.target_intervals.list/;
              $irBam =~ s/\.bam/\.$CHR\.bam/;
              my $cmd;
              unless (-s "$indelTargetList") {
                $cmd = bwaMapping->indelRealignment1($confs{'gatkBin'}, $sortedBam, $confs{'GFASTA'}, $confs{'KNOWNINDEL1'}, $confs{'KNOWNINDEL2'}, $CHR, $indelTargetList, $options{'threads'}, $options{'maxMem'}, $options{'javaTmp'}, $confs{'java'});
                RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
              }
              $cmd = bwaMapping->indelRealignment2($confs{'gatkBin'}, $sortedBam, $confs{'GFASTA'}, $indelTargetList, $confs{'KNOWNINDEL1'}, $confs{'KNOWNINDEL2'}, $CHR, $irBam, $options{'threads'}, $options{'maxMem'}, $options{'javaTmp'}, $confs{'java'});
              RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
              $cmd = bwaMapping->bamIndex($confs{'samtoolsBin'}, $options{'threads'}, $irBam); #index it
              RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
              $manager->finish;
              $processedChroms .= $chrom.',';
            }
            $manager->wait_all_children;
            print STDERR "$processedChroms have been processed!\n";
          }
        } else {
          my $cmd;
          unless (-s "$indelTargetList") {
            $cmd = bwaMapping->indelRealignment1($confs{'gatkBin'}, $sortedBam, $confs{'GFASTA'}, $confs{'KNOWNINDEL1'}, $confs{'KNOWNINDEL2'}, $CHR, $indelTargetList, $options{'threads'}, $options{'maxMem'}, $options{'javaTmp'}, $confs{'java'});
            RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
            if ($options{'skipTask'} =~ /irStep2/) {
              print STDERR "indel realign list produced. stop now\n";
              exit 0;
            }
          }
          $cmd = bwaMapping->indelRealignment2($confs{'gatkBin'}, $sortedBam, $confs{'GFASTA'}, $indelTargetList, $confs{'KNOWNINDEL1'}, $confs{'KNOWNINDEL2'}, $CHR, $irBam, $options{'threads'}, $options{'maxMem'}, $options{'javaTmp'}, $confs{'java'});
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
          $cmd = bwaMapping->bamIndex($confs{'samtoolsBin'}, $options{'threads'}, $irBam); #index it
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        }
      } else {
        my $cmd = "mv $sortedBam $irBam";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        $cmd = bwaMapping->bamIndex($confs{'samtoolsBin'}, $options{'threads'}, $irBam); #index it
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    }
    #if ($options{'chrProcess'} ne 'SRP') {  #for each chromosome
    #  print STDERR "stop for merging irBams\n";
    #  exit 0;
    #}
    if (-s "$sortedBam" and -s "$irBam") {
      my $cmd = "rm $sortedBam $sortedBam\.bai -f";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$irBam") {               #remove redundant bai
      (my $redBai = $irBam.'.bai') =~ s/\.bam\.bai/\.bai/;
      if (-s "$redBai") {
        my $cmd = "rm $redBai -f";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    }

    if ((-s "$irBam" and !(-s "$brBam")) or exists($runTask{'BaseRecalibration'})) {   #base recalibration
      my $cmd;
      if ($options{'skipTask'} !~ /BaseRecalibration/) {             #if skipped
        my $brTable = $irBam.".baseRecal.table";
        unless (-s "$brTable") {
          $cmd = bwaMapping->BaseRecalibration($confs{'gatkBin'}, $irBam, $confs{'GFASTA'}, $confs{'muTectDBSNP'}, $confs{'KNOWNINDEL1'}, $confs{'KNOWNINDEL2'}, $brTable, $options{'threads'}, $options{'maxMem'}, $options{'javaTmp'}, $confs{'java'});
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        }
        $cmd = bwaMapping->BaseRecalibrationPrint($confs{'gatkBin'}, $irBam, $confs{'GFASTA'}, $brTable, $brBam, $options{'threads'}, $options{'maxMem'}, $options{'javaTmp'}, $confs{'java'});
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      } else {
        $cmd = "mv $irBam $brBam";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
      $cmd = bwaMapping->bamIndex($confs{'samtoolsBin'}, $options{'threads'}, $brBam); #index it
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$brBam" and (-s "$irBam" or -s "$irBam\.bai")) {
      my $cmd = "rm $irBam $irBam\.bai -f";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$brBam") {          #remove redundant bai
      (my $redBai = $brBam.'.bai') =~ s/\.bam\.bai/\.bai/;
      if (-s "$redBai") {
        my $cmd = "rm $redBai -f";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    }
    if ( exists($runTask{'BaseRecalibrationOnly'}) ) {
      exit 0;
    }

    if ((-s "$brBam" and !(-s "$rmDupBam")) or exists($runTask{'MarkDuplicates'})) {  #rmDup
      my $rmDupMetric = $brBam.".rmDupMetric";
      my $cmd = bwaMapping->MarkDuplicates($confs{'MarkDuplicatesBin'}, $brBam, $rmDupBam, $rmDupMetric, $options{'maxMem'}, $options{'javaTmp'}, $confs{'java'});
      if ($options{'skipTask'} =~ /MarkDuplicates/) {
        $cmd = "mv $brBam $rmDupBam";
      }
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$brBam" and -s "$rmDupBam") {
      my $cmd = "rm $brBam $brBam\.bai -f";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    if ((-s "$rmDupBam" and !(-s "$finalBam")) or exists($runTask{'recalMD'})) {
      if ($options{'skipTask'} !~ /recalMD/) {
        my $cmd = bwaMapping->recalMD($confs{'samtoolsBin'}, $rmDupBam, $confs{'GFASTA'}, $finalBam, $options{'threads'});
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        $cmd = bwaMapping->bamIndex($confs{'samtoolsBin'}, $options{'threads'}, $finalBam); #index it
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      } else {
        my $cmd = "mv $rmDupBam $finalBam";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        $cmd = bwaMapping->bamIndex($confs{'samtoolsBin'}, $options{'threads'}, $finalBam); #index it
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    }
    if (-s "$rmDupBam" and -s "$finalBam") {
      my $cmd = "rm $rmDupBam $rmDupBam\.bai -f";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
  }

  if (-s "$finalBam" and !-s "$finalBam\.bai") {
    my $cmd = bwaMapping->bamIndex($confs{'samtoolsBin'}, $options{'threads'}, $finalBam); #index it
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";
}


###
###runlevel3: STATS
###

$runlevels = 3;
if (exists $runlevel{$runlevels}) {

  unless (-e "$options{'lanepath'}/03_STATS") {
    my $cmd = "mkdir -p $options{'lanepath'}/03_STATS";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  #basic read counting stats:
  my $mappingStats = "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.mapping.stats";
  my $finalBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.br\.rmDup\.md\.bam";
  my $statBam = ($options{'recheckBams'} eq "SRP")? $finalBam : $options{'recheckBams'};

  if ( ! -e "$finalBam" ) {
    $finalBam = $statBam;
  }

  unless (-s "$mappingStats") {
    #my $cmd = seqStats->mappingStats("$options{'bin'}/Rseq_bam_stats", $statBam, $options{'readlen'}, $mappingStats);
    my $cmd = seqStats->mappingStats($confs{'samtoolsBin'}, $statBam, $options{'readlen'}, $mappingStats);
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  #loren curve
  my $bedCover = "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.bedcoverNoDup";
  my $lorenzCover = "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.lorenzNoDup";
  unless (-s "$lorenzCover") {
    unless (-s "$bedCover") {
      my $cmd = seqStats->grepStarts("$options{'bin'}/grep_starts", $confs{'targetRegion'}, $statBam, $bedCover, $options{'chrPrefInBam'});
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    my $cmd = seqStats->getLorenz("$options{'bin'}/lorenzCurveNGS.pl", $bedCover, $lorenzCover, $options{'lorenzScaleFactor'});
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  if (-s "$bedCover" and -s "$lorenzCover") {
    my $cmd = "rm $bedCover -f";
    #RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  #for titanCNA
  my $bedCount = "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.w1k.count";
  my $wigOut = "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.wig";
  unless (-s "$wigOut") {
    unless (-s "$bedCount") {
      my $cmd = seqStats->grepStarts("$options{'bin'}/grep_starts", $confs{'w1kBed'}, $finalBam, $bedCount, $options{'chrPrefInBam'});
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    my $cmd = seqStats->bed2wig("$options{'bin'}/bed2wig.pl", $bedCount, $wigOut);
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  if (-s "$bedCount" and -s "$wigOut") {
    my $cmd = "rm $bedCount -f";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'}) unless ($options{'keepBedCount'});
  }

  #for insert size
  if ($options{'seqType'} =~ /paired/) { #do the insert size only if it is paired-end
    unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ins\.gz") {
      unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ins") {
        my $insertBam = $finalBam;
        if ( $finalBam !~ /\.bam$/ ) {  #likely a file of file names
          $insertBam = `head -1 $finalBam`;
          $insertBam =~ s/[\s\n]//;
        }
        my $cmd = seqStats->insertSize($confs{'samtoolsBin'}, $insertBam, "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ins", $options{'maxInsLine'});
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
      if (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ins") {
        my $cmd = "gzip $options{'lanepath'}/03_STATS/$options{'sampleName'}\.ins";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    }

    unless (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.insertSize\.pdf") {
      my $cmd = seqStats->plotInsertSize($confs{'RscriptBin'}, "$options{'bin'}/insertSize.R", "$options{'lanepath'}/03_STATS/", $options{'sampleName'},
                                         "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.ins\.gz", "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.insertSize\.pdf");
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
  } #insert size


  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";

}


###
###runlevel4: SNV calling
###

$runlevels = 4;
if (exists($runlevel{$runlevels}) or exists($runTask{'recheck'}) or exists($runTask{'annovar'}) or exists($runTask{'mergeCustomCallingChr'})) {

  unless (-e "$options{'lanepath'}/04_SNV") {
    my $cmd = "mkdir -p $options{'lanepath'}/04_SNV";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";


  #my $finalBam = ($options{'splitChr'})?"$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.$chrs[0]\.rmDup\.bam":"$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.rmDup\.bam";
  my $finalBam = "$options{'lanepath'}/02_MAPPING/$options{'sampleName'}\.sorted\.ir\.br\.rmDup\.md\.bam";
  my $normalBam;

  if (!exists($runlevel{$runlevels}) and exists($runTask{'recheck'})) {
    goto RECHECK;
  }

  if (!exists($runlevel{$runlevels}) and exists($runTask{'annovar'})) {
    goto ANNOVAR;
  }

  if (!exists($runlevel{$runlevels}) and exists($runTask{'mergeCustomCallingChr'})) {
    goto MERGECUSTOM;
  }

  if ($options{'somaticInfo'} eq "SRP"){
    print STDERR "ERROR: somaticInfo is not provided! Must set for somatic calling!\n";
    exit 22;
  } elsif ( !exists( $somatic{$options{'sampleName'}} ) ) {
    print STDERR "ERROR: $options{'sampleName'} is not in the somatic hash table!\n";
  } else { #get normal bam
    my $normalSampleName = $somatic{$options{'sampleName'}};
    if ($options{'samplePairNames'} eq 'redefine') {
      $options{'samplePairNames'} = $options{'sampleName'}.','.$normalSampleName;                   #redefine sample pair name, default SRP not redefine
      print STDERR "samplePairNames: $options{'samplePairNames'}\n";
    }
    $normalBam = "$options{'root'}/$normalSampleName/02_MAPPING/$normalSampleName\.sorted\.ir\.br\.rmDup\.md\.bam";
    unless (-s "$normalBam") {
      print STDERR "ERROR (warning): $normalBam is not found, please run mapping and processing for $normalSampleName!!\n";
      exit 22;
    }
  }

  if ($options{'germline'} =~ /samtoolsOnly/) {
    goto GERMLINE;
  }

  if ($options{'indel'} =~ /strelkaOnly/) {
    goto INDEL;
  }

  my $muTectOut = ($options{'chrProcess'} eq 'SRP')? "$options{'lanepath'}/04_SNV/$options{'sampleName'}\.mutect" : "$options{'lanepath'}/04_SNV/$options{'sampleName'}\.$options{'chrProcess'}\.mutect";
  my $vcfOutTmp = $muTectOut.'.vcf';
  my $vcfOut = $muTectOut.'.genome.vcf';
  my $vcfOutSorted = $muTectOut.'.genome.sorted.vcf';
  my $vcfMultiAnno = $vcfOutSorted."\.$confs{'species'}_multianno.txt";
  my $vcfMultiAnnoVCF = $vcfOutSorted."\.$confs{'species'}_multianno.vcf";
  my $vcfMultiAnnoMod = $vcfOutSorted."\.$confs{'species'}_multianno.mod.vcf";

  if (exists($runTask{'mergeMutectChr'})) {
    goto MERGE;
  }

  unless ((-s "$muTectOut" or -s "$muTectOut\.gz") or !exists( $somatic{$options{'sampleName'}} ) ) {
    my $cmd = snvCalling->muTectCalling($confs{'muTectBin'}, $finalBam, $normalBam, $confs{'GFASTA'}, $confs{'muTectCOSMIC'}, $confs{'muTectDBSNP'}, $muTectOut, $vcfOutTmp, $options{'chrProcessRegion'}, $options{'maxMem'}, $options{'javaTmp'}, $confs{'java'});
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  #annoVar annotate---------------------------------------------------------------------
  if ((-s "$muTectOut" or -s "$muTectOut\.gz") and !-s "$vcfMultiAnnoMod" and exists( $somatic{$options{'sampleName'}} ) ) {

    my $cmd = snvCalling->muTect2vcf("$options{'bin'}/mutect2vcf.pl", $muTectOut, $vcfOut, $options{'samplePairNames'});                                                #convert mutect 2 vcf
    if (-s "$muTectOut\.gz") {
      $cmd = snvCalling->muTect2vcf("$options{'bin'}/mutect2vcf.pl", $muTectOut.".gz", $vcfOut, $options{'samplePairNames'});
    }
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

    $cmd = snvCalling->vcfSort($confs{'vcfSortBin'}, $vcfOut, $vcfOutSorted);                                                                                    #sort vcf
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

    $cmd = snvCalling->runAnnovar("$confs{'ANNOVARDIR'}/table_annovar.pl", $vcfOutSorted, $confs{'ANNOVARDB'}, $confs{'species'});       #table annovar
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

    $cmd = snvCalling->convertVCFannovar("$options{'bin'}/convert_annovar_vcf.pl", $vcfMultiAnno, $vcfMultiAnnoVCF);
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

  }

  #rm temporary files
  if (-s "$vcfMultiAnnoMod" and -s "$vcfOutTmp") {
    my $cmd = "rm -rf $vcfOutTmp $vcfOutTmp\.idx";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  if (-s "$vcfMultiAnnoMod" and -s "$vcfOut") {
    my $cmd = "rm -rf $vcfOut";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  if (-s "$vcfMultiAnnoMod" and -s "$vcfOutSorted\.avinput") {
    my $cmd = "rm -rf $vcfOutSorted\.avinput";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  if (-s "$vcfMultiAnnoMod" and -s "$vcfOutSorted\.invalid_input") {
    my $cmd = "rm -rf $vcfOutSorted\.invalid_input";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  if (-s "$vcfMultiAnnoMod" and -s "$vcfOutSorted\.refGene.invalid_input") {
    my $cmd = "rm -rf $vcfOutSorted\.refGene.invalid_input";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }
  if (-s "$vcfMultiAnnoMod" and -s "$vcfMultiAnno") {
    my $cmd = "rm -rf $vcfMultiAnno $vcfMultiAnnoVCF";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  #------------------------------------------------------------------------------------

 GERMLINE:

  if ($options{'germline'} =~ /samtools/) {

    my $vcfOut = ($options{'chrProcess'} eq 'SRP')? "$options{'lanepath'}/04_SNV/$options{'sampleName'}\.samtools.genome.vcf" : "$options{'lanepath'}/04_SNV/$options{'sampleName'}\.$options{'chrProcess'}\.samtools.genome.vcf";
    (my $vcfOutSorted = $vcfOut) =~ s/\.vcf$/.sorted.vcf/;
    my $vcfMultiAnno = $vcfOutSorted."\.$confs{'species'}_multianno.txt";
    my $vcfMultiAnnoVCF = $vcfOutSorted."\.$confs{'species'}_multianno.vcf";
    my $vcfMultiAnnoMod = $vcfOutSorted."\.$confs{'species'}_multianno.mod.vcf";
    my $vcfMultiAnnoModsnv = $vcfOutSorted."\.$confs{'species'}_multianno.mod.vcf.snv";
    my $vcfMultiAnnoModindel = $vcfOutSorted."\.$confs{'species'}_multianno.mod.vcf.indel";

    if (exists($runTask{'mergeSamtoolsChr'})) {
      goto MERGESAMTOOLS;
    }

    unless (-s "$vcfOut" or -s "$vcfOutSorted" or -s "$vcfMultiAnnoMod" or -s "$vcfMultiAnnoModsnv") {
      my $cmd = snvCalling->samtoolsCalling($confs{'samtoolsBin'}, $confs{'bcftoolsBin'}, $finalBam, $normalBam, $confs{'GFASTA'}, $vcfOut, $options{'samCallmaxDepth'}, $options{'ignoreRG'}, $options{'chrProcessRegion'}, $options{'samSens'}, $options{'CfiftyOff'});
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    #annoVar annotate---------------------------------------------------------------------
    if ((-s "$vcfOut" or -s "$vcfMultiAnnoMod") and !-s "$vcfMultiAnnoModsnv") {

      unless (-s "$vcfMultiAnnoMod") {

        my $cmd = snvCalling->vcfSort($confs{'vcfSortBin'}, $vcfOut, $vcfOutSorted);      #sort vcf
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

        $cmd = snvCalling->runAnnovar("$confs{'ANNOVARDIR'}/table_annovar.pl", $vcfOutSorted, $confs{'ANNOVARDB'}, $confs{'species'}); #table annovar
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

        $cmd = snvCalling->convertVCFannovar("$options{'bin'}/convert_annovar_vcf.pl", $vcfMultiAnno, $vcfMultiAnnoVCF);
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

      }

      my $cmd = snvCalling->grepSNVvcf($vcfMultiAnnoMod, $vcfMultiAnnoModsnv, $options{'snvPattern'});
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

      $cmd = snvCalling->grepINDELvcf($vcfMultiAnnoMod, $vcfMultiAnnoModindel, $options{'indelPattern'});
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

    }

    if ( $options{'samplePairNames'} =~ /\,/ ) { #rename samples in vcf
      my @samplePariNames = split( ',', $options{'samplePairNames'} );
      my $snsreplacement = join("\t", @samplePariNames);
      foreach my $vcfsamresult (($vcfMultiAnnoModsnv,$vcfMultiAnnoModindel)) {
        open IN, "$vcfsamresult";
        open OUT, ">$vcfsamresult\.2";
        while ( <IN> ) {
          chomp;
          if (/^#CHROM\t/) {
            s/\S+\t\S+$/$snsreplacement/;
          }
          print OUT "$_\n";
        }
        close IN;
        my $cmd = "mv $vcfsamresult\.2 $vcfsamresult -f";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }  #snv and indel files
    }

    #rm temporary files
    if (-s "$vcfMultiAnnoModsnv" and -s "$vcfOut") {
      my $cmd = "rm -rf $vcfOut";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$vcfMultiAnnoModsnv" and -s "$vcfOutSorted\.avinput") {
      my $cmd = "rm -rf $vcfOutSorted\.avinput";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$vcfMultiAnnoModsnv" and -s "$vcfOutSorted\.invalid_input") {
      my $cmd = "rm -rf $vcfOutSorted\.invalid_input";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$vcfMultiAnnoModsnv" and -s "$vcfOutSorted\.refGene.invalid_input") {
      my $cmd = "rm -rf $vcfOutSorted\.refGene.invalid_input";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$vcfMultiAnnoModsnv" and -s "$vcfMultiAnno") {
      my $cmd = "rm -rf $vcfMultiAnno $vcfMultiAnnoVCF";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$vcfMultiAnnoModsnv" and -s "$vcfMultiAnnoMod") {
      my $cmd = "rm -rf $vcfMultiAnnoMod";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    #------------------------------------------------------------------------------------

  MERGESAMTOOLS:

    if (exists($runTask{'mergeSamtoolsChr'})) {
      my @samtoolsChrSnvs = bsd_glob("$options{'lanepath'}/04_SNV/*samtools*.mod.vcf.snv");
      my $samtoolsChrSnvs = join(',', @samtoolsChrSnvs);
      my @samtoolsChrIndels = bsd_glob("$options{'lanepath'}/04_SNV/*samtools*.mod.vcf.indel");
      my $samtoolsChrIndels = join(',', @samtoolsChrIndels);
      unless (-s "$vcfMultiAnnoModsnv") {
        my $cmd = "perl $options{'bin'}/mergeMutFiles.pl $samtoolsChrSnvs >$vcfMultiAnnoModsnv";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        if (-s "$vcfMultiAnnoModsnv" and -s "$samtoolsChrSnvs[0]") {
          my $cmd = "rm ".join(" ", @samtoolsChrSnvs).' -rf';
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        }
      }
      unless (-s "$vcfMultiAnnoModindel") {
        my $cmd = "perl $options{'bin'}/mergeMutFiles.pl $samtoolsChrIndels >$vcfMultiAnnoModindel";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        if (-s "$vcfMultiAnnoModindel" and -s "$samtoolsChrIndels[0]") {
          my $cmd = "rm ".join(" ", @samtoolsChrIndels).' -rf';
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        }
      }
    }

  } #germline calling with samtools


 INDEL:

  if ($options{'indel'} =~ /strelka/) {                   #do somatic small indel calling by strelka

    my $strelkaOutDir = "$options{'lanepath'}/04_SNV/strelka";
    my $vcfOutRaw = "$options{'lanepath'}/04_SNV/strelka/results/all.somatic.indels.vcf";
    my $vcfOut = "$options{'lanepath'}/04_SNV/strelka/results/all.somatic.indels.gtfix.vcf";
    (my $vcfOutSorted = $vcfOut) =~ s/\.vcf$/.sorted.vcf/;
    my $vcfMultiAnno = $vcfOutSorted."\.$confs{'species'}_multianno.txt";
    my $vcfMultiAnnoVCF = $vcfOutSorted."\.$confs{'species'}_multianno.vcf";
    my $vcfMultiAnnoMod = $vcfOutSorted."\.$confs{'species'}_multianno.mod.vcf";

    #if (exists($runTask{'mergeStrelkaChr'})) {
    #  goto MERGESTRELKA;
    #}

    unless (-s "$vcfOut" or -s "$vcfOutSorted" or -s "$vcfMultiAnnoMod") {
      unless (-s "$vcfOutRaw") {
        my $cmd = snvCalling->strelkaCalling1($confs{'strelkaBin'}, $normalBam, $finalBam, $confs{'GFASTA'}, $confs{'strelkaConfig'}, $strelkaOutDir);   #configure
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        $cmd = snvCalling->strelkaCalling2($strelkaOutDir, $options{'threads'});
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
      if (-s "$vcfOutRaw" and -e "$strelkaOutDir/chromosomes") {
        my $cmd = "rm $strelkaOutDir/chromosomes/ -rf";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
      my $cmd = snvCalling->vcfFixGT("$options{'bin'}/fixVcfGT.pl", $vcfOutRaw);
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    #annoVar annotate---------------------------------------------------------------------
    if (-s "$vcfOut" or -s "$vcfMultiAnnoMod") {

      unless (-s "$vcfMultiAnnoMod") {

        my $cmd = snvCalling->vcfSort($confs{'vcfSortBin'}, $vcfOut, $vcfOutSorted);      #sort vcf
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

        $cmd = snvCalling->runAnnovar("$confs{'ANNOVARDIR'}/table_annovar.pl", $vcfOutSorted, $confs{'ANNOVARDB'}, $confs{'species'}); #table annovar
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

        $cmd = snvCalling->convertVCFannovar("$options{'bin'}/convert_annovar_vcf.pl", $vcfMultiAnno, $vcfMultiAnnoVCF);
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

      }

    }

    #make a link
    my $cmd = "ln -s $vcfMultiAnnoMod $options{'lanepath'}/04_SNV/$options{'sampleName'}.strelka.indel.genome.sorted.vcf.$confs{'species'}_multianno.mod.vcf";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

    #rm temporary files
    if (-s "$vcfMultiAnnoMod" and -s "$vcfOut") {
      my $cmd = "rm -rf $vcfOut";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$vcfMultiAnnoMod" and -s "$vcfOutSorted\.avinput") {
      my $cmd = "rm -rf $vcfOutSorted\.avinput";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$vcfMultiAnnoMod" and -s "$vcfOutSorted\.invalid_input") {
      my $cmd = "rm -rf $vcfOutSorted\.invalid_input";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$vcfMultiAnnoMod" and -s "$vcfOutSorted\.refGene.invalid_input") {
      my $cmd = "rm -rf $vcfOutSorted\.refGene.invalid_input";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$vcfMultiAnnoMod" and -s "$vcfMultiAnno") {
      my $cmd = "rm -rf $vcfMultiAnno $vcfMultiAnnoVCF";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

  }

 ANNOVAR:

  if (exists($runTask{'annovar'}) and $options{'customCalling'} ne 'SRP' and -s "$options{'customCalling'}") { #do the annovar

    my $vcfOut = $options{'customCalling'};
    (my $vcfOutSorted = $vcfOut) =~ s/\.vcf$/.sorted.vcf/;
    my $vcfMultiAnno = $vcfOutSorted."\.$confs{'species'}_multianno.txt";
    my $vcfMultiAnnoVCF = $vcfOutSorted."\.$confs{'species'}_multianno.vcf";
    my $vcfMultiAnnoMod = $vcfOutSorted."\.$confs{'species'}_multianno.mod.vcf";
    my $vcfMultiAnnoModsnv = $vcfOutSorted."\.$confs{'species'}_multianno.mod.vcf.snv";
    my $vcfMultiAnnoModindel = $vcfOutSorted."\.$confs{'species'}_multianno.mod.vcf.indel";

    if (exists($runTask{'mergeCustomCallingChr'})) {
      goto MERGECUSTOM;
    }

    unless (-s "$vcfOut") {
      die "custome calling file $vcfOut does not exist!\n";
    }

    #annoVar annotate---------------------------------------------------------------------
    if ((-s "$vcfOut" or -s "$vcfMultiAnnoMod") and !-s "$vcfMultiAnnoModsnv") {

      unless (-s "$vcfMultiAnnoMod") {

        my $cmd = snvCalling->vcfSort($confs{'vcfSortBin'}, $vcfOut, $vcfOutSorted); #sort vcf
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

        $cmd = snvCalling->runAnnovar("$confs{'ANNOVARDIR'}/table_annovar.pl", $vcfOutSorted, $confs{'ANNOVARDB'}, $confs{'species'}); #table annovar
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

        $cmd = snvCalling->convertVCFannovar("$options{'bin'}/convert_annovar_vcf.pl", $vcfMultiAnno, $vcfMultiAnnoVCF);
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

      }

      my $cmd = snvCalling->grepSNVvcf($vcfMultiAnnoMod, $vcfMultiAnnoModsnv, $options{'snvPattern'});
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

      $cmd = snvCalling->grepINDELvcf($vcfMultiAnnoMod, $vcfMultiAnnoModindel, $options{'indelPattern'});
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});

    }

    #rm temporary files
    if (-s "$vcfMultiAnnoModsnv" and -s "$vcfOut") {
      my $cmd = "rm -rf $vcfOut";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$vcfMultiAnnoModsnv" and -s "$vcfOutSorted\.avinput") {
      my $cmd = "rm -rf $vcfOutSorted\.avinput";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$vcfMultiAnnoModsnv" and -s "$vcfOutSorted\.invalid_input") {
      my $cmd = "rm -rf $vcfOutSorted\.invalid_input";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$vcfMultiAnnoModsnv" and -s "$vcfOutSorted\.refGene.invalid_input") {
      my $cmd = "rm -rf $vcfOutSorted\.refGene.invalid_input";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$vcfMultiAnnoModsnv" and -s "$vcfMultiAnno") {
      my $cmd = "rm -rf $vcfMultiAnno $vcfMultiAnnoVCF";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    if (-s "$vcfMultiAnnoModsnv" and -s "$vcfMultiAnnoMod") {
      my $cmd = "rm -rf $vcfMultiAnnoMod";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }

    #------------------------------------------------------------------------------------

  } #annovar

 MERGECUSTOM:

  if (exists($runTask{'mergeCustomCallingChr'})) {
    my $callingString = $options{'customCalling'};
    my @customChrSnvs = bsd_glob("$options{'lanepath'}/04_SNV/*$callingString*.mod.vcf.snv");
    my $customChrSnvs = join(',', @customChrSnvs);
    my @customChrIndels = bsd_glob("$options{'lanepath'}/04_SNV/*$callingString*.mod.vcf.indel");
    my $customChrIndels = join(',', @customChrIndels);
    (my $vcfMultiAnnoModsnv = $customChrSnvs[0]) =~ s/\.chr[0-9A-Z]+//;
    (my $vcfMultiAnnoModindel =  $customChrIndels[0]) =~ s/\.chr[0-9A-Z]+//;
    unless (-s "$vcfMultiAnnoModsnv") {
      my $cmd = "perl $options{'bin'}/mergeMutFiles.pl $customChrSnvs >$vcfMultiAnnoModsnv";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      if (-s "$vcfMultiAnnoModsnv" and -s "$customChrSnvs[0]") {
        my $cmd = "rm ".join(" ", @customChrSnvs).' -rf';
        #RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    }
    unless (-s "$vcfMultiAnnoModindel") {
      my $cmd = "perl $options{'bin'}/mergeMutFiles.pl $customChrIndels >$vcfMultiAnnoModindel";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      if (-s "$vcfMultiAnnoModindel" and -s "$customChrIndels[0]") {
        my $cmd = "rm ".join(" ", @customChrIndels).' -rf';
        #RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    }
  }


 RECHECK:

  if ($options{'recheck'} ne 'SRP' and -s "$options{'recheck'}") { #do the recheck
    my $recheckBams = ($options{'recheckBams'} ne 'SRP')? $options{'recheckBams'} : $finalBam;
    my $recheckBasename = basename($options{'recheck'});
    my $recheckOut = "$options{'lanepath'}/04_SNV/$options{'sampleName'}\.$recheckBasename\.rechecked";
    my $cmd = snvCalling->rechecksnv("$options{'bin'}/novelSnvFilter_ACGT", $options{'recheck'}, $recheckBams, $recheckOut, $options{'chrPrefInBam'}, $options{'skipPileup'});
    if ($options{'recheck'} =~ /indel/) {
      $cmd = snvCalling->rechecksnv("$options{'bin'}/novelIndelFilter", $options{'recheck'}, $recheckBams, $recheckOut, $options{'chrPrefInBam'}, $options{'skipPileup'});
    }
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

 MERGE:

  if (exists($runTask{'mergeMutectChr'})) {
    my @mutectChrOuts = bsd_glob("$options{'lanepath'}/04_SNV/*.mutect");
    my $mutectChrOuts = join(',', @mutectChrOuts);
    my @mutectChrVcfs = bsd_glob("$options{'lanepath'}/04_SNV/*.mutect.genome.sorted.vcf.*_multianno.mod.vcf");
    my $mutectChrVcfs = join(',', @mutectChrVcfs);
    unless (-s "$muTectOut") {
      my $cmd = "perl $options{'bin'}/mergeMutFiles.pl $mutectChrOuts 2 >$muTectOut";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      if (-s "$muTectOut" and -s "$mutectChrOuts[0]") {
        my $cmd = "rm ".join(" ", @mutectChrOuts).' -rf';
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    }
    unless (-s "$vcfMultiAnnoMod") {
      my $cmd = "perl $options{'bin'}/mergeMutFiles.pl $mutectChrVcfs >$vcfMultiAnnoMod";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      if (-s "$vcfMultiAnnoMod" and -s "$mutectChrVcfs[0]") {
        my $cmd = "rm ".join(" ", @mutectChrVcfs).' -rf';
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    }
    if (-s "$vcfMultiAnnoMod" and -s "$muTectOut" and !-s "$muTectOut\.gz") {
      my $cmd = "gzip $muTectOut";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";

}


###
###runlevel5: CNA calling
###

$runlevels = 5;
if (exists $runlevel{$runlevels}) {

  unless (-e "$options{'lanepath'}/05_CNA") {
    my $cmd = "mkdir -p $options{'lanepath'}/05_CNA";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";
  print STDERR "transtate: $options{'transtate'} and tranclone: $options{'tranclone'}\n";

  my $tumorTitan = (-s "$options{'lanepath'}/04_SNV/$options{'sampleName'}\_titan")? "$options{'lanepath'}/04_SNV/$options{'sampleName'}\_titan" : die("$options{'sampleName'}\_titan not found!!!\n");
  my $tumorWig = (-s "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.wig")? "$options{'lanepath'}/03_STATS/$options{'sampleName'}\.wig" : die("$options{'sampleName'}\.wig not found!!!\n");
  my $normalWig;
  if ($options{'somaticInfo'} eq "SRP") {
    print STDERR "ERROR: somaticInfo is not provided! Must set for somatic calling!\n";
    exit 22;
  } elsif ( !exists( $somatic{$options{'sampleName'}} ) ){
    print STDERR "ERROR: $options{'sampleName'} is not in the somatic hash table!\n";
  } else { #get normal bam
    my $normalSampleName = $somatic{$options{'sampleName'}};
    my $aimedNormalWig = $normalSampleName.'.wig';
    if ($options{'nwigString'} ne 'SRP') {  #specific wiggle file
      $aimedNormalWig =~ s/\.wig/\.$options{'nwigString'}\.wig/;
    }
    $normalWig = (-s "$options{'root'}/$normalSampleName/03_STATS/$aimedNormalWig")? "$options{'root'}/$normalSampleName/03_STATS/$aimedNormalWig" : die("aimedNormalWig not found!!!\n");
  }

  if ($options{'seqType'} =~ /WGS/){
    $confs{'targetRegionTitan'} = "SRP";
  }

  my $segFile = "$options{'lanepath'}/05_CNA/$options{'sampleName'}\_nclones1.TitanCNA.segments.txt";
  unless (-s "$segFile"){
    my $cmd = cnaCalling->runTitan($confs{'RscriptBin'}, "$options{'bin'}/titan.R", "$options{'lanepath'}/05_CNA/", "$options{'bin'}", $options{'sampleName'}, $options{'gender'}, $tumorTitan, $tumorWig, $normalWig, $confs{'gcWigTitan'}, $confs{'mapWigTitan'},
                                   $options{'plpTitan'}, $options{'plpeTitan'}, $options{'ncTitan'}, $options{'ncmTitan'}, $options{'symmetric'}, $options{'transtate'}, $options{'tranclone'}, $confs{'targetRegionTitan'});
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  my $pngFileExample = "$options{'lanepath'}/05_CNA/$options{'sampleName'}\_nclones1_chr1.TitanCNA.png";
  unless (-s "$pngFileExample") {
    my $cmd = "perl $options{'bin'}/convertTitanPDF2PNG.pl $options{'lanepath'}/05_CNA/ nclones1";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  my $pngFileExample2 = "$options{'lanepath'}/05_CNA/$options{'sampleName'}\_nclones2_chr1.TitanCNA.png";
  unless (-s "$pngFileExample2") {
    my $cmd = "perl $options{'bin'}/convertTitanPDF2PNG.pl $options{'lanepath'}/05_CNA/ nclones2";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  my $pngFileExample3 = "$options{'lanepath'}/05_CNA/$options{'sampleName'}\_nclones3_chr1.TitanCNA.png";
  unless (-s "$pngFileExample3") {
    my $cmd = "perl $options{'bin'}/convertTitanPDF2PNG.pl $options{'lanepath'}/05_CNA/ nclones3";
    RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";

}


###
###runlevel8: merge mutations to a table for further recheck
###

$runlevels = 8;
if (exists($runlevel{$runlevels}) or exists($runTask{'mergeMutect'}) or exists($runTask{'mergeSamtools'}) or exists($runTask{'mergeStrelka'})) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  my $vcflist_mutect = "$options{'root'}/mutect.vcf.list";
  my $vcftable_mutect = "$options{'root'}/mutect.snv.table";
  my $originaltable_mutect = "$options{'root'}/mutect.snv.table.annotated";

  my $vcflist_samtools = "$options{'root'}/samtools.vcf.list";
  my $vcftable_samtools = "$options{'root'}/samtools.snv.table";
  my $originaltable_samtools = "$options{'root'}/samtools.snv.table.annotated";

  my $vcflist_strelka = "$options{'root'}/strelka.vcf.list";
  my $vcftable_strelka = "$options{'root'}/strelka.indel.table";
  my $originaltable_strelka = "$options{'root'}/strelka.indel.table.annotated";

  my $PREF;
  my $BLOOD;
  for my $eatumor (keys %somatic) {
    $PREF .= $eatumor.',';
  }
  for my $eanormal (keys %germline) {
    $PREF .= $eanormal.',';
    $BLOOD .= $eanormal.',';
  }
  $PREF =~ s/\,$//;
  $BLOOD =~ s/\,$//;

  print STDERR "PREF: $PREF\n";
  print STDERR "BLOOD: $BLOOD\n";

  #mutect merge
  if ((exists($runlevel{$runlevels}) or exists($runTask{'mergeMutectOnly'})) and !exists($runTask{'mergeSamtoolsOnly'}) and !exists($runTask{'mergeStrelkaOnly'})) {
    unless (-s "$vcflist_mutect") {
      for my $eatumor (keys %somatic) {
        my $eavcfmutect = "$options{'root'}/$eatumor/04_SNV/$eatumor\.mutect.genome.sorted.vcf.$confs{'species'}_multianno.mod.vcf";
        if (-s "$eavcfmutect") {
          my $cmd = "echo $eavcfmutect >>$vcflist_mutect";
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        } elsif (-s "$eavcfmutect".".gz") {
	  my $cmd = "echo ".$eavcfmutect.".gz"." >>$vcflist_mutect";
	  RunCommand($cmd,$options{'noexecute'},$options{'quiet'});    
        } else {
          print STDERR "warning: $eavcfmutect is not found!\n";
        }
      }
    }
    unless (-s "$originaltable_mutect") {
      unless (-s "$vcftable_mutect") {
        my $optionNosegdup = ($options{'mergeNonsegdup'} == 1)? '--nonsegdup':'';
        my $optionRare = ($options{'mergeRare'} == 1)? 'rare,muTect':'muTect';

        my $cmd = "perl $options{'bin'}/mergeMut.pl --list $vcflist_mutect --prefix $PREF --normal $BLOOD --type snv --task $optionRare --dbsnp yes $optionNosegdup >$vcftable_mutect";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
      if (-s "$vcftable_mutect") {
        my $cmd = "perl $options{'bin'}/junkAnnotate.pl --nonrepeat $confs{'repeatMasker'} --nonselfchain $confs{'selfChain'} --file $vcftable_mutect >$originaltable_mutect";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    }
  }

  #samtools merge
  if ((exists($runlevel{$runlevels}) or exists($runTask{'mergeSamtoolsOnly'})) and !exists($runTask{'mergeMutectOnly'}) and !exists($runTask{'mergeStrelkaOnly'})) {
    unless (-s "$vcflist_samtools") {
      for my $eatumor (keys %somatic) {
        my $eavcfsamtools = "$options{'root'}/$eatumor/04_SNV/$eatumor\.samtools.genome.sorted.vcf.$confs{'species'}_multianno.mod.vcf.snv";
        if (-s "$eavcfsamtools") {
          my $cmd = "echo $eavcfsamtools >>$vcflist_samtools";
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        } elsif (-s "$eavcfsamtools".".gz") {
          my $cmd = "echo ".$eavcfsamtools.".gz"." >>$vcflist_samtools";
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        } else {
          print STDERR "warning: $eavcfsamtools is not found!\n";
        }
      }
    }
    unless (-s "$originaltable_samtools") {
      unless (-s "$vcftable_samtools") {
        my $optionTask = ( $options{'rareVariants'} )? 'rare,titan':'titan';
        my $cmd = "perl $options{'bin'}/mergeMut.pl --list $vcflist_samtools --prefix $PREF --normal $BLOOD --type snv --task $optionTask --dbsnp yes --qualTitan $options{'qualTitan'} --vafTitan $options{'vafTitan'} --nonsegdup >$vcftable_samtools";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
      if (-s "$vcftable_samtools") {
        my $cmd = "perl $options{'bin'}/junkAnnotate.pl --nonrepeat $confs{'repeatMasker'} --nonselfchain $confs{'selfChain'} --file $vcftable_samtools >$originaltable_samtools";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    }
  }

  #strelka merge
  if ((exists($runlevel{$runlevels}) or exists($runTask{'mergeStrelkaOnly'})) and !exists($runTask{'mergeMutectOnly'}) and !exists($runTask{'mergeSamtoolsOnly'})) {
    unless (-s "$vcflist_strelka") {
      for my $eatumor (keys %somatic) {
        my $eavcfstrelka = "$options{'root'}/$eatumor/04_SNV/$eatumor\.strelka.indel.genome.sorted.vcf.$confs{'species'}_multianno.mod.vcf";
        if (-s "$eavcfstrelka") {
          my $cmd = "echo $eavcfstrelka >>$vcflist_strelka";
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        } elsif (-s "$eavcfstrelka".".gz") {
          my $cmd = "echo ".$eavcfstrelka.".gz"." >>$vcflist_strelka";
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        } else {
          print STDERR "warning: $eavcfstrelka is not found!\n";
        }
      }
    }
    unless (-s "$originaltable_strelka") {
      unless (-s "$vcftable_strelka") {
        my $optionTask = ( $options{'rareVariants'} )? 'rare,strelka':'strelka';
        my $cmd = "perl $options{'bin'}/mergeMut.pl --list $vcflist_strelka --prefix $PREF --normal $BLOOD --type indel --task $optionTask --dbsnp yes --nonsegdup >$vcftable_strelka";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
      if (-s "$vcftable_strelka") {
        my $cmd = "perl $options{'bin'}/junkAnnotate.pl --nonrepeat $confs{'repeatMasker'} --nonselfchain $confs{'selfChain'} --file $vcftable_strelka >$originaltable_strelka";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
    }
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";

}


###
###runlevel9: variant classification from merged and recheck files
###

$runlevels = 9;
if (exists($runlevel{$runlevels}) or exists($runTask{'MutectCallOnly'}) or exists($runTask{'SamtoolsCallOnly'})) {

  printtime();
  print STDERR "####### runlevel $runlevels now #######\n\n";

  my $originaltable_mutect = "$options{'root'}/mutect.snv.table.annotated";
  my $rechecklist_mutect = "$options{'root'}/mutect.snv.rechecked.list";
  my $realmaf_mutect = "$options{'root'}/mutect.snv.realmaf";
  my $varout_mutect = "$options{'root'}/mutect.snv.res";

  my $originaltable_samtools = "$options{'root'}/samtools.snv.table.annotated";
  my $rechecklist_samtools = "$options{'root'}/samtools.snv.rechecked.list";
  my $realmaf_samtools = "$options{'root'}/samtools.snv.realmaf";
  my $varout_samtools = "$options{'root'}/samtools.snv.res";

  my $PREF;
  my $BLOOD;
  for my $eatumor (keys %somatic) {
    $PREF .= $eatumor.',';
  }
  for my $eanormal (keys %germline) {
    $PREF .= $eanormal.',';
    $BLOOD .= $eanormal.',';
  }
  $PREF =~ s/\,$//;
  $BLOOD =~ s/\,$//;

  print STDERR "PREF: $PREF\n";
  print STDERR "BLOOD: $BLOOD\n";

  #mutect classification
  if ((exists($runlevel{$runlevels}) or exists($runTask{'MutectCallOnly'})) and !exists($runTask{'SamtoolsCallOnly'})) {
    unless (-s "$rechecklist_mutect") {   #generate recheck list for mutect
      for my $eatumor (keys %somatic) {
        my $eaRecheckmutect = "$options{'root'}/$eatumor/04_SNV/$eatumor\.mutect.snv.table.annotated.rechecked";
        if (-s "$eaRecheckmutect") {
          my $cmd = "echo $eaRecheckmutect >>$rechecklist_mutect";
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        } else {
          print STDERR "warning: $eaRecheckmutect is not found!\n";
        }
      }
      for my $eacontrol (keys %germline) {
        my $eaRecheckmutect = "$options{'root'}/$eacontrol/04_SNV/$eacontrol\.mutect.snv.table.annotated.rechecked";
        if (-s "$eaRecheckmutect") {
          my $cmd = "echo $eaRecheckmutect >>$rechecklist_mutect";
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        } else {
          print STDERR "warning: $eaRecheckmutect is not found!\n";
        }
      }
    }
    unless (-s "$realmaf_mutect") {
      my $cmd = "perl $options{'bin'}/realmaf.pl --file $rechecklist_mutect --type snv --original $originaltable_mutect --prefix $PREF --blood $BLOOD >$realmaf_mutect";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    unless (-s "$varout_mutect") {
      my $lastcolindex = `perl $options{'bin'}/columnIndex.pl cmedianav $realmaf_mutect`;
      $lastcolindex =~ s/\n$//;
      my $cmd = "perl $options{'bin'}/intersectFiles.pl -o $originaltable_mutect -m $realmaf_mutect -extraOM 4,4 -column 5-$lastcolindex >$varout_mutect";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      $cmd = "perl $options{'bin'}/columnRearrange.pl --file $varout_mutect --prefix $PREF --sep bp >$varout_mutect\.1";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      $cmd = "mv $varout_mutect\.1 $varout_mutect -f";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    unless (-s "$varout_mutect\.filtered") {
      my $cmd = "perl $options{'bin'}/recurrency.pl --file $varout_mutect --type snv --task filter --noStrandBias $options{'noStrandBias'} --cmeme $options{'cmeme'} --cmedian $options{'cmedian'} >$varout_mutect\.filtered";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    unless (-s "$varout_mutect\.filtered\.classified") {
      my $cmd = "perl $options{'bin'}/recurrency.pl --file $varout_mutect\.filtered --type snv --task somatic --somaticInfo $options{'somaticInfo'} --noStrandBias $options{'noStrandBias'} --cmeme $options{'cmeme'} --cmedian $options{'cmedian'} >$varout_mutect\.filtered\.classified";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    unless (-s "$varout_mutect\.filtered\.classified\.founds") {
      my $cmd = "perl $options{'bin'}/recurrency.pl --file $varout_mutect\.filtered\.classified --type snv --task samfounds --somaticInfo $options{'somaticInfo'} --noStrandBias $options{'noStrandBias'} --cmeme $options{'cmeme'} --cmedian $options{'cmedian'} >$varout_mutect\.filtered\.classified\.founds";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    unless (-s "$varout_mutect\.filtered\.classified\.founds\.nopara") {
      my $cmd = "perl $options{'bin'}/readsFlankingVariants.pl $confs{'GFASTA'} $varout_mutect\.filtered\.classified\.founds snv >$varout_mutect\.filtered\.classified\.founds.flanking.fa";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      unless (-s "$varout_mutect\.filtered\.classified\.founds.flanking.bam") {
         $cmd = bwaMapping->bowtieMappingSnv($confs{'bowtieBin'}, $confs{'BowtieINDEX'}, "$varout_mutect\.filtered\.classified\.founds.flanking.fa", "$varout_mutect\.filtered\.classified\.founds.flanking.sam", $options{'threads'});
         RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
         $cmd = bwaMapping->samToBam($confs{'samtoolsBin'}, "$varout_mutect\.filtered\.classified\.founds.flanking.sam", "$varout_mutect\.filtered\.classified\.founds.flanking.bam", $options{'threads'});
         RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
      $cmd = "$options{'bin'}/mappingFlankingVariants --mapping $varout_mutect\.filtered\.classified\.founds.flanking.bam --readlength $options{'readlen'} --type s >$varout_mutect\.filtered\.classified\.founds.flanking.bam.out";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      $cmd = "perl $options{'bin'}/badvariantmapping.pl $varout_mutect\.filtered\.classified\.founds.flanking.bam.out >$varout_mutect\.filtered\.classified\.founds.flanking.bam.out.bad";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      $cmd = "perl $options{'bin'}/intersectFiles.pl -o $varout_mutect\.filtered\.classified\.founds -m $varout_mutect\.filtered\.classified\.founds.flanking.bam.out.bad -count >$varout_mutect\.filtered\.classified\.founds\.1";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      my $PRC = `head -1 $varout_mutect\.filtered\.classified\.founds\.1 |awk '{print \$NF}'`;
      $PRC =~ s/\n$//;
      my $PRCI = `perl $options{'bin'}/columnIndex.pl $PRC $varout_mutect\.filtered\.classified\.founds\.1`;
      $PRCI =~ s/\n$//;
      $PRCI += 1;
      $cmd = "awk -F\"\\t\" \'\$$PRCI \!\= 1\' $varout_mutect\.filtered\.classified\.founds\.1 >$varout_mutect\.filtered\.classified\.founds\.nopara";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    unless (-s "$varout_mutect\.filtered\.classified\.founds\.nopara\.somatic") {
        my $SOMI = `perl $options{'bin'}/columnIndex.pl somatic $varout_mutect\.filtered\.classified\.founds\.1`;
        $SOMI =~ s/\n$//;
        $SOMI += 1;
        my $cmd = "awk -F\"\\t\" \'\$$SOMI \!\= \"NA\"\' $varout_mutect\.filtered\.classified\.founds\.nopara >$varout_mutect\.filtered\.classified\.founds\.nopara\.somatic";
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    unless (-s "$varout_mutect\.filtered\.classified\.founds\.nopara\.somatic.table") {
      my $cmd = "perl $options{'bin'}/mutationTable.pl --mutation $varout_mutect\.filtered\.classified\.founds\.nopara\.somatic --type snv --normal $BLOOD --prefix $PREF >$varout_mutect\.filtered\.classified\.founds\.nopara\.somatic.table";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    unless (-s "$varout_mutect\.filtered\.classified\.founds\.nopara\.somatic.table.simplified") {
      my $cmd = "perl $options{'bin'}/recurrency.pl --file $varout_mutect\.filtered\.classified\.founds\.nopara\.somatic.table --type snv --task split --somaticInfo $options{'somaticInfo'} --noStrandBias $options{'noStrandBias'} --cmeme $options{'cmeme'} --cmedian $options{'cmedian'} >$varout_mutect\.filtered\.classified\.founds\.nopara\.somatic.table.simplified";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
  }


  #samtools merge
  if ((exists($runlevel{$runlevels}) or exists($runTask{'SamtoolsCallOnly'})) and !exists($runTask{'MutectCallOnly'})) {
    unless (-s "$rechecklist_samtools") {   #generate recheck list for samtools
      for my $eatumor (keys %somatic) {
        my $eaRechecksamtools = "$options{'root'}/$eatumor/04_SNV/$eatumor\.samtools.snv.table.annotated.rechecked";
        if (-s "$eaRechecksamtools") {
          my $cmd = "echo $eaRechecksamtools >>$rechecklist_samtools";
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        } else {
          print STDERR "warning: $eaRechecksamtools is not found!\n";
        }
      }
      for my $eacontrol (keys %germline) {
        my $eaRechecksamtools = "$options{'root'}/$eacontrol/04_SNV/$eacontrol\.samtools.snv.table.annotated.rechecked";
        if (-s "$eaRechecksamtools") {
          my $cmd = "echo $eaRechecksamtools >>$rechecklist_samtools";
          RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
        } else {
          print STDERR "warning: $eaRechecksamtools is not found!\n";
        }
      }
    }
    unless (-s "$realmaf_samtools") {
      my $cmd = "perl $options{'bin'}/realmaf.pl --file $rechecklist_samtools --type snv --original $originaltable_samtools --prefix $PREF --blood $BLOOD >$realmaf_samtools";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    unless (-s "$varout_samtools") {
      my $lastcolindex = `perl $options{'bin'}/columnIndex.pl cmedianav $realmaf_samtools`;
      $lastcolindex =~ s/\n$//;
      my $cmd = "perl $options{'bin'}/intersectFiles.pl -o $originaltable_samtools -m $realmaf_samtools -extraOM 4,4 -column 5-$lastcolindex >$varout_samtools";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      $cmd = "perl $options{'bin'}/columnRearrange.pl --file $varout_samtools --prefix $PREF --sep bp >$varout_samtools\.1";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      $cmd = "mv $varout_samtools\.1 $varout_samtools -f";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    unless (-s "$varout_samtools\.filtered") {
      my $cmd = "perl $options{'bin'}/recurrency.pl --file $varout_samtools --type snv --task filter >$varout_samtools\.filtered";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    unless (-s "$varout_samtools\.filtered\.nopara") {
      my $cmd = "perl $options{'bin'}/readsFlankingVariants.pl $confs{'GFASTA'} $varout_samtools\.filtered snv >$varout_samtools\.filtered.flanking.fa";
      unless (-s "$varout_samtools\.filtered.flanking.fa") {
        RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
      unless (-s "$varout_samtools\.filtered.flanking.bam") {
	 unless (-s "$varout_samtools\.filtered.flanking.sam") {
           $cmd = bwaMapping->bowtieMappingSnv($confs{'bowtieBin'}, $confs{'BowtieINDEX'}, "$varout_samtools\.filtered.flanking.fa", "$varout_samtools\.filtered.flanking.sam", $options{'threads'});
           RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
	 }
         $cmd = bwaMapping->samToBam($confs{'samtoolsBin'}, "$varout_samtools\.filtered.flanking.sam", "$varout_samtools\.filtered.flanking.bam", $options{'threads'});
         RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      }
      $cmd = "$options{'bin'}/mappingFlankingVariants --mapping $varout_samtools\.filtered.flanking.bam --readlength $options{'readlen'} --type s >$varout_samtools\.filtered.flanking.bam.out";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      $cmd = "perl $options{'bin'}/badvariantmapping.pl $varout_samtools\.filtered.flanking.bam.out >$varout_samtools\.filtered.flanking.bam.out.bad";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      $cmd = "perl $options{'bin'}/intersectFiles.pl -o $varout_samtools\.filtered -m $varout_samtools\.filtered.flanking.bam.out.bad -count >$varout_samtools\.filtered\.1";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      my $PRC = `head -1 $varout_samtools\.filtered\.1 |awk '{print \$NF}'`;
      $PRC =~ s/\n$//;
      my $PRCI = `perl $options{'bin'}/columnIndex.pl $PRC $varout_samtools\.filtered\.1`;
      $PRCI =~ s/\n$//;
      $PRCI += 1;
      $cmd = "awk -F\"\\t\" \'\$$PRCI \!\= 1\' $varout_samtools\.filtered\.1 >$varout_samtools\.filtered\.nopara";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
    unless (-e "$options{'root'}/titan") {
      my $cmd = "perl $options{'bin'}/titanCNAprepare.pl $varout_samtools\.filtered\.nopara $options{'somaticInfo'} 1 $options{'homoThred'} $options{'ndepthThred'} $options{'germlineLOH'}";
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
      $cmd = "perl $options{'bin'}/redistributeTitan.pl $options{'root'}/titan/ $options{'root'}";                #distribute allele counts for titan
      RunCommand($cmd,$options{'noexecute'},$options{'quiet'});
    }
  }

  printtime();
  print STDERR "####### runlevel $runlevels done #######\n\n";

}




###------------###################################################################################################
##  sub-region  #################################################################################################
###------------###################################################################################################

sub RunCommand {
  my ($command,$noexecute,$quiet) = @_ ;
  unless ($quiet){
    printtime();
    print STDERR "$command\n\n";
  }
  unless ($noexecute) {
    system($command);
  }
}

sub helpm {
  print STDERR "\nGENERAL OPTIONS:\n\t--runlevel\tthe steps of runlevel, from 1-3, either rl1-rl2 or rl. See below for options for each runlevel.\n\tor\n";
  print STDERR "\t--runTask\tthe specific task. e.g., \'QC\', \'indelRealignment\', \'MarkDuplicates\', \'recalMD\' and \'BaseRecalibration\'.\n\n";
  print STDERR "\t--skipTask\tthe task to be skipped, e.g, \'BaseRecalibration\' for reduced bam processing time while not affecting the calling significantly.\n";
  print STDERR "\t--configure\tthe tab delimited file containing conf info for annotations\n";
  print STDERR "\t--sampleName\tthe name of the lane needed to be processed (must set for runlevel 1-5)\n";
  print STDERR "\t--gender\tthe gender of the subject (default: female), set it to male to correct CNA calls for X chromosome\n";
  print STDERR "\t--seqType\tcomma separated, possible arguments \'paired-end\', \'single-end\', \'WXS\' and \'WGS\' (default).\n";
  print STDERR "\t--root\t\tthe root directory of the pipeline (default is \$bin/../PIPELINE/, MUST set using other dir)\n";
  print STDERR "\t--species\tspecify the reference version of the species, such as hg19 (default), mm10.\n";
  print STDERR "\t--patient\tthe patient id, which will be written into the target file for edgeR\n";
  print STDERR "\t--tissue\tthe tissue type name (like \'normal\', \'cancer\'), for the target file for running edgeR and cuffdiff\n";
  print STDERR "\t--Rbinary\tthe name of R executable, default is \'R\'. Set if your R binary name is different.\n";
  print STDERR "\t--platform\tthe sequencing platform (default: ILLUMINA)\n\n";

  #print STDERR "CONTROL OPTIONS FOR EACH RUNLEVEL:\n";
  #print STDERR "runlevel 1: quality checking and insert size estimatiion using part of reads\n";
  print STDERR "\t--readpool\tthe directory where all the read files with names ending with \.f(ast)?q\.[gb]z2\? located.\n";
  print STDERR "\t--qcOFF\t\tturn off quality check.\n";
  print STDERR "\t--FASTQ1\tcomma separated zipped fastq file names for mate 1 (without dir name).\n";
  print STDERR "\t--FASTQ2\tcomma separated zipped fastq file names for mate 2 (without dir name).\n";
  print STDERR "\t--fastqFiles1\tcomma separated zipped fastq file names for mate 1 (with dir name).\n";
  print STDERR "\t--fastqFiles2\tcomma separated zipped fastq file names for mate 2 (with dir name).\n";
  print STDERR "\t--readlen\tthe sequenced read length (default the length of the first read in the fastq file)\n";
  print STDERR "\t--bamID\t\tthe ID for read group record in bam file, default is 1. set to the id needed to differentiate seq-experiments\n";
  print STDERR "\t--splitChr\tsplit Chromosomes. (default not set)\n";

  print STDERR "\nrunlevel 2: mapping and report of mapping statistics\n";
  print STDERR "\t--mapper\tthe mapper for read-alignment, now support \'bwa\' (default).\n";
  print STDERR "\t--gf\t\tthe graphical format in mapping report, \'png\' (default) or \'pdf\' (when a x11 window is not available)\n";

  print STDERR "\nrunlevel 3: STATS\n";

  print STDERR "\nrunlevel 4: SNV/INDEL calling (muTect for somatic, samtools for germline, stelka for indel)\n";
  print STDERR "\t--somaticInfo\tsample information for tumor and normal pair (tab delimited)\n";
  print STDERR "\t--germline\tgermline caller name, specify it to samtools if wanted\n";
  print STDERR "\t--samCallmaxDepth\tthe maximum depth for samtools mpileup to work (default 400, set to the high quantile depth)\n";
  print STDERR "\t--indel\t\tindel caller name, specify it to strelka if wanted\n";
  print STDERR "\t--recheck\trecheck bam files against a tsv file contains mutations\n";
  print STDERR "\t--recheckBams\trecheck bam file or list of files if not the one under 02_MAPPING\n";
  print STDERR "\t--chrPrefInBam\tthe prefix of chromosome names in the bam file. default is no prefix, set to the actual prefix you have in the bam.\n";

  print STDERR "\nrunlevel 5: CNA calling (TitanCNA)\n";
  print STDERR "\t--plpTitan\tploidy starting value, default \'2.0\'\n";
  print STDERR "\t--plpeTitan\twhether estimate ploidy by itself, default \'TRUE\'\n";
  print STDERR "\t--ncTitan\tnormal contamination initial value, default \'0.5\'\n";
  print STDERR "\t--ncmTitan\tmethod used for estimating normal contamination, default \'map\', use \'fixed\' to fix.\n";

  print STDERR "\nrunlevel 8: merge calls\n";
  print STDERR "\nrunlevel 9: variant classification\n";
  print STDERR "\t--homoThred\tThreshold for excluding homozygous sites in the normal, default \'0.85\'\n";
  print STDERR "\t--ndepthThred\tAt least this depth in the normal for being considered for titan, default \'8\'\n";

  print STDERR "\nOTHER OPTIONS\n";
  print STDERR "\t--noexecute\tdo not execute the command, for testing purpose\n";
  print STDERR "\t--quiet\t\tdo not print the command line calls and time information\n";
  print STDERR "\t--threads\tthe number of threads used for the mapping (default 1)\n";
  print STDERR "\t--help\t\tprint this help message\n";
  print STDERR "\t--tmpDir\ttmp dir for generating large tmp files\n";
  print STDERR "\t--bzip\t\tthe read fastq files are bziped, rather than gziped (default).\n";
  print STDERR "\t--bin\t\tbinary directories (default is where this script is executed.').\n\n";

  #print STDERR "\nSynopsis: RTrace.pl --runlevel 1 --sampleName <sample1> --runID <ID> --root <dir_root> --anno <dir_anno> 2>>run.log\n";
  #print STDERR "Synopsis: RTrace.pl --runlevel 2 --sampleName <sample1> --runID <ID> --root <dir_root> --anno <dir_anno> --patient <ID> --tissue <type> --threads <N> 2>>run.log\n";
  #print STDERR "Synopsis: RTrace.pl --runlevel 3-4 --sampleName <sample1> --runID <ID> --root <dir_root> --anno <dir_anno> --RA 1 --threads <N> 2>>run.log\n";
  #print STDERR "Synopsis: RTrace.pl --runlevel 5 --sampleName <sample1> --runID <ID> --root <dir_root> --anno <dir_anno> --threads <N> 2>>run.log\n";
  #print STDERR "Synopsis: RTrace.pl --runlevel 7 --root <dir_root> --anno <dir_anno> --priordf 1 2>>run.log\n";
  #print STDERR "remember to set --species option to a string, such as mm10, if it is not a human sample!!!\n\n";

  #print STDERR "Runlevel dependencies (->): 4->3->2->1, 6->5->2->1, 7->2->1\n\n";
  exit 0;
}

sub printtime {
  my @time = localtime(time);
  printf STDERR "\n[".($time[5]+1900)."\/".($time[4]+1)."\/".$time[3]." ".$time[2].":".$time[1].":".$time[0]."]\t";
}

sub round {
    my $number = shift;
    my $tmp = int($number);
    if ($number >= ($tmp+0.5)){
      $tmp++;
    }
    return $tmp;
}

sub uniqueArray {
   my $array = shift;
   my %arraytmp;
   foreach my $item (@{$array}){
     $arraytmp{$item} = '';
   }
   my @arraytmp = keys %arraytmp;
   return @arraytmp;
}

sub ceiling {
  my ($num) = @_;
  return int($num) + ($num > int($num));
}


sub partitionArray {
    my ($arr, $N) = @_;

    my @res;
    my $i = 0;

    while ($i + $N-1 <= $#$arr) {
        push @res, [@$arr[$i .. $i+$N-1]];
        $i += $N;
    }

    if ($i <= $#$arr) {
        push @res, [@$arr[$i .. $#$arr]];
    }
    return \@res;
}
