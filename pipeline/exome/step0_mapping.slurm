#!/bin/bash -l

# Set job time to 2 day. 
#SBATCH --time=2-00:00:00

# Set a name for the job, visible in `squeue`
#SBATCH --job-name="emapping"

# One node.
#SBATCH --nodes=1

# One task (cpu cores)
#SBATCH --ntasks=4

# One CPU/core per task, threads?
## SBATCH --cpus-per-task=1

# 4GB of RAM (pmem is not recognized)
#SBATCH --mem=32gb

# Who to send mail to. Change to your own email
#SBATCH --mail-user=ruping@umn.edu

# What type of mail to send
#SBATCH --mail-type=FAIL

#track job status
#SBATCH --profile=task

#queue name (use small for general queue)
##SBATCH -p sunpath


# module loading section START #
# the module loading section is specific to MSI cluster
# users need to make sure the following programs are installed
# and can be executed
module load perl
module load gcc/5.4.0
module load R/4.0.4
module load samtools/1.5
module load gmap/2015-09-29
module load bwa/0.7.12
module load bowtie2/2.3.4.1
module load bedtools/2.29.2
module load java/jdk1.8.0_144

#set the ITHVARIANT to your installation direction of ith.Variant
ITHVARIANT=/home/ruping/shared/tools/ith.Variant/
#set the annotation directory where HG38 related annotation files are located
DANNO=/home/ruping/shared/annotation/
#make the perl library available
export PERL5LIB=${PERL5LIB}:/home/ruping/shared/libs/perl5/
# module loading section END #


CONFIG=$1
ROOT=$2
SOMATICINFO=$3
READPOOL=$4
SAMPLE=$5
FASTQ1=$6
FASTQ2=$7
#modify this if you decide to change the number of threads used. this should match the --ntasks in the header of this script
THREADS=4
#modify this if you know the read length, which is different from the setting here.
READLEN=100
#modify the temporary folder where temporary files will be written to.
JAVATMP="/home/ruping/shared/tmp/"


# this step will map the sequencing reads of a particular sample against the reference genome.
#perl $TOOLP/seqare/DTrace.pl --configure $CONFIG --runlevel 1 --sampleName $SAMPLE --readpool $READPOOL --FASTQ1 $FASTQ1 --FASTQ2 $FASTQ2 --seqType paired-end,WXS --root $ROOT --threads $THREADS --somaticInfo $SOMATICINFO 2>>$ROOT/$SAMPLE.run.log
perl $TOOLP/seqare/DTrace.pl --configure $CONFIG --javaTmp $JAVATMP --maxMem 6g --runlevel 1-2 --sampleName $SAMPLE --readpool $READPOOL --FASTQ1 $FASTQ1 --FASTQ2 $FASTQ2 --seqType paired-end,WXS --qcOFF --root $ROOT --threads $THREADS --somaticInfo $SOMATICINFO --skipTask indelRealignment,BaseRecalibration,recalMD 2>>$ROOT/$SAMPLE.run.log
#perl $TOOLP/seqare/DTrace.pl --configure $CONFIG --javaTmp /N/projects/curtis/POS/javaTmp/ --maxMem 4g --runlevel 2 --sampleName $SAMPLE --bams $BAM --seqType paired-end,WXS --root $ROOT --threads $THREADS --somaticInfo $SOMATICINFO --runTask recalMD 2>>$ROOT/$SAMPLE.run.log
#perl $TOOLP/seqare/DTrace.pl --configure $CONFIG --javaTmp /N/projects/curtis/POS/javaTmp/ --maxMem 4g --runlevel 2 --sampleName $SAMPLE --seqType paired-end,WGS,ignore --root $ROOT --threads $THREADS --somaticInfo $SOMATICINFO --skipTask recalMD 2>>$ROOT/$SAMPLE.run.log
#perl $TOOLP/seqare/DTrace.pl --configure $CONFIG --runlevel 3 --sampleName $SAMPLE --seqType paired-end,WGS --root $ROOT --threads $THREADS --somaticInfo $SOMATICINFO --keepBedCount --lorenzScaleFactor 0.604 --maxInsLine 10000000 --chrPrefInBam chr 2>>$ROOT/$SAMPLE.run.log
