#!/bin/bash -l

# Set job time to 2 day.
#SBATCH --time=2-00:00:00

# Set a name for the job, visible in `squeue`
#SBATCH --job-name="emapping"

# One node.
#SBATCH --nodes=1

# One task (cpu cores in MSI)
#SBATCH --ntasks=4

# One CPU/core per task, threads?
## SBATCH --cpus-per-task=1

# 4GB of RAM (pmem is not recognized)
#SBATCH --mem=32gb

# Who to send mail to.
#SBATCH --mail-user=ruping@umn.edu

# What type of mail to send
#SBATCH --mail-type=FAIL

#workdir
##SBATCH --workdir=/home/ruping/ruping/MRS/crc_uchi/

#track job status
#SBATCH --profile=task

#queue name (use small for general queue)
##SBATCH -p sunpath


module load perl
module load gcc/5.4.0
module load R/4.0.0
module load samtools/1.5
module load gmap/2015-09-29
module load bwa/0.7.12
module load bowtie2/2.3.4.1
module load bedtools/2.29.2
module load java/jdk1.8.0_144


export TOOLP=/home/ruping/shared/tools/
export DANNO=/home/ruping/shared/annotation/
export PATH=$PATH:$TOOLP/trick/:$TOOLP/vcftools/current/bin/
export PERL5LIB=${PERL5LIB}:/home/ruping/shared/libs/perl5/


CONFIG=$1
ROOT=$2
SOMATICINFO=$3
SAMPLE=$4
BAM=$5
THREADS=4
READLEN=100
JAVATMP="/home/ruping/shared/tmp/"

# this step will make raw somatic variant scan using muTect
perl $TOOLP/seqare/DTrace.pl --configure $CONFIG --javaTmp $JAVATMP --maxMem 6g --runlevel 2 --sampleName $SAMPLE --bams $BAM --seqType paired-end,WXS --qcOFF --root $ROOT --threads $THREADS --somaticInfo $SOMATICINFO --skipTask indelRealignment,BaseRecalibration,recalMD --runTask remap 2>>$ROOT/$SAMPLE.run.log
#perl $TOOLP/seqare/DTrace.pl --configure $CONFIG --javaTmp /N/projects/curtis/POS/javaTmp/ --maxMem 4g --runlevel 2 --sampleName $SAMPLE --bams $BAM --seqType paired-end,WXS --root $ROOT --threads $THREADS --somaticInfo $SOMATICINFO --runTask recalMD 2>>$ROOT/$SAMPLE.run.log
#perl $TOOLP/seqare/DTrace.pl --configure $CONFIG --javaTmp /N/projects/curtis/POS/javaTmp/ --maxMem 4g --runlevel 2 --sampleName $SAMPLE --seqType paired-end,WGS,ignore --root $ROOT --threads $THREADS --somaticInfo $SOMATICINFO --skipTask recalMD 2>>$ROOT/$SAMPLE.run.log
#perl $TOOLP/seqare/DTrace.pl --configure $CONFIG --runlevel 3 --sampleName $SAMPLE --seqType paired-end,WGS --root $ROOT --threads $THREADS --somaticInfo $SOMATICINFO --keepBedCount --lorenzScaleFactor 0.604 --maxInsLine 10000000 --chrPrefInBam chr 2>>$ROOT/$SAMPLE.run.log
