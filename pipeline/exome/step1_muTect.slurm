#!/bin/bash -l

# Set job time to 2 day.
#SBATCH --time=1-12:00:00

# Set a name for the job, visible in `squeue`
#SBATCH --job-name="emuTect"

# One node.
#SBATCH --nodes=1

# One task (cpu cores in MSI)
#SBATCH --ntasks=1

# One CPU/core per task, threads?
## SBATCH --cpus-per-task=1

# 4GB of RAM (pmem is not recognized)
#SBATCH --mem=16gb

# Who to send mail to.
#SBATCH --mail-user=maese004@umn.edu

# What type of mail to send
#SBATCH --mail-type=FAIL

#workdir
##SBATCH --workdir=/home/ruping/ruping/MRS/crc_uchi/

#track job status
#SBATCH --profile=task

#queue name (use small for general queue)
##SBATCH -p sunpath

##loading modules in the super computing environment##
module load perl
module load gcc/5.4.0
module load samtools/1.5
module load gmap/2015-09-29
module load bwa/0.7.12
module load bowtie2/2.3.4.1
module load bedtools/2.29.2

## set environment variables ########
export TOOLP=/home/ruping/shared/tools/
export DANNO=/home/ruping/shared/annotation/
export PATH=$PATH:$TOOLP/trick/:$TOOLP/vcftools/current/bin/
export PERL5LIB=${PERL5LIB}:/home/ruping/shared/libs/perl5/


## input ####
CONFIG=$1
ROOT=$2
SOMATICINFO=$3
SAMPLE=$4
#CHR=$5
THREADS=1
JAVATMP="/home/ruping/shared/tmp/"

# this step will make raw somatic variant scan using muTect
perl $TOOLP/seqare/DTrace.pl --configure $CONFIG --javaTmp $JAVATMP --runlevel 4 --sampleName $SAMPLE --seqType paired-end,WXS --root $ROOT --threads $THREADS --somaticInfo $SOMATICINFO --samplePairNames redefine 2>>$ROOT/$SAMPLE.run.log.muTect
#perl $TOOLP/seqare/DTrace.pl --configure $CONFIG --runlevel 4 --sampleName $SAMPLE --runTask mergeMutectChr --seqType paired-end,WGS --root $ROOT --threads $THREADS --somaticInfo $SOMATICINFO 2>>$ROOT/$SAMPLE.run.log
