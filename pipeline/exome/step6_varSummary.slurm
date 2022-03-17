#!/bin/bash -l

# Set job time to 2 day.
#SBATCH --time=1-00:00:00

# Set a name for the job, visible in `squeue`
#SBATCH --job-name="evarcls"

# One node.
#SBATCH --nodes=1

# One task
#SBATCH --ntasks=2

# One CPU/core per task, threads?
## SBATCH --cpus-per-task=1

# 4GB of RAM
#SBATCH --mem=64gb

# Who to send mail to.
#SBATCH --mail-user=maese004@umn.edu

# What type of mail to send
#SBATCH --mail-type=FAIL

#workdir
##SBATCH --workdir=/home/ruping/ruping/SJOS/batch1/

#track job status
#SBATCH --profile=task

#queue name
#SBATCH -p sunpath


module load perl
module load R/3.6.0
module load samtools/1.5
module load bcftools/1.6
module load gmap/2015-09-29
module load bwa/0.7.12
module load bowtie2/2.3.4.1
module load bedtools/2.29.2
module load gcc/5.4.0

export TOOLP=/home/ruping/shared/tools/
export DANNO=/home/ruping/shared/annotation/
export PATH=$PATH:$TOOLP/trick/:$TOOLP/vcftools/current/bin/
export PERL5LIB=${PERL5LIB}:/home/ruping/shared/libs/perl5/

CONFIG=$1
ROOT=$2
SOMATICINFO=$3
XLOH=$4
THREADS=2
JAVATMP="/home/ruping/shared/tmp/"
READLEN=76

perl $TOOLP/seqare/DTrace.pl --configure $CONFIG --runlevel 9 --seqType paired-end,WXS --root $ROOT --threads $THREADS --somaticInfo $SOMATICINFO --readlen $READLEN --germlineLOH $XLOH 2>>$ROOT/classification.run.log
#--cmeme 6.5 --cmedian 3
#--runTask SamtoolsCallOnly or --runTask MutectCallOnly
#--germlineLOH $XLOH
