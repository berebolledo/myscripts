#!/bin/bash

#SBATCH -C new
#SBATCH --cpus-per-task=8
#SBATCH -J map3
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=berebolledo@gmail.com
#SBATCH -e slurm-%j.err

#source /galaxy/home/boris/.bash_profile
subject=${1}


###############################################################################

reference="/nfs/brubeck.bx.psu.edu/scratch1/boris/mtprojectreference/mtprojectref.fa"

f1=${subject}_1.fastq.gz
f2=${subject}_2.fastq.gz
mkdir -p mapped

bwa mem \
    -t 8 \
    -M   \
    -R "@RG\tID:${subject}\tLB:amelia\tPL:ILLUMINA\tSM:${subject}" \
    ${reference}      \
    ${f1}             \
    ${f2} |           \
    samtools view -Sb - > mapped/${subject}.mapped.bam

###############################################################################
