#!/bin/bash

#SBATCH -C new
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH -J het.pipe
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=berebolledo@gmail.com
#SBATCH -e %j.err

source /galaxy/home/boris/.bash_profile
reference="/nfs/brubeck.bx.psu.edu/scratch1/boris/mtprojectreference/mtprojectref.fa"

# SRR1310297_1.fastq.gz
# SRR1310297_2.fastq.gz

subject=${1}
f1="${subject}_1.fastq.gz"
f2="${subject}_2.fastq.gz"

bwa7 mem  \
    -t 64 \
    -M    \
    -R "@RG\tID:${subject}\tLB:WholeGenome\tPL:ILLUMINA\tSM:${subject}"\
     ${reference} \
     ${f1} \
     ${f2} 2>/dev/null|\
     samtools view -o ${subject}.mapped.bam -Sb -
