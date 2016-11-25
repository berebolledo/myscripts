#!/bin/bash

#SBATCH -C new
##SBATCH --ntasks=1
##SBATCH --cpus-per-task=1
#SBATCH -J kallisto
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=berebolledo@gmail.com
#SBATCH -e slurm-%j.err

source /galaxy/home/boris/.bash_profile

index=${1}
sampleid=${2}
f1=${sampleid}_1.fastq.gz
f2=${sampleid}_2.fastq.gz

kallisto quant      \
    --pseudobam     \
    -i ${index}     \
    -o ${sampleid}  \
    -b 100          \
    ${f1} ${f2} | samtools view  -o ${sampleid}.kallisto.bam -Sb -


