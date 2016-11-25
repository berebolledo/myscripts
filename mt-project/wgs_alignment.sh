#!/bin/bash

subject=${1}
genome=${2}
f1=${subject}_1.fastq.gz
f2=${subject}_2.fastq.gz
mkdir -p mapped

hg="/nfs/brubeck.bx.psu.edu/scratch1/boris/mtprojectreference/mtprojectref.fa"
mt="/nfs/brubeck.bx.psu.edu/scratch1/boris/mtprojectreference/rCRS.fasta"

if [ "${genome}" = "mt" ]; then
    reference="${mt}"
    outname="${subject}_ref-mt_chrM_mapped.bam"
else
    reference="${hg}"
    outname="${subject}_ref-hg_chrM_mapped.bam"
fi


bwa mem   \
    -t 16 \
    -M    \
    -R "@RG\tID:${subject}\tLB:amelia\tPL:ILLUMINA\tSM:${subject}"\
    ${reference}  \
    ${f1}         \
    ${f2}         \
    | samtools view -Sb - | bamtools filter -region chrM -out mapped/"${outname}"

    
