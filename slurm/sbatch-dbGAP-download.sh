#!/bin/bash

#SBATCH -C new
##SBATCH --ntasks=1
##SBATCH --cpus-per-task=1
#SBATCH -J dbGAP
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=berebolledo@gmail.com
#SBATCH -e slurm-%j.err

# Definitions

source /galaxy/home/boris/.bash_profile
ASPERA="/home/boris/.aspera/connect"
DOWNLOADS="/nfs/brubeck.bx.psu.edu/scratch5/boris/amelia/ncbi/sra"
SRA_ID=${1}

# Check for file already downloaded

if [ -s ${DOWNLOADS}/${SRA_ID}_1.fastq.gz ] && [ -s ${DOWNLOADS}/${SRA_ID}.sra ] ; then
    rm -f ${DOWNLOADS}/${SRA_ID}.sra
    echo "${SRA_ID} already downloaded and converted to fastq"
    exit 1
else
	prefetch -X 100G -t ascp -a "${ASPERA}/bin/ascp|${ASPERA}/etc/asperaweb_id_dsa.openssh" ${SRA_ID}
	down_exit=$?
	
    if [ -s ${DOWNLOADS}/${SRA_ID}.sra ] && [ ${down_exit} -eq 0 ] ; then
        cd ${DOWNLOADS}
        fastq-dump --gzip --split-3 ${SRA_ID}.sra
        if [ $? -eq 0 ] ; then
            rm -f ${SRA_ID}.sra
        fi
    else
        echo "${SRA_ID} download went wrong"
        exit 1
    fi
fi







