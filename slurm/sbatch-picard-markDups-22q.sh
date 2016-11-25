#!/bin/bash

#SBATCH -C new
#SBATCH --cpus-per-task=4
#SBATCH -J markDup
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=berebolledo@gmail.com
#SBATCH -e slurm-%j.err

#H53NGBBXX:1:1115:25428:4725

picard MarkDuplicates        \
    I=${1}                   \
    O=dups.${1}              \
    M=${1}.metrics           \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR=tmp
