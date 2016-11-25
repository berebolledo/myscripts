#!/bin/bash

#SBATCH -C new
#SBATCH -J sort-picard
#SBATCH --cpus-per-task=4
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=berebolledo@gmail.com
#SBATCH -e slurm-%j.err

picard SortSam    \
    I=${1}        \
    O=sorted.${1} \
    SO=coordinate \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR=tmp

