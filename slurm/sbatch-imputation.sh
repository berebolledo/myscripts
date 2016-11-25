#!/bin/bash

#SBATCH -C new
##SBATCH --ntasks=1
##SBATCH --cpus-per-task=1
#SBATCH -J impute
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=berebolledo@gmail.com
#SBATCH -e slurm-%j.err

step=$1
chrom=$2

data="/nfs/brubeck.bx.psu.edu/scratch4/boris/playground/UDD/gabriela/phase-impute/ref-data/1000GP_Phase3"
# Example:
# 1	249250621
# 2	243199373

REF="${data}/haplotypes/1000GP_Phase3_chr${chrom}.hap.gz"
MAP="${data}/genetic-map/genetic_map_chr${chrom}_combined_b37.txt"
LEGEND="${data}/legend/1000GP_Phase3_chr${chrom}.legend.gz"
PHASED="chr${chrom}.def.phased.haps"



START=$(echo ${step}e6)
END=$(echo $(expr $step + 5)e6)

impute2                    \
    -use_prephased_g       \
    -known_haps_g $PHASED  \
    -h $REF                \
    -l $LEGEND             \
    -m $MAP                \
    -int $START $END       \
    -Ne 20000              \
    -o chunks/chr${chrom}/chr${chrom}.def.phased.imputed_chunk_${START}_${END}