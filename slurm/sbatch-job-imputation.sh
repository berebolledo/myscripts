#!/bin/bash

#SBATCH -C new
##SBATCH --ntasks=1
##SBATCH --cpus-per-task=1
#SBATCH -J job-imp
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=berebolledo@gmail.com
#SBATCH -e slurm-%j.err


chrom=$1
mkdir -p chunks
mkdir -p chunks/chr${chrom}

chrlen="/nfs/brubeck.bx.psu.edu/scratch4/boris/playground/UDD/gabriela/data/chromosome_length_hg19.txt"
repo="/nfs/brubeck.bx.psu.edu/scratch4/boris/repositories/custom-scripts/slurm/"

round_up () {

	awk -v chr=$1 'BEGIN{OFMT="%.0f"} $1==chr {print $2/1000000}' $chrlen 
}

size=$(round_up $chrom)


for i in $(seq 0 5 $size)
do
    sbatch -N 3-3 -w nn[2-4] -x nn[0-1,5-11] ${repo}/sbatch-imputation.sh $i $chrom
done
