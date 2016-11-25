#!/bin/bash
#SBATCH -C new
#SBATCH -J fqdump
#SBATCH --cpus-per-task=1
#SBATCH -t 14-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=berebolledo@gmail.com
#SBATCH -e slurm-%j.err

fastq-dump --split-3 --gzip ${1}
