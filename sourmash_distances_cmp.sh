#!/bin/sh
#SBATCH --time=06-23:59:59
#SBATCH --partition=bio-compute
#SBATCH --mem=250G
#SBATCH --job-name=sourmash
cd /mnt/scratch2/users/40309916/E_coli_genomes/genomes/Prodigal_results/Genomes
module load apps/anaconda3
source activate /mnt/scratch2/users/40309916/panaroo

sourmash compare *.sig -o distances_updated.cmp -k 31
