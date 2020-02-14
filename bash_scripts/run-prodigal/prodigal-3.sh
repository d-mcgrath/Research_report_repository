#!/bin/bash
#SBATCH --partition=scavenger # Queue selection
#SBATCH --job-name=prod_3 # Job name
#SBATCH --mail-type=ALL# Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dgellermcgrath@gmail.com # Where to send mail
#SBATCH --ntasks=1 # Run on one CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180gb # Job memory request
#SBATCH --time=2-00:00:00 # Time limit hrs:min:sec
#SBATCH --output=prod_3_%j.log# Standard output/error
#SBATCH --qos=scavenger
#
cd /vortexfs1/scratch/dgellermcgrath/draft-mags-subsets/draft-mags-subset-3
#
module load anaconda/5.1
source activate genex
#
mkdir /vortexfs1/scratch/dgellermcgrath/prodigal_mag_pred_gene_output/anti-results-subset-3
#
for i in *.fa
do
prodigal -i "$i" -o /vortexfs1/scratch/dgellermcgrath/prodigal_mag_pred_gene_output/anti-results-subset-3/"$i"-genes.gbk -p single -d /vortexfs1/scratch/dgellermcgrath/prodigal_mag_pred_gene_output/anti-results-subset-3/"$i"-gene-seqs.fa
done
#
cd /vortexfs1/scratch/dgellermcgrath/prodigal_mag_pred_gene_output/anti-results-subset-3
for i in *.fa
do sed '/^>/s\ #.*\\' "$i" > tmp.fa && mv -f tmp.fa "$i"
done
