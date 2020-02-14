#!/bin/bash
#SBATCH --partition=compute # Queue selection
#SBATCH --job-name=gen_tab # Job name
#SBATCH --mail-type=ALL# Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dgellermcgrath@gmail.com # Where to send mail
#SBATCH --ntasks=1 # Run on one CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180gb # Job memory request
#SBATCH --time=2-00:00:00 # Time limit hrs:min:sec
#SBATCH --output=logs/gen_table_for_redundancy_analysis_%j.log# Standard output/error
#SBATCH --qos=unlim
#
module load python3/3.6.5
module load bio
module load anvio/4
cd /vortexfs1/scratch/dgellermcgrath
# generate a table that describes the pairwise relationships between MAGs using previous output files
#
python 03-gen-table-for-redundancy-analysis.py -m REDUNDANT-MAGs-STATS.txt \
                                               -C REDUNDANT-MAGs-PEARSON.txt \
                                               -i REDUNDANT-MAGs-ANI.txt \
                                               -o MAGs-PAIRWISE-DATA.txt
