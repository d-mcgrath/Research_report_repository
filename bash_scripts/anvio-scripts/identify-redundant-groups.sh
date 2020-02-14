#!/bin/bash
#SBATCH --partition=compute # Queue selection
#SBATCH --job-name=id_groups # Job name
#SBATCH --mail-type=ALL# Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dgellermcgrath@gmail.com # Where to send mail
#SBATCH --ntasks=1 # Run on one CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180gb # Job memory request
#SBATCH --time=2-00:00:00 # Time limit hrs:min:sec
#SBATCH --output=logs/identify_redundant_groups_%j.log# Standard output/error
#SBATCH --qos=unlim
#
module load python3/3.6.5
module load bio
module load anvio/4
cd /vortexfs1/scratch/dgellermcgrath
# resolve/identify the redundant groups
python 04-identify-redundant-groups.py -d MAGs-PAIRWISE-DATA.txt \
                                       -m REDUNDANT-MAGs-STATS.txt \
                                       -o REDUNDANT-MAGs-GROUPs.txt
