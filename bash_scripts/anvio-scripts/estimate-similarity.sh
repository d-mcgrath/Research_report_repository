#!/bin/bash
#SBATCH --partition=compute # Queue selection
#SBATCH --job-name=est_sim # Job name
#SBATCH --mail-type=ALL# Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dgellermcgrath@gmail.com # Where to send mail
#SBATCH --ntasks=1 # Run on one CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180gb # Job memory request
#SBATCH --time=2-00:00:00 # Time limit hrs:min:sec
#SBATCH --output=logs/estimate_similarity_%j.log# Standard output/error
#SBATCH --qos=unlim
#
module load python3/3.6.5
module load bio
module load anvio/4
cd /vortexfs1/scratch/dgellermcgrath
#
# processing nucmer outputs, generating the REDUNDANT-MAGs-ANI.txt file
#
python 02-estimate-similarity.py REDUNDANT-MAGs-LENGTH.txt REDUNDANT-MAGs-ANI.txt
