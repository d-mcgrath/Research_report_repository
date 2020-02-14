#!/bin/bash
#SBATCH --partition=compute # Queue selection
#SBATCH --job-name=tigs_db # Job name
#SBATCH --mail-type=ALL # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dgellermcgrath@gmail.com # Where to send mail
#SBATCH --ntasks=36 # Run on one CPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --mem=180gb # Job memory request
#SBATCH --time=2-00:00:00 # Time limit hrs:min:sec
#SBATCH --output=logs/anvi_contigs_db_final_%j.log # Standard output/error
#SBATCH --qos=unlim
#
#
# generate a CONTIGS database, and run HMMs with anvi'o
#
cd /vortexfs1/scratch/dgellermcgrath/
module load python3/3.6.5
module load bio
module load anvio/4
module load hmmer/3.1b2
module load prodigal/2.6.3
#
anvi-gen-contigs-database -f REDUNDANT-MAGs.fa -o REDUNDANT-MAGs-CONTIGS.db
anvi-run-hmms -c REDUNDANT-MAGs-CONTIGS.db --num-threads 36
