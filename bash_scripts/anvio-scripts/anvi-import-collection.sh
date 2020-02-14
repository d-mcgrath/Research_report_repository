#!/bin/bash
#SBATCH --partition=compute # Queue selection
#SBATCH --job-name=anv_imp # Job name
#SBATCH --mail-type=ALL# Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dgellermcgrath@gmail.com # Where to send mail
#SBATCH --ntasks=1 # Run on one CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180gb # Job memory request
#SBATCH --time=2-00:00:00 # Time limit hrs:min:sec
#SBATCH --output=logs/anvio_import_%j.log# Standard output/error
#SBATCH --qos=unlim
#
module load python3/3.6.5
module load bio
module load bowtie2/2.3.4.1
module load samtools/1.8
module load anvio/4
#
cd /vortexfs1/scratch/dgellermcgrath/
#
anvi-import-collection REDUNDANT-MAGs-COLLECTION.txt \
                       -c REDUNDANT-MAGs-CONTIGS.db \
                       -p REDUNDANT-MAGs-MERGED/PROFILE.db \
                       -C REDUNDANT_MAGs
