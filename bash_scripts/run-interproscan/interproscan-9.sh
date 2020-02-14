#!/bin/bash
#SBATCH --partition=scavenger # Queue selection
#SBATCH --job-name=interpro_9 # Job name
#SBATCH --mail-type=ALL# Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dgellermcgrath@gmail.com # Where to send mail
#SBATCH --ntasks=1 # Run on one CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180gb # Job memory request
#SBATCH --time=2-00:00:00 # Time limit hrs:min:sec
#SBATCH --output=interpro_9_%j.log# Standard output/error
#SBATCH --qos=scavenger
#
cd /vortexfs1/scratch/dgellermcgrath/scripts/run-interproscan
#
module load anaconda/5.1
source activate genex
#
cat gbk-AA-list-part_09 | while read i
do
./../../interproscan/interproscan-5.40-77.0/interproscan.sh -cpu 36 -f TSV -i ../../january-final-gbks-as-AA-fa/"$i" -d ../../gene-anno-antismash/
done
