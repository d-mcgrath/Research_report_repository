#!/bin/bash
#SBATCH --partition=scavenger # Queue selection
#SBATCH --job-name=mem_16 # Job name
#SBATCH --mail-type=ALL# Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dgellermcgrath@gmail.com # Where to send mail
#SBATCH --ntasks=1 # Run on one CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180gb # Job memory request
#SBATCH --time=2-00:00:00 # Time limit hrs:min:sec
#SBATCH --output=bwa_mem_16_%j.log# Standard output/error
#SBATCH --qos=scavenger
#
cd /vortexfs1/scratch/dgellermcgrath/scripts/bwa-anti-on-prod-final
#
module load anaconda/5.1
source activate genex
#
cat mt-samples-ap | while read line
do
r1=$(echo "$line" | awk '{print $2}')
r2=$(echo "$line" | awk '{print $3}')
out=$(echo "$line" | awk '{print $1}')
bwa mem -t 36 all_genes.fa ../../cariaco_metatranscriptomes/"$r1" ../../cariaco_metatranscriptomes/"$r2" > "$out".sam
done
