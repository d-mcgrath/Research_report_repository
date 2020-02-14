#!/bin/bash
#SBATCH --partition=scavenger # Queue selection
#SBATCH --job-name=get-rc-15 # Job name
#SBATCH --mail-type=ALL# Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dgellermcgrath@gmail.com # Where to send mail
#SBATCH --ntasks=1 # Run on one CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180gb # Job memory request
#SBATCH --time=1-00:00:00 # Time limit hrs:min:sec
#SBATCH --output=get_read_counts_15_%j.log# Standard output/error
#SBATCH --qos=scavenger
#
cd /vortexfs1/scratch/dgellermcgrath/scripts/bwa-anti-on-prod-final
module load anaconda/5.1
source activate genex
#
cat sam_ao | while read line
do
name=$(echo "$line" | awk '{print $1}' | sed 's\.sam\\')
picard CleanSam I="$line" O="${name}"-clean.sam MAX_RECORDS_IN_RAM=65000000
picard SortSam I="${name}"-clean.sam O="${name}"-sorted.sam SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=65000000
samtools view -F 260 "${name}"-sorted.sam | cut -f 3 | sort | uniq -c | awk '{printf("%s\t%s\n", $2, $1)}' > "${name}"-counts.txt
rm "${name}"-clean.sam
rm "${name}"-sorted.sam
done
