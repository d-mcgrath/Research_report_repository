#!/bin/bash
#SBATCH --partition=compute # Queue selection
#SBATCH --job-name=stats_mker # Job name
#SBATCH --mail-type=ALL# Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dgellermcgrath@gmail.com # Where to send mail
#SBATCH --ntasks=1 # Run on one CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180gb # Job memory request
#SBATCH --time=2-00:00:00 # Time limit hrs:min:sec
#SBATCH --output=logs/mag_stats_input_%j.log# Standard output/error
#SBATCH --qos=unlim
#
cd /vortexfs1/scratch/dgellermcgrath
B="REDUNDANT-MAGs-SUMMARY/bin_by_bin/"
echo -e "MAG\tlength\tcompletion\tredundancy" > REDUNDANT-MAGs-STATS.txt
for MAG in `ls $B | grep CarAnox`
do
    length=`cat $B/$MAG/$MAG-total_length.txt`
    completion=`cat $B/$MAG/$MAG-percent_completion.txt`
    redundancy=`cat $B/$MAG/$MAG-percent_redundancy.txt`
    echo -e "$MAG\t$length\t$completion\t$redundancy"
done >> REDUNDANT-MAGs-STATS.txt
