#!/bin/bash
#SBATCH --partition=compute # Queue selection
#SBATCH --job-name=prof_11 # Job name
#SBATCH --mail-type=ALL# Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dgellermcgrath@gmail.com # Where to send mail
#SBATCH --ntasks=1 # Run on one CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180gb # Job memory request
#SBATCH --time=2-00:00:00 # Time limit hrs:min:sec
#SBATCH --output=logs/profile_11_%j.log# Standard output/error
#SBATCH --qos=unlim
#
module load python3/3.6.5
module load bio
module load bowtie2/2.3.4.1
module load samtools/1.8
module load anvio/4
cd /vortexfs1/scratch/dgellermcgrath
#
# profile each BAM file
for sample in `awk '{print $1}' sample-partitions/samples_11.txt`
do
 #
    anvi-profile -c REDUNDANT-MAGs-CONTIGS.db \
                 -i "$sample"-in-RMAGs.bam \
                 --skip-SNV-profiling \
                 --num-threads 1 \
                 -o "$sample"-in-RMAGs
#
done
