#!/bin/bash
#SBATCH --partition=compute # Queue selection
#SBATCH --job-name=split_n # Job name
#SBATCH --mail-type=ALL# Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dgellermcgrath@gmail.com # Where to send mail
#SBATCH --ntasks=1 # Run on one CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180gb # Job memory request
#SBATCH --time=2-00:00:00 # Time limit hrs:min:sec
#SBATCH --output=logs/split_name_%_j.log# Standard output/error
#SBATCH --qos=unlim
#
module load python3/3.6.5
module load bio
module load bowtie2/2.3.4.1
module load samtools/1.8
module load anvio/4
cd /vortexfs1/scratch/dgellermcgrath
#
for split_name in `sqlite3 REDUNDANT-MAGs-CONTIGS.db 'select split from splits_basic_info'`
do
    # in this loop $split_name goes through names like this: ANW_MAG_00001_000000000001,
    # ANW_MAG_00001_000000000002, ANW_MAG_00001_000000000003, ...; so we can extract
    # the MAG name it belongs to:
    MAG=`echo $split_name | awk 'BEGIN{FS="_"}{print $1"_"$2"_"$3}'`

    # print it out with a TAB character
    echo -e "$split_name\t$MAG"
done > REDUNDANT-MAGs-COLLECTION.txt
