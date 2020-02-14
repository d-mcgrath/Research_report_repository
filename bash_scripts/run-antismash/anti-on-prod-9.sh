#!/bin/bash
#SBATCH --partition=compute # Queue selection
#SBATCH --job-name=anti_9 # Job name
#SBATCH --mail-type=ALL # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dgellermcgrath@gmail.com# Where to send mail
#SBATCH --ntasks=36 # Run on one CPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --mem=180gb # Job memory request
#SBATCH --time=2-00:00:00 # Time limit hrs:min:sec
#SBATCH --output=antiProd_9_%j.log # Standard output/error
#SBATCH --qos=unlim
#
#
# run antismash on each individual genome in a directory
#
mkdir /vortexfs1/scratch/dgellermcgrath/prodigal_mag_pred_gene_output/anti-results-subset-9
cd /vortexfs1/scratch/dgellermcgrath/prodigal_mag_pred_gene_output/prodigal-res-subset-9
module load anaconda/5.1
source activate antismash_prod
#
for i in *.fa
do antismash $i --outputfolder /vortexfs1/scratch/dgellermcgrath/prodigal_mag_pred_gene_output/anti-results-subset-9/"$i" -c 36 --knownclusterblast --asf --smcogs --clusterblast
done
