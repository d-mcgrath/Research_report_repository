#! /usr/bin/python3

import pandas as pd
import glob
import sys

Usage = "USAGE: python3 hw3_counts-to-table.py 'FILE_PATTERN_TO_MATCH' ALIGNER"

# Permitted Aligner names:
# salmon
# kallisto
# bowtie2
# BWA

# Function to return aligner-specific info for reading in files in pandas (whether or not there's a header, columns to select)

def aligner_dets(aligner):
	if (aligner == "bowtie2") or (aligner == "BWA"): # No headers, just contig names & counts
		return(None, [0,1])
	elif aligner == "salmon": # file columns: ["Name", "Length", "EffectiveLength", "TPM", "NumReads"]
		return(0, [0,4])
	elif aligner == "kallisto":	# file columns: ["target_id", length", "eff_length", "est_counts", "tpm"]
		return(0, [0,3])
	else:	# End things now if an incorrect aligner name has been entered
		print("Aligner name must be one of the following: salmon, kallisto, bowtie2, BWA")
		print(Usage)
		exit()

all_files = glob.glob(sys.argv[1])
aligner = sys.argv[2]
head, cols = aligner_dets(aligner)

OutFileName = "read_count_table_" + aligner + ".tsv"

first = 1	# Switch to 'seed' the dataframe with the contig names from the first file, and rename count header to filename (= sample name)
master = []

for filename in all_files:
	if (first == 1):	# Read in first file, and make it into a pandas dataframe
		master = pd.read_csv(filename, header = head, sep='\t', usecols=cols)
		master.columns = names=["Contig", filename]
		first = 0
	else:	# Add count column to master file, and rename count header to filename (= sample name)
		df = pd.read_csv(filename, header = head, sep='\t', usecols=cols)
		df.columns = names=["Contig", filename]
		master = pd.merge(master, df, how = "outer", on="Contig")

master = master.fillna(0) # Turn any blank cells into zero counts rather than NA

master.to_csv(OutFileName, sep = "\t", float_format='%.f', index = False)    # write the table as a tab-delimited file
