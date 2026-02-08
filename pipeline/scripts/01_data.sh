#!/bin/bash
# After loading all the containers, run this file to output your main read data into the data folder

cd ../data

# Location: pipeline/data

module load apptainer

# Using the singualirty to run the SRA Toolkit to acquire the Salmonella Enterica Data 
singularity exec ../containers/sra-toolkit.sif prefetch SRR32410565

cd SRR32410565

# Converting the SRA data from prefetch into a FASTQ file
singularity exec ../../containers/sra-toolkit.sif fasterq-dump SRR32410565.sra

# Zipping the FASTQ file 
gzip SRR32410565.fastq

# You can use zcat to further view the file to validate that the data is correct
