#!/bin/bash
# These are the commands to filter the raw data using Filtlong

cd /scratch/$USER/salmonella

mkdir filtlong
cd filtlong

module load StdEnv/2023
module load filtlong/0.2.1

filtlong ../raw/SRR32410565/SRR32410565.fastq.gz \
--min_length 1000 \
--keep_percent 90 \
--target_bases 800000000 | gzip > filteredSRR32410565.fastq.gz 

