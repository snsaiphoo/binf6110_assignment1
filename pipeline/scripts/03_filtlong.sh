#!/bin/bash
# These are the commands to filter the raw data using Filtlong

cd ../data

# Location: pipeline/data

module load apptainer

apptainer exec ../containers/filtlong.sif filtlong \
SRR32410565.fastq.gz \
--min_length 1000 \
--keep_percent 90 \
--target_bases 800000000 | gzip > filteredSRR32410565.fastq.gz

# Once the data has been filtered, NanoPlot will be ran again

cd ../nanoplot

apptainer exec ../containers/nanoplot.sif NanoPlot \
--fastq ../data/filteredSRR32410565.fastq.gz \
-o 02_nanoplot
