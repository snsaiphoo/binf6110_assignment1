#!/bin/bash
# This file will run NanoPlot on the raw reads

cd ..
mkdir -p nanoplot
cd nanoplot

module load apptainer

apptainer exec ../containers/nanoplot.sif NanoPlot \
--fastq ../data/SRR32410565.fastq.gz \
-o 01_nanoplot

