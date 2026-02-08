#!/bin/bash
# This method is going to run the Flye assembly on the filtered data

module load apptainer

cd ..
mkdir -p assembly
cd assembly

apptainer exec ../containers/flye.sif flye \
-o . \
--nano-hq ../data/filteredSRR32410565.fastq.gz \
--genome-size 5m \
--asm-coverage 145 \
-t 8



