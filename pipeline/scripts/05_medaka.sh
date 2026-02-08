#!/bin/bash
# Run medaka using the filtered data and the assembly output from Flye

cd ..
mkdir -p polishing
cd polishing

module load apptainer

# Running medaka to further polish our assembly
# Using the filtered Filtlong data and the Flye assembly as input
# Using 8 threads

apptainer exec ../containers/medaka.sif medaka_consensus \
-i ../data/filteredSRR32410565.fastq \
-d ../assembly/assembly.fasta \
-o . \
-t 8
