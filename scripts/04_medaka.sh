#!/bin/bash
# These commands should be used to move into the $SCRATCH environment, make the new folder for the polishing output files
# Run medaka using the filtered data and the assembly output from Flye
cd /scratch/$USER/salmonella || exit 1

mkdir polishing
cd polishing

module load apptainer

apptainer exec ../containers/medaka_2.1.1--py311h1d3aea1_0.sif medaka_consensus \
-i ../filtlong/filteredSRR32410565.fastq \
-d ../assembly/assembly.fasta \
-o . \
-t 8
