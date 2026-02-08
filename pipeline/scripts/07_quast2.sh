#!/bin/bash
# This script compares the polished assembly to the reference

cd ../quast

module load apptainer

apptainer exec ../containers/quast.sif quast.py \
../polishing/consensus.fasta \
-R ../data/GCF_000006945.2_ASM694v2_genomic.fna \
-o quast_02
