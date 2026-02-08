#!/bin/bash
# This script compares the base assembly to the reference

cd ..
mkdir -p quast
cd quast

module load apptainer

apptainer exec ../containers/quast.sif quast.py ../assembly/assembly.fasta -R
../data/GCF_000006945.2_ASM694v2_genomic.fna -o quast_01
