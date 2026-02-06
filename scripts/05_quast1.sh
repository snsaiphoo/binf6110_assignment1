#!/bin/bash
# This script compares the base assembly to the reference
cd /scratch/$USER/salmonella

mkdir quast

module load apptainer

apptainer exec containers/quast_5.3.0--py313pl5321h5ca1c30_2.sif quast.py assembly/assembly.fasta -R
raw/GCF_000006945.2_ASM694v2_genomic.fna -o quast/quast_01
