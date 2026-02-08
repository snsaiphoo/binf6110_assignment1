#!/bin/bash
# This script compares the polished assembly to the reference
cd /scratch/$USER/salmonella

module load apptainer

apptainer exec containers/quast_5.3.0--py313pl5321h5ca1c30_2.sif quast.py \
polishing/consensus.fasta \
-R raw/GCF_000006945.2_ASM694v2_genomic.fna \
-o quast/quast_02
