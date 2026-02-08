#!/bin/bash
# This script generates a Circos plot comparing assembly vs reference genome
# Software installation instructions can be found below:
# https://github.com/PombertLab/SYNY?tab=readme-ov-file#Circos-plots

cd ..
mkdir -p prokka_asm
cd prokka_asm

# Annotate polished assembly with Prokka

apptainer exec ../containers/prokka.sif prokka \
  ../polishing/consensus.fasta \
  --outdir .

cd ..
mkdir -p prokka_ref
cd prokka_ref

# Annotate reference genome with Prokka

apptainer exec ../containers/prokka.sif prokka \
  ../data/GCF_000006945.2_ASM694v2_genomic.fna \
  --outdir .

# Run SYNY to generate Circos plot
# Visualization to compare annotated genomes

cd ..
mkdir -p SYNY_run
cd SYNY_run

perl ~/SYNY/run_syny.pl \
  -a ../prokka_asm/asm.gbk ../prokka_ref/ref.gbk \
  -r ref \
  --custom_file ../data/custom_color_2.txt \ #providing a color theme
  --circos all \
  --outdir SYNY_custom


