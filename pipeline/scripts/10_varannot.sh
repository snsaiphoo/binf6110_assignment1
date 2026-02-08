#!/bin/bash

# Genome annotation with Prokka
# Annotates the Salmonella reference genome
# Output a prokka_out directory with .gff + gene features
# This will be use to annotate the Clair3 variant calls

cd ..
mkdir -p prokka_out
cd prokka_out

module load apptainer

apptainer exec ../containers/prokka.sif prokka \
  --outdir . \
  --prefix salmonella \
  ../data/GCF_000006945.2_ASM694v2_genomic.fna

# Indexing reference genome (just in case it did not work previously)

apptainer exec ../containers/samtools.sif samtools faidx \
  ../data/GCF_000006945.2_ASM694v2_genomic.fna

# Variant annotation with bcftools
# Use .gff from Prokka to annotate variants
# Output a vcf file

cd ../clair3_out

apptainer exec ../containers/bcftools.sif bcftools csq \
  -f ../data/GCF_000006945.2_ASM694v2_genomic.fna \
  -g ../prokka_out/salmonella.gff \
  highconf.vcf \
  -Ov -o annotated.vcf # output in vcf format


