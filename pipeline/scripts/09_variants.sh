#!/bin/bash

module load apptainer

cd ..
mkdir -p clair3_out
cd clair3_out

# Run Clair3 variant calling
# input:
# reads.sorted.bam - aligned ONT reads from minimap2
# reference genome - Salmonella reference fasta
# CPU threads count 32
# ont = Ocford Nanopore sequencing
# default model
# include all contigs in variant calling
# output to current directory

apptainer exec ../containers/clair.sif run_clair3.sh \
  --bam_fn=../minimap/reads.sorted.bam \
  --ref_fn=../data/GCF_000006945.2_ASM694v2_genomic.fna \
  --threads=32 \
  --platform=ont \
  --model_path=/opt/models/r1041_e82_400bps_sup_v500 \
  --include_all_ctgs \
  --output=.

