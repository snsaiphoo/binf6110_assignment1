#!/bin/bash

module load apptainer

cd ..
mkdir -p clair3_out
cd clair3_out

# Run Clair3 variant calling
# I: reads.sorted.bam - aligned ONT reads from minimap2
#    reference genome - Salmonella reference fasta

apptainer exec ../containers/clair.sif run_clair3.sh \
  --bam_fn=../minimap/reads.sorted.bam \
  --ref_fn=../data/GCF_000006945.2_ASM694v2_genomic.fna \
  --threads=32 \          # CPU threads count
  --platform=ont \          # Oxford Nanopore sequencing
  --model_path=/opt/models/r1041_e82_400bps_sup_v500 \  # default model
  --include_all_ctgs \     # include all contigs in variant calling
  --output=.  # output to current directory

