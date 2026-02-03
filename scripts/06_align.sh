#!/bin/bash

module load StdEnv/2023
module load minimap2/2.28

# Docker image can be found in the buildcontainers.sh file

minimap2 -ax asm5 ../raw/GCF_000006945.2_ASM694v2_genomic.fna \
../polishing/consensus.fasta > aln.sam

module load samtools/1.22.1

samtools view -bS aln.sam -o aln.bam
samtools sort aln.bam -o aln.sorted.bam
samtools index aln.sorted.bam

# These files will be used in IGV for visualizations
