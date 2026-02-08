#!/bin/bash

cd ..
mkdir -p minimap
cd minimap

module load apptainer

REF=../data/GCF_000006945.2_ASM694v2_genomic.fna
READS=../data/SRR32410565.fastq.gz
ASM=../polishing/consensus.fasta

# Indexing reference for variant calling
apptainer exec ../containers/samtools.sif samtools faidx $REF

# Sorted raw read BAM file for variant calling

apptainer exec ../containers/minimap2.sif minimap2 -ax map-ont $REF $READS | \
apptainer exec ../containers/samtools.sif samtools sort -o reads.sorted.bam
apptainer exec ../containers/samtools.sif samtools index reads.sorted.bam

# Sorted assembly BAM file for IGV visualizations

echo "Aligning assembly..."

apptainer exec ../containers/minimap2.sif minimap2 -ax asm5 $REF $ASM | \
apptainer exec ../containers/samtools.sif samtools sort -o aln.sorted.bam
apptainer exec ../containers/samtools.sif samtools index aln.sorted.bam

