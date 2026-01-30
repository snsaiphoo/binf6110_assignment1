#!/bin/bash
#SBATCH --job-name=nanoplotqc
#SBATCH --time=00:10:00
#SBATCH --mem=8GB
#SBATCH --output=/scratch/%u/salmonella/logs/01_%j.out
#SBATCH --error=/scratch/%u/salmonella/logs/01_%j.err

cd /scratch/$USER/salmonella || exit 1

module load apptainer

apptainer exec containers/nanoplot_1.46.2--pyhdfd78af_0.sif NanoPlot \
--fastq raw/SRR32410565/SRR32410565.fastq.gz \
-o nanoplot/01_nanoplot 



