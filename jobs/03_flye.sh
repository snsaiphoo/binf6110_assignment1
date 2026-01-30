#!/bin/bash 
#SBATCH --job-name=flye
#SBATCH --time=04:00:00
#SBATCH --mem=16GB
#SBATCH --output=/scratch/%u/salmonella/logs/03_%j.out
#SBATCH --error=/scratch/%u/salmonella/logs/03_%j.err

cd /scratch/$USER/salmonella || exit 1

module load apptainer

apptainer exec containers/flye_2.9.6--py311h2de2dd3_0.sif flye \
-o assembly \
--nano-hq filtlong/filteredSRR32410565.fastq.gz \
--genome-size 5m \
--asm-coverage 145 \
-t 16
