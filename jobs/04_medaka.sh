#!/bin/bash 
#SBATCH --job-name=medakapolish
#SBATCH --time=06:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8
#SBATCH --output=/scratch/%u/salmonella/logs/04_%j.out
#SBATCH --error=/scratch/%u/salmonella/logs/04_%j.err

set -e 
set -o pipefail

cd /scratch/$USER/salmonella || exit 1

module load apptainer

apptainer exec /scratch/$USER/salmonella/containers/medaka_2.1.1--py311h1d3aea1_0.sif medaka_consensus \
-i filtlong/filteredSRR32410565.fastq \
-d assembly/assembly.fasta \
-o polishing \
-t 8 \
-m r1041_e82_400bps_sup_v5.2.0
