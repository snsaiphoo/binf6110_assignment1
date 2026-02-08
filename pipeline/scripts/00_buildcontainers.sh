#!/bin/bash

# This script will contain the commands needed to obtain the .sif files needed to run the pipeline
# Creates the containers folder to hold the .sif files
# The main jobs will pull from these directories to run the commands
# Apptainer is required, the DRAC cluster module was used

cd ..
mkdir containers
cd containers

module load apptainer

# Location: pipeline/containers

# SRA Toolkit version 3.2.1
singularity pull docker://quay.io/biocontainers/sra-tools:3.2.1--h4304569_1
mv sra-tools_3.2.1--h4304569_1.sif sra-toolkit.sif

# NanoPlot version 1.46.2
singularity pull docker://quay.io/biocontainers/nanoplot:1.46.2--pyhdfd78af_0
mv nanoplot_1.46.2--pyhdfd78af_0.sif nanoplot.sif

# Filtlong version 0.3.1
singularity pull docker://quay.io/biocontainers/filtlong:0.3.1--h077b44d_0
mv filtlong_0.3.1--h077b44d_0.sif filtlong.sif

# Flye version 2.9.6
singularity pull docker://quay.io/biocontainers/flye:2.9.6--py311h2de2dd3_0
mv flye_2.9.6--py311h2de2dd3_0.sif flye.sif

# Medaka version 2.1.1
singularity pull docker://quay.io/biocontainers/medaka:2.1.1--py311h1d3aea1_0
mv medaka_2.1.1--py311h1d3aea1_0.sif medaka.sif

# QUAST version 5.3.0
singularity pull docker://quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2
mv quast_5.3.0--py313pl5321h5ca1c30_2.sif quast.sif

# minimap2 version 2.28
singularity pull docker://quay.io/biocontainers/minimap2:2.28--h577a1d6_4
mv minimap2_2.28--h577a1d6_4.sif minimap2.sif

# samtools version 1.22.1
singularity pull docker://quay.io/biocontainers/samtools:1.22.1--h96c455f_0
mv samtools_1.22.1--h96c455f_0.sif samtools.sif

# Clair3 version 1.2.0
singularity pull docker://hkubal/clair3
mv clair3_latest.sif clair.sif

# Prokka version 1.15.6
singularity pull docker://quay.io/biocontainers/prokka:1.15.6--pl5321hdfd78af_0
mv prokka_1.15.6--pl5321hdfd78af_0.sif prokka.sif

# bcftools version 1.5
singularity pull docker://biocontainers/bcftools:v1.5_cv3
mv bcftools_v1.5_cv3.sif bcftools.sif
