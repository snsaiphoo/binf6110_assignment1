#!/bin/bash

# This script will contain the commands needed to obtain the .sif files needed to run the pipeline
# The code below moves you to the $SCRATCH environment
# Creates the containers folder to hold the .sif files
# The main jobs will pull from these directories to run the commands 

cd /scratch/$USER/salmonella
mkdir containers
cd containers

# Location: scratch/$USER/salmonella/containers

# NanoPlot version 1.46.2
singularity pull docker://quay.io/biocontainers/nanoplot:1.46.2--pyhdfd78af_0

# Flye version 2.9.6
singularity pull docker://quay.io/biocontainers/flye:2.9.6--py311h2de2dd3_0

# Medaka version 2.1.1
singularity pull docker://quay.io/biocontainers/medaka:2.1.1--py311h1d3aea1_0

# QUAST version 5.3.0
singularity pull docker://quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2
