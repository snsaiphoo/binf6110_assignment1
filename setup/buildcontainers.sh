#!/bin/bash

# This script will contain the commands needed to obtain the .sif files needed to run the pipeline
# These commands should be ran in the $SCRATCH environment, specifically in the containers directory 
# The jobs will pull from these directories to run the commands 

# NanoPlot version 1.46.2
# scratch/$USER/salmonella/containers

module load apptainer
singularity pull docker://quay.io/biocontainers/nanoplot:1.46.2--pyhdfd78af_0
