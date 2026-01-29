#!/bin/bash

# This code is meant to be ran in your home directory to execute specific jobs in your $SCRATCH environment
# This is with the assumption you are using a DRAC cluster

module load sra-toolkit
prefetch SRR32410565
cd SRR32410565
fasterq-dump SRR32410565.sra
gzip SRR32410565.fastq

# You can use zcat to further view the file to validate the data is correct
