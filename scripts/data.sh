#!/bin/bash

# This code is meant to be ran in the raw directory of your $SCRATCH environment
# This is with the assumption you are using a DRAC cluster

module load sra-toolkit/3.0.9
prefetch SRR32410565
cd SRR32410565
fasterq-dump SRR32410565.sra
gzip SRR32410565.fastq

# You can use zcat to further view the file to validate that the data is correct
