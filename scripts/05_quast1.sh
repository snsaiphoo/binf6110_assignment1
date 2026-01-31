#!/bin/bash

cd /scratch/$USER/salmonella

mkdir quast

module load apptainer

apptainer exec containers/quast_5.3.0--py313pl5321h5ca1c30_2.sif quast.py polishing/consensus.fasta -R
assembly/assembly.fasta -o quast/quast_01
