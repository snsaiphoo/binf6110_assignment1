# BINF6110 - Genome Assembly: Salmonella Enterica

## Introduction
This pipeline aims to assemble a bacterial genome, Salmonella enterica, and subsequently compare it to the reference genome. The raw reads are from Oxford Nanopore (ONT) sequencing that have R10 chemistry, and were acquired from NCBI [1]. The general workflow will include data quality checking, genome assembly, assembly analysis, alignment with the reference genome, and concluding visualizations to demonstrate differences. 

The de novo genome assembly process starts with taking raw reads, either short or long, and assembling them using an algorithm to reconstruct a complete genome [2]. This would produce a consensus sequence that would need to be polished before being used in further analysis or alignment. Currently, genome assembly poses specific issues, such as sequencing bias, errors from sequencing method/platform, repetitive regions, and high computational load [2,3]. These limitations need to be considered when analyzing assembly results. 


### References
[1] https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR32410565
[2] https://onlinelibrary.wiley.com/doi/10.1007/s40484-019-0166-9
[3] https://link.springer.com/article/10.1186/s12859-015-0818-3



