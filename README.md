# BINF6110 - Genome Assembly: _Salmonella Enterica_

## Introduction
This pipeline aims to assemble a bacterial genome, Salmonella enterica, and subsequently compare it to the reference genome. The raw reads are from Oxford Nanopore (ONT) sequencing that have R10 chemistry, and were acquired from NCBI [1]. The general workflow will include data quality checking, genome assembly, assembly analysis, alignment with the reference genome, variant calling, and concluding visualizations to demonstrate differences. 

The de novo genome assembly process starts with raw sequencing reads, either short or long, and assembling them using an algorithm to reconstruct a complete genome [2]. This would produce a consensus sequence that would need to be polished before being used in further analysis or alignment [3]. Currently, genome assembly poses specific issues, such as sequencing bias, errors from sequencing method/platform, repetitive regions, and high computational loads [2,4]. These limitations need to be considered when analyzing assembly results.

Additional challenges come from the alignment process when the complete reference genome is unavailable; in this case, the assembly cannot be fully validated [3]. The main goal of alignment is to see how well the assembly worked in comparison to the reference genome [4]. These challenges can be addressed with the use of the software QUAST [5]. This is a tool to analyze the results of a genome assembly with or without a reference genome, which will be used in the pipeline [5].

Many assemblers attempt to address the mentioned concerns, but here the focus will be on Flye, Canu, and Shasta [3]. Canu is an older assembler that has a very heavy computational load and long run time for error correction and overlap detection [6]. It has been shown to have performed a good assembly for Salmonella Enterica, but it was only with very specific parameters set [6]. Shasta is another assembler that is similar to Flye in terms of performance, and can be seen to be more computationally efficient [7]. Shasta, though, has been used a lot for human genomes, whereas Flye has been used for both [7]. Flye is a widely used assembler that has shown to be robust at assembling complex genomic regions and works for both human and bacterial genomes [6-8]. Through the comparison of literature, Flye was chosen for this pipeline due to its robustness and versatility amongst different genome types. While this might be more computationally heavy than the Shasta assembler, it will still be faster than Canu [6, 7].

The decisions for the data quality checking, variant calling, and alignment visualizations were chosen through evaluating the literature as well. The software specified in the methods below was shown to be appropriate for handling long read data at this quantity. NanoPlot has been used in a Nanopore quality check pipeline, which is why it has been chosen [9]. Variant calling will be performed with Clair3, as it is to be used for long-read data [10]. Lastly, minimap2, samtools, and IGV have been used in the past for alignment and visualization [3]. 

## Methods 
### DRAC Cluster 
The proposed pipeline will be run on a DRAC cluster through a SLURM job. While the individual commands are listed here, everything will be consolidated into one script. 

### Acquiring the data from NCBI
The raw reads are to be downloaded from NCBI SRA. The reads will be in SRA format, and converted to FASTQ format. This will be done through the SRA-toolkit to continue with quality checking. 

### Nanoplot for Quality Checking 
Nanoplot v 1.42.0 will be used to evaluate the quality of the long-read raw data. The results will list out read length distributions, quality scores, etc. [11]. The results will be evaluated to determine if Porechop_ABI will be used to further clean the data [9]. 
```
NanoPlot -t 4 -o OUTDIR --fastq file [11]
```
### Running the Genome Assembly with Flye
Flye v 2.9.6 will be used for this pipeline for reasons specified above. It has shown promising results in the literature, and it is versatile and robust [6-8].  
```
flye --nano-hq FILE --genome-size 5m --asm-coverage 160 –t 16
```
The --nano-hq option will be selected to account for higher-quality ONT R10 reads, and --asm-coverage 160 will be used to limit excessive coverage.
The assembly will be run using 16 CPU cores with 64 GB of memory, and an estimated runtime of two hours, as listed in the ONT documentation [12, 13].

### Polishing Assembly Results with Medaka
The consensus assembly will be further polished using Medaka v 1.2.0. This is commonly used for ONT assemblies from Flye [14, 15, 12].   
```
medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC} --bacteria
```

### Assembly Analysis with QUAST
Assembly analysis will be performed using QUAST v 5.20. This analysis will provide standard assembly metrics such as total length, N50, and misassembly detection [5].
```
./quast.py assembly.fasta -R reference.fasta -o OUTDIR [16]
```
### Alignment to Reference
The polished assembly from the previous step will be aligned to the reference genome taken from NCBI. The software minimap2 v 2.28 will be used as it is used for long-read alignment [15]. These results will produce a SAM file that will be converted into a BAM file using samtools 1.22.1, this is the input for the visualizations [17].
```
./minimap2 -ax lr:hq reference.fa assembly.fa > aln.sam [18]
```
```
samtools view -bS aln.sam
samtools sort -o aln.sorted.bam
samtools index aln.sorted.bam [17]
```
### Variant Calling 
The variant calling step will be performed using Clair3, which is run within an Apptainer container on the DRAC cluster. Clair3 is a long-read-based variant caller used to identify SNPs and small deletions or insertions between the assembly and reference [10]. 

### Alignment Visualizations with IGV
The alignment results will be visualized with IGV version 2.19.7 [19]. Import BAM file and Assembly information and visualize the coverage patterns and structural differences. 

## References
[1] NCBI Sequence Read Archive. “SRA Archive: NCBI,” Nih.gov, 2026. [Online]. Available: https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR32410565 (accessed Jan. 17, 2026).  <br/> 
[2] X. Liao, M. Li, Y. Zou, F.-X. Wu, Yi-Pan, and J. Wang, “Current challenges and solutions of de novo assembly,” Quantitative Biology, vol. 7, no. 2, pp. 90–109, Jun. 2019, doi: https://doi.org/10.1007/s40484-019-0166-9. <br/>
[3] M. G. Khrenova et al., “Nanopore Sequencing for De Novo Bacterial Genome Assembly and Search for Single-Nucleotide Polymorphism,” vol. 23, no. 15, pp. 8569–8569, Aug. 2022, doi: https://doi.org/10.3390/ijms23158569.<br/> 
[4] X. Zhu et al., “misFinder: identify mis-assemblies in an unbiased manner using reference and paired-end reads,” BMC Bioinformatics, vol. 16, no. 1, Nov. 2015, doi: https://doi.org/10.1186/s12859-015-0818-3. <br/> 
[5] A. Gurevich, V. Saveliev, N. Vyahhi, and G. Tesler, “QUAST: quality assessment tool for genome assemblies,” Bioinformatics, vol. 29, no. 8, pp. 1072–1075, Feb. 2013, doi: https://doi.org/10.1093/bioinformatics/btt086. <br/> 
[6] A. Schiavone et al., “Factors Affecting the Quality of Bacterial Genomes Assemblies by Canu after Nanopore Sequencing,” Applied Sciences, vol. 12, no. 6, p. 3110, Mar. 2022, doi: https://doi.org/10.3390/app12063110. <br/> 
[7] K. Shafin et al., “Nanopore sequencing and the Shasta toolkit enable efficient de novo assembly of eleven human genomes,” Nature Biotechnology, May 2020, doi: https://doi.org/10.1038/s41587-020-0503-6.<br/> 
[8] M. Kolmogorov, J. Yuan, Y. Lin, and P. A. Pevzner, “Assembly of long, error-prone reads using repeat graphs,” Nature Biotechnology, vol. 37, no. 5, pp. 540–546, Apr. 2019, doi: https://doi.org/10.1038/s41587-019-0072-8. <br/> 
[9] NIAID Nephele Pipeline Documentation. BCBB/NIAID/NIH, “NIAID Nephele,” Nih.gov, 2026. [Online]. Available: https://nephele.niaid.nih.gov/user-guide/pipeline-descriptions/nanopore-qc (accessed Jan. 17, 2026). <br/> 
[10] R. R. Wick, L. M. Judd, T. P. Stinear, and I. R. Monk, “Are reads required? High-precision variant calling from bacterial genome assemblies,” Access Microbiology, vol. 7, no. 5, May 2025, doi: https://doi.org/10.1099/acmi.0.001025.v3. <br/> 
[11] W. D. Coster, “NanoPlot,” GitHub, May 16, 2022. [Online]. Available: https://github.com/wdecoster/NanoPlot <br/> 
[12]  “wf-bacterial-genomes,” Oxford Nanopore Technologies, 2025. [Online]. Available: https://nanoporetech.com/document/epi2me-workflows/wf-bacterial-genomes (accessed Jan. 17, 2026). <br/> 
[13] Pasteur.fr, 2026. [Online]. Available: https://gensoft.pasteur.fr/docs/Flye/2.9/USAGE.html (accessed Jan. 17, 2026). <br/> 
[14] “nanoporetech/medaka,” GitHub, Apr. 15, 2021. [Online]. Available: https://github.com/nanoporetech/medaka <br/> 
[15] J. Y. Lee et al., “Comparative evaluation of Nanopore polishing tools for microbial genome assembly and polishing strategies for downstream analysis,” Scientific Reports, vol. 11, no. 1, Oct. 2021, doi: https://doi.org/10.1038/s41598-021-00178-w. <br/> 
[16] “QUAST 2.3 manual,” Umanitoba.ca, 2026. [Online]. Available: https://home.cc.umanitoba.ca/~psgendb/birchhomedir/doc/local/pkg/rnaQUAST/quast23/manual.html#sec2 (accessed Jan. 17, 2026). <br/>  
[17] “samtools(1) manual page,” www.htslib.org. [Online]. Available: https://www.htslib.org/doc/samtools.html <br/> 
[18] H. Li, “lh3/minimap2,” GitHub, Feb. 20, 2022. [Online]. Available: https://github.com/lh3/minimap2 <br/> 
[19] “Integrated Genomics Viewer - SciLifeLab Courses,” Github.io, 2026. [Online]. Available: https://scilifelab.github.io/courses/rnaseq/labs/IGV (accessed Jan. 17, 2026). <br/> 



