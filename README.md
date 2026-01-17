# BINF6110 - Genome Assembly: _Salmonella Enterica_

## Introduction
This pipeline aims to assemble a bacterial genome, Salmonella enterica, and subsequently compare it to the reference genome. The raw reads are from Oxford Nanopore (ONT) sequencing that have R10 chemistry, and were acquired from NCBI [1]. The general workflow will include data quality checking, genome assembly, assembly analysis, alignment with the reference genome, and concluding visualizations to demonstrate differences. 

The de novo genome assembly process starts with taking raw reads, either short or long, and assembling them using an algorithm to reconstruct a complete genome [2]. This would produce a consensus sequence that would need to be polished before being used in further analysis or alignment. Currently, genome assembly poses specific issues, such as sequencing bias, errors from sequencing method/platform, repetitive regions, and high computational load [2,3]. These limitations need to be considered when analyzing assembly results. 

Some issues come from the process of alignment. These include not having a complete reference genome; in this case, the assembly can’t be fully validated. The main goal of alignment is to see how well the assembly worked in comparison to the reference genome [3]. These challenges can be addressed with the use of the software QUAST [4]. This is a tool to analyze the results of a genome assembly with or without a reference genome [4]. 

For this pipeline, three assemblers that attempt to address these concerns were researched. These were Flye, Canu, and Shasta [5]. Canu is an older assembler that has a very heavy computational load and long run time [6]. It has seen to have performed a good assembly for Salmonella Enterica, but it was only with very specific parameters set [6]. Shasta is another assembler that is similar to Flye in terms of performance, but can be seen to be cheaper [7]. Flye is a widely used assembler that has shown to be robust and work for both human and bacterial genomes [6, 7, 8]. There is also literature to show that Flye was performing at the level of Shasta when compared to a human genome example [7]. Through the comparison of literature, Flye was chosen for this pipeline due to its robustness and versatility amongst different genome types. While this might be more computationally heavy than the Shasta assembler, it will still be faster than Canu. 

## Methods 
1. Acquiring the data. The raw reads are to be downloaded from trace.ncbi. The download will be in SRA format, so this will need to be converted to FASTQ format.
2. Nanoplot v 1.42.0 was chosen to evaluate the quality of the long-read raw data. [9] Here the data will be investigated, and depending on the output Porechop_ABI will be 3. used to further clean the data. Default settings 
4. Flye:2.9.6 → genome size 5Mbp, 16 CPUs, memory 64GB [12, 13] –nano-hq, based off description of our data, asm-coverage 160 CPU time 2 hours
5. Medaka 1.2.0→ polishing step, this is commonly used of ONT data and in conjunction with Flye [10, 11]   
6. medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC} --bacteria
7. QUAST 5.20 → is then performed to analyze the genome assembly quality [4]
8. Alignment →  
   - minimap2 2.28 [14] ./minimap2 -ax lr:hq ref.fa ont-Q20.fq.gz > aln.sam       # Nanopore Q20 genomic reads (v2.27+) [14]
   - + samtools 1.22.1, this is to generate the BAM file that can be used for alignment in IGV of the assembly and the reference genome [5, 15]
9. IGV → visualizations, version 2.19.7 [5, 16]


### References
[1] “SRA Archive: NCBI,” Nih.gov, 2026. https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR32410565 (accessed Jan. 17, 2026).<br/>
[2] X. Liao, M. Li, Y. Zou, F.-X. Wu, Yi-Pan, and J. Wang, “Current challenges and solutions of de novo assembly,” Quantitative Biology, vol. 7, no. 2, pp. 90–109, Jun. 2019, doi: https://doi.org/10.1007/s40484-019-0166-9.<br/>
[3] X. Zhu et al., “misFinder: identify mis-assemblies in an unbiased manner using reference and paired-end reads,” BMC Bioinformatics, vol. 16, no. 1, Nov. 2015, doi: https://doi.org/10.1186/s12859-015-0818-3.<br/>
[4] A. Gurevich, V. Saveliev, N. Vyahhi, and G. Tesler, “QUAST: quality assessment tool for genome assemblies,” Bioinformatics, vol. 29, no. 8, pp. 1072–1075, Feb. 2013, doi: https://doi.org/10.1093/bioinformatics/btt086.<br/>
[5] M. G. Khrenova et al., “Nanopore Sequencing for De Novo Bacterial Genome Assembly and Search for Single-Nucleotide Polymorphism,” vol. 23, no. 15, pp. 8569–8569, Aug. 2022, doi: https://doi.org/10.3390/ijms23158569.<br/>
[6] A. Schiavone et al., “Factors Affecting the Quality of Bacterial Genomes Assemblies by Canu after Nanopore Sequencing,” Applied Sciences, vol. 12, no. 6, p. 3110, Mar. 2022, doi: https://doi.org/10.3390/app12063110.<br/>
[7] K. Shafin et al., “Nanopore sequencing and the Shasta toolkit enable efficient de novo assembly of eleven human genomes,” Nature Biotechnology, May 2020, doi: https://doi.org/10.1038/s41587-020-0503-6.<br/>
[8] M. Kolmogorov, J. Yuan, Y. Lin, and P. A. Pevzner, “Assembly of long, error-prone reads using repeat graphs,” Nature Biotechnology, vol. 37, no. 5, pp. 540–546, Apr. 2019, doi: https://doi.org/10.1038/s41587-019-0072-8.<br/>
[9] BCBB/NIAID/NIH, “NIAID Nephele,” Nih.gov, 2026. https://nephele.niaid.nih.gov/user-guide/pipeline-descriptions/nanopore-qc (accessed Jan. 17, 2026).<br/>
[10] “nanoporetech/medaka,” GitHub, Apr. 15, 2021. https://github.com/nanoporetech/medaka<br/>
[11] J. Y. Lee et al., “Comparative evaluation of Nanopore polishing tools for microbial genome assembly and polishing strategies for downstream analysis,” Scientific Reports, vol. 11, no. 1, Oct. 2021, doi: https://doi.org/10.1038/s41598-021-00178-w.<br/>
[12] “wf-bacterial-genomes,” Oxford Nanopore Technologies, 2025. https://nanoporetech.com/document/epi2me-workflows/wf-bacterial-genomes (accessed Jan. 17, 2026).<br/>
[13] Pasteur.fr, 2026. https://gensoft.pasteur.fr/docs/Flye/2.9/USAGE.html (accessed Jan. 17, 2026).<br/>
[14] H. Li, “lh3/minimap2,” GitHub, Feb. 20, 2022. https://github.com/lh3/minimap2<br/>
[15] “samtools(1) manual page,” www.htslib.org. https://www.htslib.org/doc/samtools.html<br/>
[16] “Integrated Genomics Viewer - SciLifeLab Courses,” Github.io, 2026. https://scilifelab.github.io/courses/rnaseq/labs/IGV (accessed Jan. 17, 2026).<br/>



