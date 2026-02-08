# BINF6110 - Genome Assembly: _Salmonella Enterica_

## Introduction
This pipeline aims to assemble a bacterial genome, _Salmonella enterica_, and subsequently compare it to the reference genome. The raw reads are from Oxford Nanopore (ONT) sequencing that have R10 chemistry, and were acquired from NCBI [1]. The general workflow will include data quality checking, genome assembly, assembly analysis, alignment with the reference genome, variant calling, and concluding visualizations to demonstrate differences. 

The _de novo_ genome assembly process starts with raw sequencing reads, either short or long, and assembling them using an algorithm to reconstruct a complete genome [2]. This would produce a consensus sequence that would need to be polished before being used in further analysis or alignment [3]. Currently, genome assembly poses specific issues, such as sequencing bias, errors from sequencing method/platform, repetitive regions, and high computational loads [2,4]. These limitations need to be considered when analyzing assembly results.

Additional challenges come from the alignment process when the complete reference genome is unavailable; in this case, the assembly cannot be fully validated [3]. The main goal of alignment is to see how well the assembly worked in comparison to the reference genome [4]. These challenges can be addressed with the use of the software QUAST [5]. This is a tool to analyze the results of a genome assembly with or without a reference genome, which will be used in the pipeline [5].

Many assemblers attempt to address the mentioned concerns, but here the focus will be on Flye, Canu, and Shasta [3]. Canu is an older assembler that has a very heavy computational load and long run time for error correction and overlap detection [6]. It has been shown to have performed a good assembly for _Salmonella Enterica_, but it was only with very specific parameters set [6]. Shasta is another assembler that is similar to Flye in terms of performance, and can be seen to be more computationally efficient [7]. Shasta, though, has been used a lot for human genomes, whereas Flye has been used for both [7]. Flye is a widely used assembler that has shown to be robust at assembling complex genomic regions and works for both human and bacterial genomes [6-8]. Through the comparison of literature, Flye was chosen for this pipeline due to its robustness and versatility amongst different genome types. While this might be more computationally heavy than the Shasta assembler, it will still be faster than Canu [6, 7].

The decisions for the data quality checking, variant calling, and alignment visualizations were chosen through evaluating the literature as well. The software specified in the methods below was shown to be appropriate for handling long read data at this quantity. NanoPlot has been used in a Nanopore quality check pipeline, which is why it has been chosen [9]. Variant calling will be performed with Clair3, as it is to be used for long-read data [10]. Lastly, minimap2, samtools, and IGV have been used in the past for alignment and visualization [3]. 

## Methods 
### Overview
The following section describes the detailed methods used for the _Salmonella enterica_ genome assembly. To execute the pipeline, download the [`pipeline`](./pipeline) directory into your workspace, go to the scripts directory and run each script sequentially, from [`00_ to 11_.sh`](./pipeline/scripts). Output directories will be created automatically at each stage of the process to visualize the results.

### 1.1 - Data Acquisition through SRA 
The Salmonella Enterica read data were obtained using the SRA toolkit version 3.2.1, with the specific accession number SRR32410565 [20]. Once the .sra file was downloaded to the terminal, the same toolkit was used to convert it to the fastq format. The specific commands can be found in the [`00_data.sh`](./pipeline/scripts/01_data.sh) file. 

### 2.1 - Nanoplot for Quality Checking 
NanoPlot v1.46.2 was used to evaluate data quality before and after filtering. The results provided length distributions, quartile values, and plots to visualize the data quality [11]. NanoPlot analysis was run on the raw fastq file. The commands for the NanoPlot run can be found in [`02_nanoplot.sh`](./pipeline/scripts/02_nanoplot.sh)

### 2.2 - Filtlong for Data Filtering
Filtlong v0.3.1 was used to filter the raw sequencing reads [21], retaining only reads longer than 1000 bp with a quality score above 90%. These thresholds follow the recommended baseline settings described in the Filtlong documentation. NanoPlot was then rerun to evaluate the quality of the filtered dataset. The exact commands used for this step are provided in the script [`03_filtlong.sh`](./pipeline/scripts/03_filtlong.sh).

### 3.1 - Flye for Genome Assembly
Flye v2.9.6 was selected for this pipeline for the reasons described above, as it has demonstrated strong performance in the literature and is considered both versatile and robust [6–8]. The commands used are provided in the script [`04_flye.sh`](./pipeline/scripts/04_flye.sh). The `--nano-hq` option was specified to accommodate higher-quality ONT R10 reads, while `--asm-coverage 145` was chosen based on the coverage estimates from NanoPlot. The `--genome-size` parameter was set to 5 Mb based on the NanoPlot output and the known genome size. The assembly was executed using 8 CPU threads, with an estimated runtime of approximately 2 hours as reported in the ONT documentation [12, 13].

### 4.1 - Polishing with Medaka
The consensus assembly will be further polished using Medaka v2.1.1. This is commonly used for ONT assemblies from Flye [12, 14, 15]. The commented commands can be found in the [`05_medaka.sh`](./pipeline/scripts/05_medaka.sh) file.

### 5.1 - QUAST for Assembly Analysis
Assembly analysis will be performed using QUAST v5.3.0. This analysis will provide standard assembly metrics such as total length, N50, and misassembly detection [5]. QUAST was run twice to assess the polishing quality against the reference genome from NCBI. The commented commands for both runs can be found in the [`06_quast1.sh`](./pipeline/scripts/06_quast1.sh) and [`07_quast2.sh`](./pipeline/scripts/07_quast2.sh) files.

### 6.1 - Alignment to Reference
The polished assembly from 4.1 was aligned to the reference. The software minimap2 v2.28 will be used, as it is used for long-read alignment [15]. These results produced a SAM file, which was converted to a BAM file using samtools; this served as the input for the visualizations [17]. The specific commands are in the script file [`08_align.sh1`](./pipeline/scripts/08_align.sh). The parameters chosen were to output the results in SAM format (- a) and to use preset high-quality long reads (-x lr:hq).

### 7.1 - Clair3 for Variant Calling 
Clair3 v1.2.0 was used for variant calling; it is a long-read-based variant caller that identifies SNPs and small insertions or deletions relative to the reference genome [10]. samtools v1.22.1 was used to index the reference genome and to convert the raw reads into a sorted BAM file. Genome annotation was performed using Prokka v1.15.6 to generate gene models from the reference assembly. These annotations were then used with bcftools v1.5 to perform functional consequence prediction with bcftools csq, enabling classification of variants as synonymous, missense, in-frame indels, or frameshift mutations and allowing inference of predicted amino acid changes. The annotated variant calls were written to the file annotated.vcf and were used for downstream visualization. The commented commands can be found in the following files [`09_variants.sh`](./pipeline/scripts/09_variants.sh) and [`10_varannot.sh`](./pipeline/scripts/10_varannot.sh).

### 8.1 -Alignment Visualizations with IGV
The alignment results were visualized with IGV version 2.19.7 [19]. The BAM file was imported to visualize coverage patterns. Variant calls generated in 7.1 were loaded alongside the alignment to identify sequence differences between the assembly and the reference genome. The functionally annotated reference genome was also imported, allowing variants to be interpreted in their gene context and enabling direct visualization of predicted coding consequences.

### 9.1 - Assembly Visualizations - Circos Plot
A dedicated Conda environment was created to install and manage the dependencies required for genome comparison and visualization. Prokka was run on both the reference and assembled genomes to generate standardized gene annotations in GFF format. These annotations were used as input for SYNY to identify syntenic regions between the genomes [22]. The resulting alignment and feature files were formatted for Circos, which was used to generate a circular visualization of the genome structure. The installation instructions and script to generate the plot can be found in [`11_circos.sh`](./pipeline/scripts/11_circos.sh).

## Results 

### NanoPlot


<p align="center">
  <img src="https://github.com/user-attachments/assets/7488d953-bfb8-4eb4-851a-fe370be5d244" width="80%">
  <br>
  <strong>Figure 1: </strong><em> NanoPlot read length vs quality distribution before filtering.</em>
</p>
<p align="center">
  <img src="https://github.com/user-attachments/assets/efc672f3-924a-4e6d-840a-d72668987773" width="80%" />
  <br>
  <strong>Figure 2: </strong><em> NanoPlot read length vs quality distribution after filtering.</em>
</p>

<p align="center"><strong>Table 1.</strong> Summary statistics of raw and filtered sequencing reads.</p>

<table align="center">
<tr><th>Metric</th><th>Raw Data</th><th>Filtered Data</th></tr>
<tr><td>Number of reads</td><td>196,031</td><td>157,479</td></tr>
<tr><td>Mean read length (bp)</td><td>4,128.4</td><td>4,625.2</td></tr>
<tr><td>N50 (bp)</td><td>4,683</td><td>4,903</td></tr>
<tr><td>Median read quality</td><td>23.7</td><td>23.9</td></tr>
<tr><td>Mean quality</td><td>18.9</td><td>20.0</td></tr>
<tr><td>Reads > Q20</td><td>150,827 (76.9%) — 625.5 Mb</td><td>126,072 (80.1%) — 577.9 Mb</td></tr>
</table>

The raw reads had an N50 of about 4,700 bp (Table 1), which is consistent with expectations for R10 chemistry. The median read quality was 23.7, and roughly 77% of reads were above Q20, showing that the dataset was already fairly accurate. However, the quality distribution had a small low-quality tail (Figure 1), with the median higher than the mean, suggesting that a portion of reads were pulling the average down. Because of this, light filtering was applied using Filtlong with a minimum read length of 1,000 bp and a minimum accuracy of 90%.

After filtering, 157,479 reads totaling 728 Mb were retained. The mean read length (4,625 bp) and N50 (4,903 bp) show that read length was mostly preserved while quality improved. The median read quality increased slightly to 23.9, and 80.1% of reads exceeded Q20, indicating enrichment for higher-accuracy reads. While a few long reads still showed lower quality, the overall distribution shifted toward more reliable base calls. This balance between keeping long reads and improving accuracy made the filtered dataset more suitable for downstream Flye assembly, where higher read confidence supports better graph construction and consensus.

The filtered dataset retained 728 Mb of sequence. With an assembled genome size of approximately 5.1 Mb, this corresponds to an average sequencing depth of ~145× coverage.

## Flye Assembly



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
[16] “QUAST 2.3 manual,” Umanitoba.ca, 2026. [Online]. Available: https://home.cc.umanitoba.ca/~psgendb/birchhomedir/doc/local/pkg/rnaQUAST/quast23/manual.html#sec2 (accessed Jan. 17, 2026).
[17] “samtools(1) manual page,” www.htslib.org. [Online]. Available: https://www.htslib.org/doc/samtools.html <br/> 
[18] H. Li, “lh3/minimap2,” GitHub, Feb. 20, 2022. [Online]. Available: https://github.com/lh3/minimap2 <br/> 
[19] “Integrated Genomics Viewer - SciLifeLab Courses,” Github.io, 2026. [Online]. Available: https://scilifelab.github.io/courses/rnaseq/labs/IGV (accessed Jan. 17, 2026). <br/> 
[20] A. Klymenko , “08. prefetch and fasterq dump,” GitHub, Sep. 12, 2023. https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump <br/> 
[21] R. Wick, “rrwick/Filtlong,” GitHub, Nov. 12, 2022. https://github.com/rrwick/Filtlong <br/> 
[22] PombertLab, “GitHub - PombertLab/SYNY: The SYNY pipeline investigates synteny between species by reconstructing protein clusters from gene pairs.,” GitHub, 2025. https://github.com/PombertLab/SYNY?tab=readme-ov-file#Circos-plots (accessed Feb. 05, 2026). <br/> 




