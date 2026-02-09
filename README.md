# BINF6110 - Genome Assembly: _Salmonella Enterica_

## Introduction
This pipeline aims to assemble a bacterial genome, _Salmonella enterica_, and subsequently compare it to the reference genome. The raw reads are from Oxford Nanopore (ONT) sequencing that have R10 chemistry, and were acquired from NCBI [1]. The general workflow will include data quality checking, genome assembly, assembly analysis, alignment with the reference genome, variant calling, and concluding visualizations to demonstrate differences. 

The _de novo_ genome assembly process starts with raw sequencing reads, either short or long, and assembling them using an algorithm to reconstruct a complete genome [2]. This would produce a consensus sequence that would need to be polished before being used in further analysis or alignment [3]. Currently, genome assembly poses specific issues, such as sequencing bias, errors from sequencing method/platform, repetitive regions, and high computational loads [2,4]. These limitations need to be considered when analyzing assembly results.

Additional challenges come from the alignment process when the complete reference genome is unavailable; in this case, the assembly cannot be fully validated [3]. The main goal of alignment is to see how well the assembly worked in comparison to the reference genome [4]. These challenges can be addressed with the use of the software QUAST [5]. This is a tool to analyze the results of a genome assembly with or without a reference genome, which will be used in the pipeline [5].

Many assemblers attempt to address the mentioned concerns, but here the focus will be on Flye, Canu, and Shasta [3]. Canu is an older assembler that has a very heavy computational load and long run time for error correction and overlap detection [6]. It has been shown to perform well in assembling _Salmonella enterica_, but only with very specific parameters [6]. Shasta is another assembler similar to Flye in terms of performance and is considered more computationally efficient [7]. Shasta, though, has been used a lot for human genomes, whereas Flye has been used for both [7]. Flye is a widely used assembler that has shown to be robust at assembling complex genomic regions and works for both human and bacterial genomes [6-8]. Based on a literature review, Flye was chosen for this pipeline due to its robustness and versatility across different genome types. While this might be more computationally heavy than the Shasta assembler, it will still be faster than Canu [6, 7].

The decisions for the data quality checking, variant calling, and alignment visualizations were chosen through evaluating the literature as well. The software specified in the methods below was shown to be appropriate for handling long read data at this quantity. NanoPlot has been used in a Nanopore quality check pipeline, which is why it has been chosen [9]. Variant calling will be performed with Clair3, as it is to be used for long-read data [10]. Lastly, minimap2, samtools, and IGV have been used in the past for alignment and visualization [3]. 

## Methods 
### Overview
The following section describes the detailed methods used for the _Salmonella enterica_ genome assembly. To execute the pipeline, download the [`pipeline`](./pipeline) directory into your workspace, go to the scripts directory and run each script sequentially, from [`00_ to 11_.sh`](./pipeline/scripts). Output directories will be created automatically at each stage of the process to visualize the results.

### 1.1 - Data Acquisition through SRA 
The _Salmonella enterica_ read data were obtained using the SRA toolkit version 3.2.1, with the specific accession number SRR32410565 [11]. Once the .sra file was downloaded to the terminal, the same toolkit was used to convert it to the fastq format. The specific commands can be found in the [`00_data.sh`](./pipeline/scripts/01_data.sh) file. 

### 2.1 - Nanoplot for Quality Checking 
NanoPlot v1.46.2 was used to evaluate data quality before and after filtering. The results provided length distributions, quartile values, and plots to visualize the data quality [12]. NanoPlot analysis was run on the raw fastq file. The commands for the NanoPlot run can be found in [`02_nanoplot.sh`](./pipeline/scripts/02_nanoplot.sh)

### 2.2 - Filtlong for Data Filtering
Filtlong v0.3.1 was used to filter the raw sequencing reads [13], retaining only reads longer than 1000 bp with a quality score above 90%. These thresholds follow the recommended baseline settings described in the Filtlong documentation. NanoPlot was then rerun to evaluate the quality of the filtered dataset. The exact commands used for this step are provided in the script [`03_filtlong.sh`](./pipeline/scripts/03_filtlong.sh).

### 3.1 - Flye for Genome Assembly
Flye v2.9.6 was selected for this pipeline for the reasons described above, as it has demonstrated strong performance in the literature and is considered both versatile and robust [6–8]. The commands used are provided in the script [`04_flye.sh`](./pipeline/scripts/04_flye.sh). The `--nano-hq` option was specified to accommodate higher-quality ONT R10 reads, while `--asm-coverage 145` was chosen based on the coverage estimates from NanoPlot. The `--genome-size` parameter was set to 5 Mb based on the NanoPlot output and the known genome size. The assembly was executed using 8 CPU threads, with an estimated runtime of approximately 2 hours as reported in the ONT documentation [14, 15].

### 4.1 - Polishing with Medaka
The consensus assembly will be further polished using Medaka v2.1.1. This is commonly used for ONT assemblies from Flye [14, 16, 17]. The commented commands can be found in the [`05_medaka.sh`](./pipeline/scripts/05_medaka.sh) file.

### 5.1 - QUAST for Assembly Analysis
Assembly analysis will be performed using QUAST v5.3.0. This analysis will provide standard assembly metrics such as total length, N50, and misassembly detection [5, 16]. QUAST was run twice to assess the polishing quality against the reference genome from NCBI. The commented commands for both runs can be found in the [`06_quast1.sh`](./pipeline/scripts/06_quast1.sh) and [`07_quast2.sh`](./pipeline/scripts/07_quast2.sh) files.

### 6.1 - Alignment to Reference
The polished assembly from 4.1 was aligned to the reference. The software minimap2 v2.28 will be used, as it is used for long-read alignment [17]. These results produced a SAM file, which was converted to a BAM file using samtools; this served as the input for the visualizations [19]. The specific commands are in the script file [`08_align.sh1`](./pipeline/scripts/08_align.sh). The parameters chosen were to output the results in SAM format (- a) and to use preset high-quality long reads (-x lr:hq) [20].

### 7.1 - Clair3 for Variant Calling 
Clair3 v1.2.0 was used for variant calling; it is a long-read-based variant caller that identifies SNPs and small insertions or deletions relative to the reference genome [10]. samtools v1.22.1 was used to index the reference genome and to convert the raw reads into a sorted BAM file. Genome annotation was performed using Prokka v1.15.6 [21] to generate gene models from the reference assembly. These annotations were then used with bcftools v1.5 to perform functional consequence prediction with bcftools csq [22], enabling classification of variants as synonymous, missense, in-frame indels, or frameshift mutations and allowing inference of predicted amino acid changes. The annotated variant calls were written to the file annotated.vcf and were used for downstream visualization. The commented commands can be found in the following files [`09_variants.sh`](./pipeline/scripts/09_variants.sh) and [`10_varannot.sh`](./pipeline/scripts/10_varannot.sh).

### 8.1 - Alignment Visualizations with IGV
The alignment results were visualized with IGV version 2.19.7 [23]. The BAM file was imported to visualize coverage patterns. Variant calls generated in 7.1 were loaded alongside the alignment to identify sequence differences between the assembly and the reference genome. The functionally annotated reference genome was also imported, allowing variants to be interpreted in their gene context and enabling direct visualization of predicted coding consequences.

### 9.1 - Assembly Visualizations - Circos Plot
A dedicated Conda environment was created to install and manage the dependencies required for genome comparison and visualization. Prokka was run on both the reference and assembled genomes to generate standardized gene annotations in GFF format. These annotations were used as input for SYNY to identify syntenic regions between the genomes [24]. The resulting alignment and feature files were formatted for Circos, which was used to generate a circular visualization of the genome structure. The installation instructions and script to generate the plot can be found in [`11_circos.sh`](./pipeline/scripts/11_circos.sh).

## Results 

### NanoPlot
<p align="center">
  <img src="https://github.com/user-attachments/assets/7488d953-bfb8-4eb4-851a-fe370be5d244" width="80%">
  <br>
  <strong>Figure 1: </strong><em> Read length versus quality score distribution of the original sequencing dataset before filtering, visualized with NanoPlot to assess initial data quality.</em>
</p>
<p align="center">
  <img src="https://github.com/user-attachments/assets/efc672f3-924a-4e6d-840a-d72668987773" width="80%" />
  <br>
  <strong>Figure 2: </strong><em> NanoPlot visualization of sequencing reads after Filtlong filtering, showing the relationship between read length and quality score and highlighting improved read quality.</em>
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

The raw reads had an N50 of about 4,700 bp (Table 1), which is consistent with expectations for R10 chemistry. The median read quality was 23.7, and roughly 77% of reads were above Q20, showing that the dataset was already fairly accurate. However, the quality distribution had a small low-quality tail (Figure 1), with the median higher than the mean, suggesting that a portion of reads was pulling the average down. Because of this, minor filtering was applied using Filtlong with a minimum read length of 1,000 bp and a minimum accuracy of 90%.

After filtering, 157,479 reads totaling 728 Mb were retained. Looking at the results in Table 1, the mean read length (4,625 bp) and N50 (4,903 bp) indicate that read length was mostly preserved while quality improved. The median read quality increased slightly to 23.9, and 80.1% of reads exceeded Q20, indicating enrichment for higher-accuracy reads. The median and mean also became closer together, minimizing the effects of the tail. While a few long reads still showed lower-quality base calls, the overall distribution shifted toward higher-quality base calls. Many long reads were retained, and improvements in accuracy made the filtered dataset more suitable for downstream Flye assembly, where higher read confidence supports better graph construction and consensus.

The filtered dataset retained 728 Mb of sequence. With an assembled genome size of approximately 5.1 Mb, this corresponds to an average sequencing depth of approximately 145×. The output files from both NanoPlot runs are in the [`results/nanoplot`](./results/nanoplot) directory. 

### Flye Assembly
<p align="center">
  <img src="https://github.com/user-attachments/assets/1a58bf15-f3d5-4bcc-b8b3-639c7f7852c9" width="80%" />
  <br></br>
  <strong>Figure 3: </strong><em> Summary of the Flye assembly generated from Filtlong-filtered reads, showing three contigs. Two large contigs represent the main chromosome, while a smaller high-coverage circular contig is consistent with a plasmid.</em>
</p>


The Flye assembly produced three contigs totaling approximately 5.1 Mb. There were two large contigs, with a size of approximately 3.32 Mb and 1.68 Mb, respectively. They each showed similar sequencing coverage (142× and 153×, respectively), suggesting they both come from the same genome. Both contigs were classified as having no repeats, which indicates a stable assembly. A third, smaller contig of 109 kb showed higher coverage (218×) and was predicted to be circular, indicating a plasmid [25]. The absence of repeat flags and alternative graph paths suggests that the assembly is structurally clean and well-supported by reads. Together, these results display a high-confidence bacterial assembly with one possible plasmid and a chromosome represented by two large contigs [25]. The output files from the assembly are in the [`results/assembly`](./results/assembly) directory. 

### Polishing Assembly
<p align="center">
  <img src="https://github.com/user-attachments/assets/32871ef7-25cb-481c-b249-7c08ccab3cca"80%" />
  <br>
  <strong>Figure 4: </strong><em> Cumulative length plot comparing contigs from the consensus assembly and the reference genome.</em>
</p>

<p align="center"><strong>Table 2.</strong> QUAST analysis results comparing the Flye assembly to the polished consensus.</p>

<table align="center">
<tr><th>Metric</th><th>Assembly</th><th>Consensus (Polished)</th></tr>
<tr><td>Genome fraction (%)</td><td>95.669</td><td>95.669</td></tr>
<tr><td>Duplication ratio</td><td>1.002</td><td>1.002</td></tr>
<tr><td>N50</td><td>3,318,776</td><td>3,318,770</td></tr>
<tr><td>Misassemblies</td><td>25</td><td>25</td></tr>
<tr><td>Type</td><td>All relocations</td><td>All relocations</td></tr>
</table>

QUAST analysis comparing the assembly to the _Salmonella enterica serovar Typhimurium LT2_ reference genome (GCF_000006945.2) [26]  showed strong overall agreement [5]. Both the assembly and consensus sequences achieved a genome fraction of approximately 95.7% and a duplication ratio of 1.002, indicating close to complete coverage and low redundancy. The N50 value of 3.3 Mb reflects high contiguity when compared to the 5 Mb reference genome. A total of 25 misassemblies were detected, all classified as relocations. Because the isolate is not identical to the LT2 strain, these events likely represent biological strain variation rather than major assembly errors. The near-identical QUAST metrics between the assembly and polished consensus indicate that polishing preserved global genome structure while maintaining assembly quality. In both cases, there are negligible differences, but the polished consensus was used moving forward. The output files from both QUAST runs are in the [`results/quast`](./results/quast) directory. 

### Variant Calling
<p align="center">
  <img src="https://github.com/user-attachments/assets/02532cfb-e416-418c-b137-8c94dfe38955"80%" />
  <br>
  <strong>Figure 5: </strong><em> Raw variant call set produced by Clair3 prior to downstream filtering and annotation, representing candidate genomic differences from the reference.</em>
</p>

Clair3 was used for variant calling, identifying high-confidence SNPs and small indels relative to the _Salmonella enterica LT2 reference genome_ [10]. Most variants passed filtering with high read depth and genotype quality, and allele frequencies were close to 1, 0.9–1.0. Both substitutions and indels were observed, reflecting normal strain-level differences from the laboratory reference.

Variant summary statistics generated with `bcftools stats` [22] detected 10,502 total variants, mostly SNPs (9,353) and fewer indels (1,166). The transition/transversion ratio was 1.08. Quality score distributions showed many high-confidence calls, including a large number of variants with QUAL ≥30, indicating strong support from the sequencing reads. Indels and substitutions are expected in long-read sequencing because these technologies have known error patterns, especially in low-complexity or homopolymer regions [27]. Overall, the variant patterns are consistent with expected long-read behavior and support the reliability of the Clair3 call set.

<p align="center">
  <img src="https://github.com/user-attachments/assets/30499c71-015c-4eb7-bc24-7b143e410192"80%" />
  <br>
  <strong>Figure 6: </strong><em> Summary of annotated variant calls generated with bcftools, highlighting the number of variants in each functional class. </em>
</p>

Functional annotation using `bcftools csq` with Prokka gene models  [28, 29] showed that most variants were synonymous (n = 2139), followed by missense mutations (n = 865). A smaller number of high-impact variants were present, including frameshifts and stop-gained mutations. Nearly all variants were homozygous; the single heterozygous flag was likely due to noise rather than true biological variation, since bacterial genomes are haploid [30].

The annotated VCF was filtered with bcftools [25] to select representative variants for visualization. Selected variants all had maximum QUAL scores (110). Visualizations with IGV can be seen below. 

### Variant Visualizations
<p align="center">
  <img src="https://github.com/user-attachments/assets/74cabb70-2cc6-424c-b115-74072c06e606" width="80%" />
  <br>
  <strong>Figure 7: </strong><em> Synonymous mutation - IGV visualization of a variant in the polished assembly at NC_003197.2:1101072. The A→T substitution occurs within STM1008, encoding a Gifsy-2 prophage protein, and preserves the amino acid sequence. </em>
</p>
  <br>
<p align="center">
  <img src="https://github.com/user-attachments/assets/bf44e24c-a150-4939-84a7-4d52267a0e77" width="80%" />
  <br>
  <strong>Figure 8: </strong><em> Frameshift mutation - IGV visualization of variant at NC_003197.2:1789147 in the polished assembly, showing a deletion that shifts the coding frame and alters the downstream amino acid sequence. </em>
</p>
  <br>

### Assembly Visualizations
<p align="center">
  <img src="https://github.com/user-attachments/assets/047d140b-84f2-4a2e-95d8-e0fcfec9f6a2" width="80%" />
  <br>
  <strong>Figure 9: </strong><em> Genome-wide Circos comparison between the polished assembly and reference sequence, showing large-scale structural agreement. </em>
</p>
  <br>
  
## Discussion
### Variant Discussion
#### Synonymous Mutation
IGV was used to visualize the variants on the polished assembly. Results for one variant are shown in Figure 7. This image demonstrates a synonymous variant corresponding to an A→T substitution within gene STM1008, which encodes a Gifsy-2 prophage protein [31]. The substitution changes the codon silently, leaving the encoded alanine unchanged, as shown in Figure 7. This confirms a synonymous mutation and shows no change in the primary protein sequence. In the IGV view, the variant shows no evidence of local misalignment and has a Phred-scaled quality score of 110, indicating very high confidence. In this case, the mutation appears neutral.

However, synonymous mutations are not always biologically silent. Recent studies on _Salmonella enterica_ have shown that certain synonymous substitutions can alter mRNA structure, leading to measurable increases in growth rate [32]. This aligns with the broader understanding that Gifsy prophage-associated genes contribute to bacterial fitness and survival [33]. Therefore, even when a specific synonymous variant appears neutral at the protein level, variation within prophage regions could be biologically linked to growth and adaptation.

#### Frameshift Mutation
The frameshift variant observed in Figure 8 occurs within the gene STM1697. This gene is annotated as a putative diguanylate cyclase/phosphodiesterase domain protein. It can be seen that the deletion creates a gap in the consensus relative to the reference sequence. Functional annotation predicts that the frameshift begins at amino acid 179 due to a 2-nucleotide deletion, after which the downstream protein sequence is altered. This shift changes the reading frame and could substantially modify the encoded protein.

STM1697 has been experimentally shown to regulate motility, biofilm formation, and virulence in _Salmonella_ [34]. Because the frameshift disrupts the normal protein structure, this mutation could hypothetically alter how _Salmonella_ regulates environmental adaptation. Since STM1697 operates within a regulatory network that influences bacterial behavior [34], disruption of this gene may affect how the organism responds to environmental conditions. However, these effects remain speculative, and experimental validation would be required to confirm any phenotypic consequences.

Together, these variants illustrate how both neutral and high-impact mutations coexist within the genome, with regulatory disruptions having greater potential to influence bacterial behavior than silent substitutions.

### Assembly Discussion 
#### Circos Plot 
Figure 9 is the assembly analysis that compares the assembled _Salmonella enterica_ genome to the reference. Initial impressions show strong overall structural agreement between the two. In the plot, the reference genome is displayed as two contigs labeled with numbers, while the consensus assembly is labeled using Roman numerals. This allows a direct comparison between both sequences. The two halves of the circos plot are connected by large continuous ribbons, indicating that long genomic regions align in the same order and that the assembly preserved chromosome structure well. 

These localized differences reflect the conserved core genome of _Salmonella enterica_, with most variation restricted to genomic islands associated with mobile elements [35]. The outer composition tracks follow very similar patterns in both the assembly and the reference, supporting the idea that the overall nucleotide structure was maintained. Together, the figure demonstrates that the assembly is globally accurate while still capturing localized variation consistent with natural strain diversity.

## Conlcusion 
This assignment resulted in a reproducible pipeline for genome assembly and variant analysis of the bacterium _Salmonella enterica_. The Circos plot showed that the assembly largely preserved the genome's structure relative to the reference, with only localized relocations, suggesting accurate reconstruction at the whole-genome scale. The variant calling revealed many mutations ranging from synonymous changes to frameshifts. These variants were inspected in IGV to reinforce the importance of visually confirming computational results. 

Looking at the real-world impacts, this assignment highlighted the biological importance of assembling the bacterium _Salmonella enterica_. This process allows researchers to visualize conserved gene regions and identify where strains evolve and differ over time. _Salmonella enterica_ is a foodborne pathogen, and accurate genome reconstruction is beneficial for public health initiatives such as outbreak tracking and long-term genomic surveillance [35]. This was performed on a small scale, and more tools should be investigated to continue improving the assembly process going forward. 

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
[11] A. Klymenko, “08. prefetch and fasterq dump,” GitHub, Sep. 12, 2023. https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump <br/> 
[12] W. D. Coster, “NanoPlot,” GitHub, May 16, 2022. [Online]. Available: https://github.com/wdecoster/NanoPlot <br/> 
[13] R. Wick, “rrwick/Filtlong,” GitHub, Nov. 12, 2022. https://github.com/rrwick/Filtlong <br/> 
[14]  “wf-bacterial-genomes,” Oxford Nanopore Technologies, 2025. [Online]. Available: https://nanoporetech.com/document/epi2me-workflows/wf-bacterial-genomes (accessed Jan. 17, 2026). <br/> 
[15] Pasteur.fr, 2026. [Online]. Available: https://gensoft.pasteur.fr/docs/Flye/2.9/USAGE.html (accessed Jan. 17, 2026). <br/> 
[16] “nanoporetech/medaka,” GitHub, Apr. 15, 2021. [Online]. Available: https://github.com/nanoporetech/medaka <br/> 
[17] J. Y. Lee et al., “Comparative evaluation of Nanopore polishing tools for microbial genome assembly and polishing strategies for downstream analysis,” Scientific Reports, vol. 11, no. 1, Oct. 2021, doi: https://doi.org/10.1038/s41598-021-00178-w. <br/> 
[18] “QUAST 2.3 manual,” Umanitoba.ca, 2026. [Online]. Available: https://home.cc.umanitoba.ca/~psgendb/birchhomedir/doc/local/pkg/rnaQUAST/quast23/manual.html#sec2 (accessed Jan. 17, 2026).
[19] “samtools(1) manual page,” www.htslib.org. [Online]. Available: https://www.htslib.org/doc/samtools.html <br/> 
[20] H. Li, “lh3/minimap2,” GitHub, Feb. 20, 2022. [Online]. Available: https://github.com/lh3/minimap2 <br/> 
[21] tseemann, “tseemann/prokka,” GitHub, Nov. 11, 2019. [Online]. Available: https://github.com/tseemann/prokka <br/>
[22] “bcftools(1),” Github.io, 2025. [Online]. Available: https://samtools.github.io/bcftools/bcftools-man.html#csq (accessed Feb. 08, 2026).<br/>
[23] “Integrated Genomics Viewer - SciLifeLab Courses,” Github.io, 2026. [Online]. Available: https://scilifelab.github.io/courses/rnaseq/labs/IGV (accessed Jan. 17, 2026). <br/> 
[24] PombertLab, “GitHub - PombertLab/SYNY: The SYNY pipeline investigates synteny between species by reconstructing protein clusters from gene pairs.,” GitHub, 2025. [Online]. Available: https://github.com/PombertLab/SYNY?tab=readme-ov-file#Circos-plots (accessed Feb. 05, 2026). <br/> 
[25] J. Johnson, M. Soehnlen, and H. M. Blankenship, “Long read genome assemblers struggle with small plasmids,” Microbial genomics, vol. 9, no. 5, May 2023, doi: https://doi.org/10.1099/mgen.0.001024. <br/>
[26] “Index of /genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2,” Nih.gov, 2025. https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/ (accessed Feb. 08, 2026). <br/>
[27] S. L. Amarasinghe, S. Su, X. Dong, L. Zappia, M. E. Ritchie, and Q. Gouil, “Opportunities and challenges in long-read sequencing data analysis,” Genome Biology, vol. 21, no. 1, Feb. 2020, doi: https://doi.org/10.1186/s13059-020-1935-5. <br/>
[28] T. Seemann, “Prokka: rapid prokaryotic genome annotation,” Bioinformatics, vol. 30, no. 14, pp. 2068–2069, Mar. 2014, doi: https://doi.org/10.1093/bioinformatics/btu153. <br/>
[29] I. Bassano et al., “Evaluation of variant calling algorithms for wastewater-based epidemiology using mixed populations of SARS-CoV-2 variants in synthetic and wastewater samples,” Microbial Genomics, vol. 9, no. 4, Apr. 2023, doi: https://doi.org/10.1099/mgen.0.000933. <br/>
[30] R. K. Holmes and M. G. Jobling, “Genetics,” Nih.gov, 2014. https://www.ncbi.nlm.nih.gov/books/NBK7908/ <br/>
[31] J. S. Bosma, V. Noireaux, and S. D. Bowden, “A Gifsy prophage-encoded protein confers broad phage resistance in Salmonella enterica and is widely distributed across Enterobacteriaceae,” Applied and Environmental Microbiology, vol. 91, no. 12, Dec. 2025, doi: https://doi.org/10.1128/aem.01384-25.<br/>
[32] J. Kristofich et al., “Synonymous mutations make dramatic contributions to fitness when growth is limited by a weak-link enzyme,” PLoS Genetics, vol. 14, no. 8, p. e1007615, Aug. 2018, doi: https://doi.org/10.1371/journal.pgen.1007615. <br/>
[33] S. F. Bailey, L. A. Alonso Morales, and R. Kassen, “Effects of Synonymous Mutations beyond Codon Bias: The Evidence for Adaptive Synonymous Substitutions from Microbial Evolution Experiments,” Genome Biology and Evolution, vol. 13, no. 9, Jun. 2021, doi: https://doi.org/10.1093/gbe/evab141. <br/>
[34] I. Ahmad et al., “The EAL-like protein STM1697 regulates virulence phenotypes, motility and biofilm formation in Salmonella typhimurium,” Molecular Microbiology, vol. 90, no. 6, pp. 1216–1232, Nov. 2013, doi: https://doi.org/10.1111/mmi.12428. <br/>
[35] A. Jacobsen, R. S. Hendriksen, F. M. Aaresturp, D. W. Ussery, and C. Friis, “The Salmonella enterica Pan-genome,” Microbial Ecology, vol. 62, no. 3, pp. 487–504, Jun. 2011, doi: https://doi.org/10.1007/s00248-011-9880-1.
