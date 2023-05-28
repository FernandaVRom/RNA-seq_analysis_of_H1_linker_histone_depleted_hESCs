# RNA-seq_analysis_of_H1_linker_histone_depleted_hESCs
This repository contains the code used to perform a comparative gene expression profiling analysis of RNA-seq data for hESCs in its wild-type and H1 linker histone KO derivatives.




## ABOUT THE PIPELINES

The repository is divided into 2 Rmd file and 4 R scripts. Each pipeline serves a specific purpose, as described below:


1. **RNA-seq_quality_assessment_and_processing.Rmd** (for either SKO or DKO analysis)

   Previously, RNAseq read quality was assessed using FASTQC quality control tool. The reads were mapped to human GRCh38.p13 primary assembly using STAR (version 2.6.1). 

   Within this file, you will find instructions for:
   - Deriving gene counts from the number of uniquely aligned unambiguous reads using Subread:featureCount (version v2.0.2).
   - Scale-normalizing library sizes using the trimmed mean of M values (TMM) method with EdgeR software.
   - Calculating weighted likelihoods for all samples based on the observed mean-variance relationship of each gene and sample using the voomWithQualityWeights function in the R package limma.

   * Genes with fold change greater than 1.5 and false discovery rate (FDR) corrected p-value (Benjamini-Hochberg procedure) < 0.05 were considered to be differentially expressed.



2. **RNA-seq_Power_analysis.r**

   This R script performs power analysis for RNA-seq experiments. It estimates the sample size required to achieve a desired level of statistical power per condition. The script utilizes the R package ssizeRNA, which simulates data based on user-defined parameters to assess the experiment's ability to detect differential gene expression.



3. **RNA-seq_Downstream_analysis_of_DEGs.r**

   This R script focuses on the downstream analysis steps following the identification of differentially expressed genes (DEGs). The analysis includes:
   - Generating heat maps to visualize changes in gene expression across different conditions.
   - Conducting functional analysis through gene ontology tools such as gProfiler and ClusterProfiler.



4. **RNA-seq_Synergy_analysis.r** and **functions.r**

   The main R script RNA-seq_Synergy_analysis.r examines the synergy of different H1 linker histones based on the transcriptional output of single knockout (SKO) and double knockout (DKO) conditions. The script takes gene expression data from SKO and DKO samples as input and performs a comparative analysis to evaluate the combined effect of H1 linker histones on gene expression. In order to successfully run, an R script (functions.r) has been incorporated into the same repository, which contains all the functions used for the synergy analysis.

   The R package limma is used to contrast the log fold-changes (LFCs) computed in DKO against the sum of the SKO LFCs, referred to as the additive LFC. The resulting difference is termed the synergy LFC. The obtained p-values for assessing statistical significance provided by limma are further corrected to estimate the false discovery rate (FDR). 

   The resulting genes are categorized based on their potential synergistic behavior, including various scenarios of upregulation and downregulation. Genes with non-significant synergistic LFC are categorized as "SAME" or "undefined," depending on their effect size relative to the average standard errors computed across all LFCs.
   




## ABOUT THE DATA

This section provides information about the data used in the RNA-seq analysis.

The RNA-seq data utilized in this analysis consists of transcriptomic profiles obtained from wild-type hESCs and the different H1 linker histone knockout derivatives. The samples were prepared following established protocols for RNA extraction, library preparation, and sequencing.

The original raw data files generated and used for our RNA-seq analysis have been uploaded to NCBI's Sequence Read Archive (SRA) database following all mandatory requirements. These data can be accessed through Gene Expression Omnibus (GEO) after publication. Additionally, a folder named "Data" has been added to this repository. Once the work is ready for publication, the "Data" folder will include the processed and normalized gene expression data, such as count matrices and metadata.

Please note that access to the raw data files and the "Data" folder may be subject to certain restrictions and data usage policies. For further details on accessing and using the data, please refer to the documentation provided by NCBI and GEO.

For any inquiries or requests regarding the data, please contact the repository maintainers or the corresponding authors.