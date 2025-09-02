# Transcriptional analysis of ASCL1 and the effects of neurogenesis in GABAergic Neurons
 
![Powerpoint  Transcriptional analysis of ASCL1 and the effects of neurogenesis in GABAergic Neurons_page-0001](https://github.com/user-attachments/assets/96f9c6f7-b374-49fa-8edd-6f5091df033a)

This project focuses on the transcriptional analysis of the ASCL1 gene and its role in neurogenesis, particularly in the development of GABAergic neurons. Through the integration of multiple genomic datasets and bioinformatics tools, this project analyzes gene expression patterns, identifies key transcription factors, and examines co-expression groups that influence neurogenesis. The study leverages various data engineering techniques for dataset acquisition, processing, and integration, along with data analysis and visualization methods to draw biological insights.

## Project Overview
**Institution**: University of Portsmouth  
**Program**: BSc (Hons) Biochemistry  
**Year**: 2020/2021  
**Supervisor**: Dr. Frank Schubert

## Table of Contents

1. [Data Engineering](#data-engineering)
    - [Data Acquisition and Storage](#data-acquisition-and-storage)
    - [Data Processing](#data-processing)
    - [Data Integration](#data-integration)
2. [Data Analysis & Visualization](#data-analysis--visualization)
    - [Differential Gene Expression Analysis](#differential-gene-expression-analysis)
    - [Gene Ontology (GO) Term Analysis](#gene-ontology-go-term-analysis)
    - [STRING Analysis](#string-analysis)
3. [Tools and Technologies](#tools-and-technologies)
4. [Results](#results)
5. [Future Research](#future-research)

## Data Engineering

### Data Acquisition and Storage
![Figure 6](https://github.com/user-attachments/assets/ad12dc16-640b-43c4-8eac-3ab846aa1aa0)
- **Datasets**: The study utilised five publicly available datasets from the Gene Expression Omnibus (GEO) and ArrayExpress, including RNA-Seq and microarray data. These datasets were selected based on relevance to ASCL1 and GABAergic neuron development.
- **Repositories**: 
  - Gene Expression Omnibus (GEO): GSE29985, GSE31635, GSE78949, GSE46791
  - ArrayExpress: E-MTAB-4840
- **Data Storage**: Raw data was collected and stored using platforms such as ArrayExpress and GEO, ensuring proper organization and versioning.

### Data Processing

The project required extensive preprocessing of genomic data in Galaxy Server:
- **Quality Control**: FASTQC was used to assess the quality of RNA-Seq data, identifying and removing low-quality sequences with tools like trimmomatic.
- **Alignment**: HISAT2 was used to align the cleaned RNA-Seq reads to the reference genome (GRCh38.p13) to generate aligned read counts.
- **Quantification**: Transcript expression was measured using RPKM (Reads Per Kilobase of transcript per Million mapped reads) and quantified with tools like StringTie.
- **Differential Gene Expression Analysis**: DESeq2 was used to compute differential expression between samples. This step transformed raw counts into meaningful biological insights by identifying upregulated and downregulated genes.

### Data Integration

- **Merged Gene Lists**: The project employed Venn diagrams generated via InteractiVenn to compare overlapping genes between datasets. This enabled identification of key genes involved in GABAergic neuron development.
- **GO Term Classification**: Gene lists from multiple datasets were analyzed using g:Profiler, a tool for statistical GO Term enrichment, providing insights into molecular function, biological processes, and cellular components.


## Data Analysis & Visualization

### Differential Gene Expression Analysis

- **Volcano Plots**: GEO2R and Galaxy were used to generate volcano plots, visualizing statistically significant gene expression changes between samples. Upregulated and downregulated genes were identified based on log2 fold change and p-adjusted values.

### Gene Ontology (GO) Term Analysis

- **GO Term Enrichment**: Overlapping gene lists were analyzed with g:Profiler to classify genes into molecular function, biological process, and cellular component categories. This provided a biological context for the identified genes.
  - Example GO terms: Neurogenesis, neuron differentiation, calcium ion binding, and synaptic communication.

### STRING Analysis
- **Protein-Protein Interaction (PPI) Networks**: STRING was used to visualize networks of interacting proteins. This analysis highlighted interactions among key proteins such as GAD2, CALB1/2, and NOTCH2, which are involved in neurogenesis and synaptic communication in GABAergic neurons.


## Tools and Technologies

- **Data Acquisition**: ArrayExpress, Gene Expression Omnibus (GEO), GEO2R
- **Data Processing**: FASTQC, Trimmomatic, HISAT2, DESeq2, StringTie
- **Data Visualization**: GEO2R, Galaxy, InteractiVenn, g:Profiler, STRING
- **Programming Languages**: Python (for Galaxy), R (for differential expression analysis)
- **Operating Systems**: Web-based bioinformatics tools and Linux for local data processing

---

## Results

- **Key Findings**: The analysis identified critical genes such as GAD2, CALB1, CALB2, and NOTCH2, which play roles in GABAergic neurogenesis. GO term analysis revealed associations with neurogenesis and neuron differentiation, while STRING analysis uncovered protein-protein interactions related to synaptic communication.
- **Data Visualization**: Principal Component Analysis (PCA) and Uniform Manifold Approximation and Projection (UMAP) plots were generated to identify clustering patterns and relationships among different samples.

---

## Future Research

- **Adult Neurogenesis**: Future studies could focus on conducting RNA-Seq analysis on adult human brain samples to compare embryonic and adult neurogenesis.
- **Slit2 Function**: Further research is needed to clarify the role of Slit2 in axonal guidance and its interaction with ASCL1 in neurogenesis.

## Authors

- [@garcane](https://github.com/garcane)


## Acknowledgements
* [Galaxy](https://galaxyproject.org/) is an open, web-based platform for data-intensive biomedical research. It allows users to perform accessible, reproducible, and transparent computational research.
* [GEO2R](https://www.ncbi.nlm.nih.gov/geo/geo2r/) is an interactive web tool that allows users to compare two or more groups of samples in the Gene Expression Omnibus (GEO) database to identify differentially expressed genes.
* [g:Profiler](https://biit.cs.ut.ee/gprofiler/) is a web server for functional enrichment analysis and conversions of gene lists, identifiers, and various other information related to gene function and pathways.
* [STRING](https://string-db.org/) is a database of known and predicted protein-protein interactions. It includes direct (physical) and indirect (functional) associations.
* [InteractiVenn](http://www.interactivenn.net/) is a web-based tool for the visualization of Venn diagrams for multiple sets with a user-friendly interface.
* [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/) is a public database of gene expression data and a repository for high-throughput functional genomics data.
* [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/) is a public repository that archives and freely distributes comprehensive sets of gene expression data.
* [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a quality control tool for high-throughput sequence data. It provides summaries and visualizations to assess sequencing quality.
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) is a flexible read trimming tool for Illumina NGS data. It removes low-quality bases and adapter sequences.
* [HISAT2](https://daehwankimlab.github.io/hisat2/) is a fast and sensitive alignment tool for mapping next-generation sequencing reads to a reference genome.
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) is a statistical tool for analyzing count-based RNA-seq data, allowing differential gene expression analysis.
* [StringTie](https://ccb.jhu.edu/software/stringtie/) is a transcript assembler that reconstructs full-length transcripts and estimates their expression levels from RNA-seq data.




