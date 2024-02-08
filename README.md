# DECAFE: Differential Expression for Count Data Analysis with Functional Enrichment

**Authors:** [@AudreyBeaufils](https://github.com/AudreyBeaufils) & [@CamillePignolet](https://github.com/CamillePignolet) (GeNeHetX)  
**Date:** February 2024  
**Contact:** [audrey1.beaufils@inserm.fr](mailto:audrey1.beaufils@inserm.fr), [camille.pignolet@inserm.fr](mailto:camille.pignolet@inserm.fr)  

## Overview
DECAFE is a tool for analyzing RNA-Seq data, focusing on identifying differentially expressed genes and annotating them with functional information.

  **Principal Component Analysis (PCA)**
  Principal Component Analysis (PCA) is a statistical method used to explore relationships between variables in a dataset and reduce its dimensionality.

  **Differential Expression Analysis by DESeq2**
  DESeq2 is a method used in genomics for identifying genes with significantly different expression levels between sample groups, particularly suited for RNA-Seq datasets.

  **Gene Set Enrichment Analysis (GSEA)**
  GSEA is a bioinformatics method used to interpret functional genomics studies by comparing sets of functionally related genes for coordinated expression changes.


## Upload data
**Count matrix:** You need to upload the entire RNA-Seq count matrix with Sample_ID as column names and GeneName as row names.
**Annotation file:** The annotation file should contain only the samples to be studied, with Sample_ID in the first column and other columns for annotations used to create groups in the analysis.


### DECAFE GitHub Download 
```bash
  git clone https://github.com/GeNeHetX/DECAFE.git

```
### First-time use
```bash
cd DECAFE/
Rscript -e 'install.packages(c("shiny", "shinydashboard", "shinycssloaders", "plotly", "DT")); shiny::runApp()'
```

#### Run DECAFE 
```bash
  cd DECAFE/
  Rscript -e 'shiny::runApp()'
```
