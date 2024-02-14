# DECAFE: Differential Expression for Count Data Analysis with Functional Enrichment

**Authors:** [@AudreyBeaufils](https://github.com/AudreyBeaufils) & [@CamillePignolet](https://github.com/CamillePignolet) (GeNeHetX)  
**Date:** February 2024  
**Contact:** [audrey1.beaufils@inserm.fr](mailto:audrey1.beaufils@inserm.fr), [camille.pignolet@inserm.fr](mailto:camille.pignolet@inserm.fr)  

## Overview
DECAFE is a tool for analyzing RNA-Seq data, focusing on identifying differentially expressed genes and annotating them with functional information.<br><br>

  <ins>Principal Component Analysis (PCA)<ins><br>
  
  PCA is a statistical method used to explore relationships between variables in a dataset and reduce its dimensionality.<br><br>

  <ins>Differential Expression Analysis by DESeq2<ins><br>
  
  DESeq2 is a method used in genomics for identifying genes with significantly different expression levels between sample groups, particularly suited for RNA-Seq datasets.<br><br>

  <ins>Gene Set Enrichment Analysis (GSEA)<ins><br>
  
  GSEA is a bioinformatics method used to interpret functional genomics studies by comparing sets of functionally related genes for coordinated expression changes.<br><br>


## Upload data
**Count matrix:** You need to upload the entire RNA-Seq count matrix with Sample_ID as column names and GeneName as row names.<br>

**Annotation file:** it should contain only the samples to be studied, with Sample_ID in the first column and other columns for annotations used to create groups in the analysis.<br>

___________________________________
### Prerequisites : 
DECAFE requires the R language (at least version 4.0).<br>

If R is installed, you can launch the application directly via a command terminal or work on Rstudio

**install R:** [downlaod here](https://cran.r-project.org/)
you can find **Rstudio** here : [download](https://posit.co/download/rstudio-desktop/)

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
