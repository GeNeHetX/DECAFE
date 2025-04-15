# DECAFE: Differential Expression for Count Data Analysis with Functional Enrichment

**Authors:** [@AudreyBeaufils](https://github.com/AudreyBeaufils) & [@CamillePignolet](https://github.com/CamillePignolet) (GeNeHetX)  
**Date:** February 2024  
**Contact:** [audrey1.beaufils@inserm.fr](mailto:audrey.beaufils@inserm.fr), [camille.pignolet@inserm.fr](mailto:camille.pignolet@inserm.fr)  

## Overview
DECAFE is a tool for analyzing RNA-Seq data, focusing on identifying differentially expressed genes and annotating them with functional information.<br><br>

  - **Principal Component Analysis (PCA)**<br>
  PCA is a statistical method used to explore relationships between variables in a dataset and reduce its dimensionality.<br><br>

  - **Differential Expression Analysis by DESeq2**<br>
  DESeq2 is a method used in genomics for identifying genes with significantly different expression levels between sample groups, particularly suited for RNA-Seq datasets.<br><br>

  - **Gene Set Enrichment Analysis (GSEA)**<br>
  GSEA is a bioinformatics method used to interpret functional genomics studies by comparing sets of functionally related genes for coordinated expression changes.<br><br>



## Prerequisites : 
DECAFE requires the R language (at least version 4.0).<br>
If R is installed, you can launch the application directly via a command terminal or work on Rstudio.

- install R: [download](https://cran.r-project.org/)

- You can find Rstudio here : [download](https://posit.co/download/rstudio-desktop/)
<br>


## Installation 

With internet : 

1 - First-time use DECAFE, run this command in a R terminal
```R
    install.packages(c("shiny", "shinydashboard", "shinycssloaders", "plotly", "DT", "shinyBS", "devtools"))
    devtools::install_github('GeNeHetX/CancerRNASig')
```

2- Then, run this command :
```R
   shiny::runGitHub('DECAFE', 'GeNeHetX', subdir='DECAFE' ,ref='main')
```
ps : you can precise the version thanks to ref='', for example : ref='v.1.0.0'
___________________________________________________

Without internet, (use just to download and run the following commands without a connection) : 

1- If you are a git user, clone the DECAFE folder, otherwise download the DECAFE code zip via the green "<>Code" button.

```bash
  git clone https://github.com/GeNeHetX/DECAFE.git
```

2- Open an R terminal or Rstudio where the DECAFE codes are stored
     
- First-time use DECAFE
```R
  install.packages(c("shiny", "shinydashboard", "shinycssloaders", "plotly", "DT", "shinyBS", "devtools"))
  devtools::install_github('GeNeHetX/CancerRNASig')
  shiny::runApp()
```
or

- Use DECAFE 
```R
  shiny::runApp()
```
<br>

## Upload data

**Count matrix:**  It must be a .tsv file. You need to upload the entire RNA-Seq count matrix with Sample_ID as column names and genes as row names.<br>

**Annotation file:** It must be a .tsv file or .txt file with tabulation separator. It should contain only the samples to be studied, with Sample_ID in the first column and other columns for annotations used to create groups in the analysis.
