# DECAFE: Differential Expression for Count Data Analysis with Functional Enrichment

**Authors:** [@AudreyBeaufils](https://github.com/AudreyBeaufils) & [@CamillePignolet](https://github.com/CamillePignolet) (GeNeHetX)  
**Date:** February 2024  
**Contact:** [audrey.beaufils@inserm.fr](mailto:audrey.beaufils@inserm.fr), [camille.pignolet@inserm.fr](mailto:camille.pignolet@inserm.fr)  

## Overview
DECAFE is a tool for analyzing RNA-Seq data, focusing on identifying differentially expressed genes and annotating them with functional information.<br><br>

  - **Principal Component Analysis (PCA)**<br>
  PCA is a statistical method used to explore relationships between variables in a dataset and reduce its dimensionality.<br><br>
  
  - **Heatmap Visualization**<br>
  Heatmaps display gene expression across samples to reveal patterns and clusters. They can be customized with sample annotations.<br><br>

  - **Differential Expression Analysis by DESeq2**<br>
  DESeq2 is a method used in genomics for identifying genes with significantly different expression levels between sample groups, particularly suited for RNA-Seq datasets.<br><br>

  - **Gene Set Enrichment Analysis (GSEA)**<br>
  GSEA is a bioinformatics method used to interpret functional genomics studies by comparing sets of functionally related genes for coordinated expression changes.<br><br>

  - **Over-Representation Analysis (ORA)**<br>
  ORA is used to determine whether predefined sets of genes (e.g., pathways or functional categories) are over-represented among a list of differentially expressed genes, helping to reveal enriched biological processes.<br><br>

  - **MCP-counter Analysis**<br>
  MCP-counter is a method for quantifying the abundance of eight immune and two stromal cell populations in tissues using transcriptomic data. It has been validated through mRNA mixtures and immunohistochemistry, outperforming previous methods. MCP-counter helps analyze immune infiltrates in healthy tissues and tumors, supporting 
  patient stratification and survival prediction in certain cancers.<br><br>


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
    devtools::install_github("nicolash2/gggsea")
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
  devtools::install_github("nicolash2/gggsea")
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
