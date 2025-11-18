# DECAFE: Differential Expression for Count Data Analysis with Functional Enrichment

**Authors:** [@AudreyBeaufils](https://github.com/AudreyBeaufils) & [@CamillePignolet](https://github.com/CamillePignolet) (GeNeHetX)  
**Date:** February 2024  
**Contact:** [audrey.beaufils@inserm.fr](mailto:audrey.beaufils@inserm.fr), [camille.pignolet@inserm.fr](mailto:camille.pignolet@inserm.fr)  

![](https://github.com/GeNeHetX/DECAFE/DECAFE/main/pres_readme.png)

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

- **REQUIRED** install R: [download](https://cran.r-project.org/)
- **REQUIRED** install the Rtools compiler : [download](https://cran.r-project.org/bin/windows/Rtools/)
- Not required, but if you want an IDE (integrated development environment), you can install Rstudio : [download](https://posit.co/download/rstudio-desktop/)
<br>
<br> 
   

## Installation 

### YOU HAVE AN INTERNET ACCES : 

1 - Install all the tools of DECAFE, run this command in a R terminal :
```R
    install.packages(c("shiny", "shinydashboard", "shinycssloaders", "plotly", "DT", "shinyBS", "devtools"))
    devtools::install_github('GeNeHetX/CancerRNASig')
    devtools::install_github("nicolash2/gggsea")
```

2- Then, run this command to open the interactive application :
```R
   shiny::runGitHub('DECAFE', 'GeNeHetX', subdir='DECAFE' ,ref='main')
```
ps : you can precise the version of DECAFE thanks to ref='', for example : ref='v.1.0.0'
<br> 
<br> 
   

### OTHER OPTION
If you are familiar with GIT, you can also download DECAFE with : 
```bash
  git clone https://github.com/GeNeHetX/DECAFE.git
```
Then run :
```R
    install.packages(c("shiny", "shinydashboard", "shinycssloaders", "plotly", "DT", "shinyBS", "devtools"))
    devtools::install_github('GeNeHetX/CancerRNASig')
    devtools::install_github("nicolash2/gggsea")
    setwd('<put the path of you git clone repository>/DECAFE/')
    shiny::runApp()
```
<br> 
<br> 

### YOU WANT TO USE DECAFE WITHOUT AN INTERNET ACCES :  

You need to download everything in advance : 

1- install the R packages required for the app : "shiny", "shinydashboard", "shinycssloaders", "plotly", "DT", "shinyBS", "devtools", "GeNeHetX/CancerRNASig", "nicolash2/gggsea", "BiocManager", "ggplot2", ,"reshape2", "factoextra", "FactoMineR", "devtools", "ggupset", "fgsea", "DESeq2", "ggpubr", "stringr", "ggrepel", "UpSetR", "ggdendro", "dendextend","gplots","svglite", "grid","gridExtra","ROTS", "circlize","scales".

2-Then you need to download DECAFE thanks to green button "<> Code" then "Download ZIP" ![](https://github.com/GeNeHetX/DECAFE/DECAFE/main/install.png)

3-Run the app without internet connexion :
```R
  setwd('<put the path where the zip are>/DECAFE/')
  shiny::runApp()
```
<br> 
<br> 
   

## DESIGN YOUR DATA TO UPLOAD THEM IN THE APP

### **Count matrix :**  
Upload the **FULL** RNA-seq count matrix as a **.tsv.gz** file (compress format).  
- Columns must be **Sample_IDs**.  
- Rows must be **genes**.  
- Upload the **entire** matrix — you don't need to subset it to include only the samples you plan to analyze.

---

### **Annotation txt file :**  
Provide a **.tsv** or **.txt** file using **tabulation separator**.  
- The **first column** must contain *only* the IDs of the samples you want to study. **Do not add a header to this first column.**  
- The following column(s) should contain the annotations you want to use for grouping samples in the analysis (DECAFE will use them to subset and create comparison groups).
