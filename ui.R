#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

suppressMessages(library(shiny))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinycssloaders))

library(plotly)

# Define UI for application that draws a histogram

  dashboardPage(skin = "purple",
  dashboardHeader(
            
    title = "DECAFE"

  
),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home",icon=shiny::icon('home')),
     
      menuItem("Differential Analysis", tabName = "tab1",icon=icon('arrows-rotate'))
      
    )
  ),
  dashboardBody(
    tags$style(HTML("
       
      .main-header .logo {
        font-family: 'Georgia', Times, 'Times New Roman', serif;
        font-weight: bold;
        font-size: 24px;
      }

      .box.box-solid.box-info>.box-header {
        color:#fff;
        background:#262686
      }
      .box.box-solid.box-info{
        border-bottom-color:#262686;
        border-left-color:#262686;
        border-right-color:#262686;
        border-top-color:#262686;
      }

      .box.box-solid.box-success>.box-header {
        color:#fff;
        background:#CDCDE6;
      }

      .box.box-solid.box-success{
        border-bottom-color:#CDCDE6;
        border-left-color:#CDCDE6;
        border-right-color:#CDCDE6;
        border-top-color:#CDCDE6;
      }
      
      .progress-bar {
        background-color: #CDCDE6;
      }
              
      .btn-default {
        background-color: #CDCDE6;
        color: #262686;
        border-color: #ddd;
      }

      .nav-tabs-custom>.nav-tabs>li.active {
        border-top-color: #262686;
      }
              
    ")),      tabItems(
  tabItem(tabName ="home",

          fluidRow(
             tags$div(style="text-align:center;",                        
                       h1(HTML("Welcome to <strong>D</strong>ifferential <strong>E</strong>xpression for <strong>C</strong>ount data <strong>A</strong>nalysis with <strong>F</strong>unctional <strong>E</strong>nrichment"))),br(),

          box(width = 12, title = h1('Presentation', icon('display')), status = 'success', solidHeader = TRUE, collapsible = TRUE,
              fluidRow(
                column(width = 12,
                       strong("Authors :"), "Audrey Beaufils & Camille Pignolet (GeNeHetX)", br(),
                       strong("Date :"), "February 2024", br(),
                       strong("Contact :"), "audrey1.beaufils@inserm.fr, camille.pignolet@inserm.fr", br(),
                       strong("GitHub :"), a("GeNeHetX", href = "https://github.com/GeNeHetX/"), br()))),

        
                       
        box(width = 12, title = h2('RNA-Seq analysis', icon('dna')), status = 'success', solidHeader = TRUE, collapsible = TRUE,
              fluidRow(
                column(width = 12,                 
                       "RNA-Seq analysis is a multi-step process used to extract meaningful information from raw sequencing data. It begins with quality control to assess 
                       the quality of the Fastq files, followed by preprocessing steps to trim adapters or low-quality bases. The reads are then aligned to a reference genome 
                       or transcriptome, and their abundance is quantified by counting the aligned reads. Differential expression analysis is performed to identify genes with 
                       significant expression changes between different conditions. These differentially expressed genes are annotated with functional information using databases 
                       like Gene Ontology (GO) or Kyoto Encyclopedia of Genes and Genomes (KEGG). Finally, functional enrichment analysis is conducted to identify overrepresented 
                       biological processes, pathways, or gene sets among the differentially expressed genes. The results are visualized and interpreted to gain insights into the biological significance of the gene expression changes.", align="justify"),
                        column(width = 12, align = "center",
                        imageOutput("rna_Image")))),
                     

        box(width = 12, title = h2('Principal Component Analysis (PCA)', icon('slack')), status = 'success', solidHeader = TRUE, collapsible = TRUE,
              fluidRow(
                column(width = 12,                 
                       "Principal Component Analysis (PCA) is a statistical method used to explore and summarize the relationships between variables in a multivariate dataset. 
                       Its main objective is to transform a set of correlated variables into a set of uncorrelated variables, called principal components. 
                       These components are ordered based on the importance of their contribution to the total variability of the data. 
                       PCA is often used to reduce the dimensionality of data by selecting the most significant principal components. 
                       This allows for visualizing the data in a lower-dimensional space, making interpretation and analysis easier.", align="justify"),
                        column(width = 12, align = "center",
                        imageOutput("pca_Image")))),
                     
                       
         box(width = 12, title = h2('Differential Expression Analysis by DESeq2', icon('dna')), status = 'success', solidHeader = TRUE, collapsible = TRUE,
              fluidRow(
                column(width = 12,         
                       "Differential Expression Analysis by DESeq2, also known as DESeq2, is a method of differential analysis primarily used in the field of genomics to 
                       identify genes whose expression varies significantly between different sample groups. This method is widely used in RNA sequencing (RNA-Seq) studies 
                       to detect genes that are differentially regulated in response to different treatments, experimental conditions, or biological states.
                       DESeq2 uses a statistical model based on the negative binomial distribution to model the variability of gene expression data and is particularly 
                       suited for RNA-Seq datasets with a small number of samples.", align="justify"), br(),
                       column(width = 12, align = "center",
                       imageOutput("anaDiff_Image")))),
                       

           box(width = 12, title = h2('Gene Set Enrichment Analysis (GSEA)', icon('circle-nodes')), status = 'success', solidHeader = TRUE, collapsible = TRUE,
              fluidRow(
                column(width = 12,                      
                       "Gene Set Enrichment Analysis (GSEA) is a widely used bioinformatics method for interpreting the results of functional genomics studies, 
                       such as DNA microarrays or RNA sequencing, focusing on functionally related sets of genes rather than individual genes. This method compares a 
                       predefined set of genes, often associated with a specific biological pathway or cellular function, with genes ranked by their differential expression 
                       across different experimental conditions. GSEA identifies sets of genes that show coordinated and significant changes in expression rather than individual genes, 
                       allowing for the detection of subtle biological alterations that may be overlooked by other differential analysis methods.", align="justify"), br(),
                        br(),
                        column(width = 12, align = "center",
                        imageOutput("gsea_Image")))),
                     

           box(width = 12, title = h1('What data do you need ?', icon('upload')), status = 'success', solidHeader = TRUE, collapsible = TRUE,
              fluidRow(
                column(width = 12,             
                       h3(""),
                        strong("Count matrix : "),
                       "You have to upload the entire matrix of RNAseq, the name of the columns is the Sample_ID and the row name is the GeneName. 
                       This matrix is created with the sequencing output. It is mapped onto a reference genome to identify which genes are present in the sample, then the number found for each gene is counted.", br(),
                        br(),
                        strong('Annotation file:'),
                       "The annotation file must contain only the samples to be studied , the first column is the Sample_ID and other columns are the annotations you want to use to create your groups in the analysis. Be careful to select only your samples of interest and your corresponding annotations. Any superfluous information will be taken into account as a comparison modality in the analysis.", br(),
                         br(),
                         column(width = 12, align = "center",
                         imageOutput("annot_Image")))))),
                 



            
          
        
          

    ),
    tabItem(tabName = "tab1",
        box(width = 12, status = 'info', title = h1("Settings", icon('cogs')), solidHeader = TRUE, collapsible=TRUE,
          fluidRow(
            column(width = 3,
              fileInput('file', 'Load Count Matrix'),
            ),
            column(width = 3,
              fileInput('annot-file',"Load the annot file")
            ),
            column(width = 3,

              uiOutput('cond_1'),
              uiOutput('cond_2')
             

            )

          )
        ),
        fluidRow(
        tabBox(
    tabPanel("Overwiew",height=1000,
        fluidRow(
        column(width=4,
        box(width=12, status = 'info', solidHeader = TRUE, title = h3("Table of Annotion", icon('table')),
          withSpinner( DT::DTOutput("tableAnnot"), type = 8, color = "#CDCDE6", size = 1)
        )),
        column(width=8,
                box(width=12, status = 'success', solidHeader = TRUE, title = h3("Overview of condition", icon('chart-simple')),
                br(), br(),
          withSpinner(plotOutput('barplot'), type = 8, color = "#CDCDE6", size = 1), br(), br(),
        )
        )

    )),
    tabPanel("Boxplot",
      fluidRow(column(width=3,uiOutput('geneTarget'))),
        fluidRow(

        column(width=6,
        box(width=12, status = 'success', solidHeader = TRUE, title = h3("Gene Target Histogram", icon('chart-simple')),
          withSpinner(plotOutput('hist_geneTarget'), type = 8, color = "#CDCDE6", size = 1)
        )
        ),
        column(width=6,
        box(width=12, status = 'success', solidHeader = TRUE, title = h3("Gene Target Boxplot", icon('chart-simple')),
            
          withSpinner(plotOutput('bp_geneTarget'), type = 8, color = "#CDCDE6", size = 1)
        )
        )

      
    )),
    tabPanel("PCA",
        fluidRow(
        box(width=12,status='success',title = h2('Graph of PCA',icon('chart-simple')),solidHeader = TRUE,
            column(width=3,
        numericInput("nb_gene", "Number of most variable Gene", min = 1, step = 1, value = 1000),
             selectInput("dim1", label = "Choose your first PCA dimension",
                  choices = list("Dim1" = 1, "Dim2" = 2,"Dim3" = 3, "Dim4" =4, "Dim5" = 5),selected = 1),
              selectInput("dim2", label = "Choose your secondPCA dimension",
                  choices = list( "Dim2" = 2,"Dim3" = 3, "Dim4" =4, "Dim5" = 5),selected = 2),
            
            #selectInput("color.pca", label = "Color by ",
            #      choices = list("Dim1" = "1", "Dim2" = "2","Dim3" = "3", "Dim4" ="4", "Dim5" = "5"),selected = "2")
            #
            ),
            column(width=9, 
          withSpinner(plotlyOutput("pcaVsd"), type = 8, color = "#CDCDE6", size = 1)
        )
        ),

            column(width=6,
                box(width=12, status = 'info', solidHeader = TRUE, title = h3("Gene most contributed to first dim", icon('table')),
          withSpinner( DT::dataTableOutput("geneTop1"), type = 8, color = "#CDCDE6", size = 1)
        ),
                 box(width=12, status = 'info', solidHeader = TRUE, title = h3("Sample most contributed to first dim", icon('table')),
          withSpinner( DT::dataTableOutput("sampleTop1"), type = 8, color = "#CDCDE6", size = 1)
        )
        
            ),
            column(width=6,
                box(width=12, status = 'info', solidHeader = TRUE, title = h3("Gene most contributed to second dim", icon('table')),
          withSpinner( DT::dataTableOutput("geneTop2"), type = 8, color = "#CDCDE6", size = 1)
        ),
                 box(width=12, status = 'info', solidHeader = TRUE, title = h3("Sample most contributed to second dim", icon('table')),
          withSpinner( DT::dataTableOutput("sampleTop2"), type = 8, color = "#CDCDE6", size = 1)
        )

            )
             
          
        

      
      )),
    tabPanel("Differential Analysis",
        fluidRow(
        column(width=12,
      box(width=NULL,status='success',solidHeader=TRUE,title=h1("Volcano Plot",icon('chart-simple')),
        fluidRow(
          column(width=4,
            sliderInput("zero_threshold",label = "Frequency of zero per gene to remove",
                        min = 0, max =0.8, value = 0.5,step=0.1),
              sliderInput("ts_padj",label = "p-Value cutoff from output ",
                        min = 0, max =0.1, value = 0.01,step=0.01)
            
          ),
          column(width=8,
            withSpinner(plotlyOutput("volcano"), type = 8, color = "#CDCDE6", size = 1)
          )
        )
      )
    )),
    fluidRow(
    column(width=12,
      box(width=NULL,status='info',title = h1('Table of results',icon('table')),solidHeader = TRUE, 
        withSpinner(DT::dataTableOutput("degTable"), type = 8, color = "#CDCDE6", size = 1)
      )
    )



    )),
    tabPanel("GSEA",
      fluidRow(
        column(width=12,
        box(width=NULL,status='info',title = h1('Table of GSEA results',icon('table')),solidHeader = TRUE, 
        withSpinner(DT::dataTableOutput("gsea"), type = 8, color = "#CDCDE6", size = 1)
      ))

    )),id="tabBox",

    width = 12

    )
      
)))))
