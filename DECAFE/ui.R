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
library(htmltools)
library(shinyBS)


# Define UI for application that draws a histogram

  dashboardPage(skin = "purple",
  dashboardHeader(
            
    title = "DECAFE",titleWidth=300
 
    #textOutput('title')
    

  
),
  dashboardSidebar(
    width=300,
    sidebarMenu(
      selectInput('lcms', 'Choose your data type', choices = list(RNASeq='rna', "LC-MS/MS"='lcms')),

      menuItem("Home", tabName = "home",icon=shiny::icon('home')),
     
      menuItem("Analysis", tabName = "tab1",icon=icon('arrows-rotate')),

      menuItem("Tools", tabName = "tab2",icon=icon('screwdriver-wrench'))
      
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

      .nav-tabs-custom>.nav-tabs {
      margin: 0;
      border-bottom-color: #f4f4f4;
      border-top-right-radius: 3px;
      border-top-left-radius: 3px;
      background: #ddd;
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
        
              
    ")),uiOutput('theme'),

    htmlDependency(
      "font-awesome", "5.3.1", "www/shared/fontawesome", package = "shiny",
               stylesheet = c("css/all.min.css", "css/v4-shims.min.css")
    ),
    tags$script(HTML('Shiny.addCustomMessageHandler("changetitle", function(x) {$(".logo").html(x);document.title=x});')),

    
    tags$head(tags$script('
      

      // Define function to set height of "map" and "map_container"
      setHeight = function() {
        var window_height = $(window).height();
        var header_height = $(".main-header").height();

        var boxHeight = (window_height - header_height) * 0.6;

        $(".map_container").height(boxHeight);
        $("#pcaVST").height(boxHeight - 20);
        $("#volcano").height(boxHeight - 20);
        $("#treePlot").height(boxHeight - 10);
      };

      // Set input$box_height when the connection is established
      $(document).on("shiny:connected", function(event) {
        setHeight();
      });

      // Refresh the box height on every window resize event    
      $(window).on("resize", function(){
        setHeight();

      $(window).resize(function(event){
      var w = $(this).width();
      var h = $(this).height();
      var obj = {width: w, height: h};
      Shiny.onInputChange("pltChange", obj);
    
      });

      

      ')),

       tabItems(
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
                       strong("GitHub :"), a("GeNeHetX", href = "https://github.com/GeNeHetX/"), br(), br(),

                      "RNA-Seq analysis is a multi-step process used to extract meaningful information from raw sequencing data:", br(), br(),
                      column(width = 12, align = "center",
                        imageOutput("rna_Image")),
                       
                      strong("Quality Control (QC):"),"This initial step is akin to performing a thorough inspection of the raw sequencing data, known as Fastq files. 
                      QC procedures assess various parameters such as sequencing depth, read quality scores, and potential biases to ensure the reliability and integrity of the data. 
                      Any anomalies or discrepancies detected during this phase are flagged for subsequent correction.", br(),br(),
                      
                      strong("Preprocessing:"),"these steps are implemented to prepare the data for downstream analyses. This involves trimming adapter sequences 
                      and removing low-quality bases that may skew the accuracy of subsequent analyses. By refining the data in this manner, 
                      researchers mitigate the risk of false-positive results and enhance the robustness of the analysis pipeline.", br(),br(),

                      strong("Alignment:"),"the next step involves aligning the sequencing reads to a reference genome or transcriptome. This alignment process 
                      involves mapping each read to its corresponding genomic location, facilitating the identification of gene transcripts and their regulatory elements.
                      Accurate alignment is crucial for accurately quantifying gene expression levels and discerning subtle variations across different experimental conditions.", br(), br(),

                      strong("Quantification:"),"Once the reads have been aligned, their abundance is quantified by enumerating the number of reads that map to each gene or transcript. 
                      This quantitative measure serves as a proxy for gene expression levels, allowing researchers to assess the relative abundance of transcripts across samples. 
                      By quantifying gene expression in this manner, researchers can identify genes that exhibit significant changes in expression under different experimental conditions.", br(),br(),

                      strong("Analysis:"),"Differential expression analysis is help to identify genes exhibiting statistically significant changes in expression levels between different experimental conditions.
                      Following the identification of differentially expressed genes, functional annotation is performed to elucidate the biological roles and significance of these genes. This involves annotating genes with
                      functional information derived from specialized databases such as Gene Ontology (GO) or Kyoto Encyclopedia of Genes and Genomes (KEGG). By annotating differentially expressed genes 
                      with functional metadata, researchers gain deeper insights into the molecular pathways and biological processes implicated in the observed expression changes. In addition, functional enrichment analysis
                       is conducted to identify overrepresented biological processes, pathways, or gene sets among the differentially expressed genes.", br(),br(),

                      strong("Visualization and Interpretation:"),"The final step of RNA-Seq analysis involves visualizing and interpreting the results. This may involve generating 
                      visual representations such as heatmaps, volcano plots, or pathway diagrams to illustrate the patterns and relationships within the data. Through 
                      careful interpretation of the results in the context of existing biological knowledge, researchers can derive valuable insights into the underlying molecular 
                      mechanisms driving the observed expression changes.", br(),br(), align="justify"))),

        
                       
        # box(width = 12, title = h2('RNA-Seq analysis', icon('dna')), status = 'success', solidHeader = TRUE, collapsible = TRUE,
        #       fluidRow(
        #         column(width = 12,                 
        #                "RNA-Seq analysis is a multi-step process used to extract meaningful information from raw sequencing data. It begins with quality control to assess 
        #                the quality of the Fastq files, followed by preprocessing steps to trim adapters or low-quality bases. The reads are then aligned to a reference genome 
        #                or transcriptome, and their abundance is quantified by counting the aligned reads. Differential expression analysis is performed to identify genes with 
        #                significant expression changes between different conditions. These differentially expressed genes are annotated with functional information using databases 
        #                like Gene Ontology (GO) or Kyoto Encyclopedia of Genes and Genomes (KEGG). Finally, functional enrichment analysis is conducted to identify overrepresented 
        #                biological processes, pathways, or gene sets among the differentially expressed genes. The results are visualized and interpreted to gain insights into the biological significance of the gene expression changes.", align="justify"),br(),
        #                 column(width = 12, align = "center",
        #                 imageOutput("rna_Image")))),
                     
        box(width = 12, title = h1('What data do you need ?', icon('upload')), status = 'success', solidHeader = TRUE, collapsible = TRUE,
              fluidRow(
                column(width = 12,             
                       h3(""),
                        strong("Count matrix : "),
                       "You have to upload the RNA-Seq count matrix with Sample_ID as column names and genes as row names, it must be a .tsv file.
                       This matrix is created with the sequencing output. It is mapped onto a reference genome to identify which genes are present in the sample, then the number found for each gene is counted.", br(),
                        br(),
                        strong('Annotation file:'),
                       "It must be a .tsv file or .txt file with tabulation separator. The annotation file must contain only the samples to be studied, the first column is the Sample_ID and other columns are the annotations you want to use to create your groups in the analysis. Be careful to select only your samples of interest and your corresponding annotations. Any superfluous information will be taken into account as a comparison modality in the analysis.", br(),
                         br(),br(),
                         column(width = 12, align = "center",
                         imageOutput("annot_Image"))))), 


        box(width = 12, title = h2('Principal Component Analysis (PCA)', icon('slack')), status = 'success', solidHeader = TRUE, collapsible = TRUE,
              fluidRow(
                column(width = 12,                 
                       "Principal Component Analysis (PCA) is a statistical method used to explore and summarize the relationships between variables in a multivariate dataset. 
                       Its main objective is to transform a set of correlated variables into a set of uncorrelated variables, called principal components. 
                       These components are ordered based on the importance of their contribution to the total variability of the data. 
                       PCA is often used to reduce the dimensionality of data by selecting the most significant principal components. 
                       This allows for visualizing the data in a lower-dimensional space, making interpretation and analysis easier.", br(),br(),
                            
                      strong("Tips:"),"PCA reduces the dimensions of your data. These compact dimensions will group your samples by similarity. The idea here is to see if your samples 
                      separate according to a new dimension. If so, you'll have the contribution value of each sample and gene in the new dimension and know which ones contribute most to your dimension.", align="justify"),br(),br(),
                        column(width = 12, align = "center",
                        imageOutput("pca_Image"))), 

                        tags$details(
                            tags$summary(strong("More information")),

                      column(width = 12, 
                            "PCA is based on fundamental mathematical concepts such as the eigenvalue and eigenvector decomposition of a covariance matrix. 
                            Firstly, the covariance matrix is calculated from the input data, representing the correlation relationships between the different variables. 
                            Next, the eigenvalue decomposition of this matrix yields the principal components, which are the directions in which the data vary most. 
                            The eigenvectors associated with these eigenvalues determine the linear combinations of the original variables that form the principal components. 
                            Finally, the data are projected into the space defined by these principal components, thus reducing dimensionality while preserving maximum variance.",br(), 
                            HTML("<p>The decomposition of the covariance matrix Σ is done as follows:</p>
                         <p>Σ = (1/(n-1)) * (X - X̄)^T * (X - X̄)</p>
                         <p>Where X is the centered data matrix, X̄ is the vector of means of each variable, and n is the number of observations. The eigenvalues λ_i and eigenvectors v_i are obtained by solving the characteristic equation:</p>
                         <p>Σv_i = λ_i * v_i</p>
                         <p>The principal components are then calculated as linear combinations of the original variables:</p>
                         <p>Y = X * V</p>
                         <p>Where Y is the matrix of principal components and V is the matrix whose columns are the corresponding eigenvectors v_i.</p>
                         "), align="justify")) 
                                            
                        ),
                     
                       
         box(width = 12, title = h2('Differential Expression Analysis by DESeq2', icon('dna')), status = 'success', solidHeader = TRUE, collapsible = TRUE,
              fluidRow(
                column(width = 12,         
                       "Differential Expression Analysis by DESeq2, also known as DESeq2, is a method of differential analysis primarily used in the field of genomics to 
                       identify genes whose expression varies significantly between different sample groups. This method is widely used in RNA sequencing (RNA-Seq) studies 
                       to detect genes that are differentially regulated in response to different treatments, experimental conditions, or biological states.
                       DESeq2 uses a statistical model based on the negative binomial distribution to model the variability of gene expression data and is particularly 
                       suited for RNA-Seq datasets with a small number of samples.", align="justify"), br(),br(),
                       column(width = 12, align = "center",
                       imageOutput("anaDiff_Image"))),
                       
                       tags$details(
                            tags$summary(strong("More information")),

                      column(width = 12, 
                            "Differential analysis, as performed by DESeq2, involves the use of advanced statistical models based on the negative binomial distribution to model 
                            the variability of gene expression data. These models take into account several factors, such as covariates and sampling bias, 
                            to identify genes whose expression varies significantly between sample groups. Statistical calculations involve fitting these models to each individual 
                            gene to estimate expression differences and assess their statistical significance. The test scores obtained enable genes to be ranked according to their 
                            degree of variation in expression between experimental conditions.",br(), 
                            HTML("<p>In the case of DESeq2, modeling gene expression variability is based on a negative binomial distribution. For each gene, a regression model is fitted as follows:</p>
                         <p>Y_i ∼ NegBin(size = μ_i, prob = θ)</p>
                         <p>Where Y_i represents the number of reads associated with gene i, μ_i is the average expression level of gene i, and θ is a global dispersion parameter. Differences in expression between groups are estimated from this model, and a Wald test is used to assess their significance.</p>
                         "), align="justify")) 

                       ),
                       

           box(width = 12, title = h2('Gene Set Enrichment Analysis (GSEA)', icon('circle-nodes')), status = 'success', solidHeader = TRUE, collapsible = TRUE,
              fluidRow(
                column(width = 12,                      
                       "Gene Set Enrichment Analysis (GSEA) is a widely used bioinformatics method for interpreting the results of functional genomics studies, 
                       such as DNA microarrays or RNA sequencing, focusing on functionally related sets of genes rather than individual genes. This method compares a 
                       predefined set of genes, often associated with a specific biological pathway or cellular function, with genes ranked by their differential expression 
                       across different experimental conditions. GSEA identifies sets of genes that show coordinated and significant changes in expression rather than individual genes, 
                       allowing for the detection of subtle biological alterations that may be overlooked by other differential analysis methods.", align="justify"), br(),
                        br(),br(),
                        column(width = 12, align = "center",
                        imageOutput("gsea_Image"))), 
                        
                        tags$details(
                            tags$summary(strong("More information")),

                      column(width = 12, 
                            "GSEA uses rank statistic calculations to assess the enrichment of gene sets associated with specific biological pathways or cellular functions. First, genes are 
                            ranked according to their difference in expression between experimental groups. In our case, they are ranked by their results in differential analysis. 
                            Next, a richness score is calculated for each predefined gene set, reflecting the collective tendency of genes in that set to be positioned at the top or bottom 
                            of the ranked list. Finally, permutation tests are performed to assess the significance of these richness scores and identify gene sets whose expression is 
                            significantly enriched under the experimental conditions studied.",br(), 
                            HTML("<p>The Enrichment Score (ES) for a gene set is calculated by evaluating the overall tendency of genes in that set to be positioned at the top or bottom of the ranked list of genes based on their expression differences between groups. The calculation is performed as follows:</p>
                         <p>ES = max(Σ_{i=1}^{k} (R_i / N_1), Σ_{i=1}^{k} (-S_i / N_2))</p>
                         <p>Where R_i and S_i represent the cumulative positive and negative expression differences up to gene i, respectively, in the ranked list. N_1 and N_2 are the total numbers of genes in the two compared groups.</p>
                         "), align="justify")) 

                        )
                     

           ),
                 



            
          
        
          

    ),
    tabItem(tabName = "tab1",
        box(width = 12, status = 'info', title = h1("Settings", icon('cogs')), solidHeader = TRUE, collapsible=TRUE,
          fluidRow(
            column(width = 3,
              fileInput('file', 'Load Count Matrix .tsv'),
              fileInput('annot-file',"Load the annot file .tsv")
              
            ),
            column(width = 3,
              selectInput('org', 'Choose your species', choices = list(Human='hs', Mouse='mm')),#, Other='oth'))
              conditionalPanel("input.lcms=='rna'",
              radioButtons('coding', label= div('Use only coding genes',icon('circle-info')), choices = list(YES=TRUE, NO=FALSE), inline=TRUE, selected = FALSE),
              bsTooltip("coding",title="Remove all non-coding genes from analysis"),
              radioButtons('sex', label=div('Sex-independent analysis',icon('circle-info')), choices = list(YES=TRUE, NO=FALSE), inline=TRUE, selected = FALSE),
              bsTooltip('sex',title="Remove all gene in chrs X/Y from analysis"),
              conditionalPanel(condition="input.sex=='TRUE'", fileInput('sexAnnot','Load sex information'))),
              conditionalPanel("input.lcms=='lcms'",
                radioButtons('normalized', 'Matrix already normalized', choices = list(YES=TRUE, NO=FALSE), inline=TRUE, selected = TRUE),
                radioButtons('adlc', 'Choose type of differential analysis', choices = list("ROTS"="ROTS", "t.test"="t.test"), inline=TRUE, selected = "ROTS")
                )
              
              # conditionalPanel("input.org == 'oth'", fileInput('genefile', 'Load Gene Annotation'))
            ),
            
            column(width = 3,
              numericInput('nb_thread', 'Provide the number of CPU to use', 1, min = 1, max = (parallel::detectCores()-1)), 
              p(icon('circle-info'),paste0(" You have ",parallel::detectCores(), " CPU")),
              radioButtons('format', 'Choose saving format', choices = list("PNG/PDF"='png', "SVG (modifiable)"='svg'), inline=TRUE, selected = 'svg')
              
            ),
            column(width =3,
              uiOutput('cond1'),
              uiOutput('cond2'),
              uiOutput('control')
             
                   
            )

          )
        ),
        fluidRow(
        tabBox(
    tabPanel("Overview", height = 1000,
         fluidRow(
           column(width = 6,
                  box(width = 12, status = 'info', solidHeader = TRUE, title = h3("Check : Sample_ID and Annotations", icon('file')),
                      withSpinner(DT::DTOutput("tableAnnot"), type = 8, color = "#CDCDE6", size = 1), style = 'overflow-x: scroll'
                  )
           ),
          column(width = 6,
                  box(width = 12, status = 'info', solidHeader = TRUE, title = h3("Check : Sample_ID and Genes", icon('dna')),
                      withSpinner(DT::DTOutput("counthead"), type = 8, color = "#CDCDE6", size = 1)
                  )
           ),
           fluidRow(
           column(width = 12,
                  box(width = 12, status = 'success', solidHeader = TRUE, title = h3("Overview of Condition", icon('chart-simple')),
                      downloadButton("downloadUpsetPlot", "Download UpsetPlot", icon('download')), br(), br(),
                      withSpinner(plotOutput('upsetPlot',  width = 1500, height = 800, inline=F), type = 8, color = "#CDCDE6", size = 1)
                  )
           ))
         )
    ),
    tabPanel("Heatmap", height = 1000,
          fluidRow(
        box(class = "map_container",width=12,status='success',title = h2('Heatmap of normalize count',icon('chart-simple')),solidHeader = TRUE,
            column(width=2,
        numericInput("nb_gene_heat", "Number of most variable Gene", min = 10, step = 1, max = 1000,value = 1000),
        uiOutput('nbGene2'),
             selectInput("data_heat", label = "Choose sample to visualize",
                  choices = list("All" = "all", "Just the two conditions" = "cond"),selected = "cond"),
              
              downloadButton("downloadHeatmap", "Download heatmap", icon('download')), br(), br(),
          
            ),
           column(width=1),
            column(width=9, 
          withSpinner(plotOutput("heatMap", width = "95%", height=1200), type = 8, color = "#CDCDE6", size = 1), style = 'display:block;width:100%'
        ),
        )),
        
    ),
    tabPanel("PCA",
        fluidRow(
        box(class = "map_container",width=12,status='success',title = h2('Graph of PCA',icon('chart-simple')),solidHeader = TRUE,
            column(width=3,
        numericInput("nb_gene", "Number of most variable", min = 1, step = 1, value = 1000),
        uiOutput('nbGene'),
             selectInput("dim1", label = "Choose your first PCA dimension",
                  choices = list("Dim1" = 1, "Dim2" = 2,"Dim3" = 3, "Dim4" =4, "Dim5" = 5),selected = 1),
              # p(icon('circle-info'),paste0(" You have ",dim_matrix, " genes in your matrix")),
              selectInput("dim2", label = "Choose your second PCA dimension",
                  choices = list( "Dim2" = 2,"Dim3" = 3, "Dim4" =4, "Dim5" = 5),selected = 2),
              downloadButton("downloadPCAPlot", "Download PCA plot with all conditions", icon('download')),
          
            ),
            column(width=1),
            column(width=8, 
          withSpinner(plotlyOutput("pcaVST"), type = 8, color = "#CDCDE6", size = 1)
        )
        )),
        fluidRow(
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
      box(class = "map_container", width=NULL,status='success',solidHeader=TRUE,title=h1("Volcano Plot",icon('chart-simple')),
        fluidRow(
          column(width=3,
            sliderInput("zero_threshold",label = "Frequency of zero per gene to remove",
                        min = 0, max =0.8, value = 0.5,step=0.1),
            sliderInput("ts_padj",label = "p-Value cutoff from output ",
                        min = 0, max =0.1, value = 0.05,step=0.01),br(),
            downloadButton("downloadVolcanoPlot", "Download VolcanoPlot with gene label", icon('download')),
              uiOutput('sens1')
            
          ),
          column(width=1),
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


    tabPanel("Boxplot",
      fluidRow(column(width=3,uiOutput('geneTarget')),
              column(width=7),
              column(width = 1, downloadButton("downloadboxplot", "Download Plots", icon('download')))),br(),
        fluidRow(
        column(width=6,
        box(width=12, status = 'success', solidHeader = TRUE, title = h3("Gene Target DensityPlot", icon('chart-simple')),
          uiOutput('sampleTarget2'),
          withSpinner(plotOutput('densityPlotgene'), type = 8, color = "#CDCDE6", size = 1)
        )
        ),
        column(width=6,
        box(width=12, status = 'success', solidHeader = TRUE, title = h3("Gene Target Boxplot", icon('chart-simple')),
          withSpinner(plotOutput('bpGeneTarget'), type = 8, color = "#CDCDE6", size = 1)
        )
        ),
        conditionalPanel("input.lcms=='lcms'",
          column(width=6,
        box(width=12, status = 'success', solidHeader = TRUE, title = h3("Gene Target Boxplot (pval t.test)", icon('chart-simple')),
          withSpinner(plotOutput('bpGeneTarget_Ttest'), type = 8, color = "#CDCDE6", size = 1)
        )
        ),


          )

      
    )),


    tabPanel("GSEA",
      fluidRow(
        column(width=7,
        imageOutput('path_Image')),
        
        column(width = 3,
        uiOutput('dbpath')),

        column(width = 2,
        radioButtons('gsea_level', label= div('Choose type of GSEA',icon('circle-info')), choices = list("Simple Level"=TRUE, "Multi Level"=FALSE), inline=TRUE, selected = TRUE),
              bsTooltip("gsea_level",title="Simple level is faster but less accurate while Multi level is slower and more acuurate (resulted more significate pathways) "),
        actionButton('gogsea', label='Run GSEA', icon('play')))
      ),
      fluidRow(
        column(width=12,br(),br(),br(),br(),
        box(width=12,status='info',title = h1('Table of GSEA results',icon('table')),solidHeader = TRUE, 
        actionButton('browsebutton', 'More info about selected pathway', icon('globe')),
          uiOutput('sens2'),
        textOutput('comment'),br(),
        withSpinner(DT::dataTableOutput("gsea"), type = 8, color = "#CDCDE6", size = 1)
      ))),
      fluidRow(
      column(width=12,
        box(width=6,status='success',title = h1('Dotplot',icon('chart-simple')),solidHeader = TRUE, collapsible=TRUE,
          fluidRow(column(width = 8, sliderInput("nbDotplot",label = "Number of pathways to display", min = 1, max =20, value = 10,step=1))),
         fluidRow( column(width=12,withSpinner(plotlyOutput("dotplot",height=1000), type = 8, color = "#CDCDE6", size = 1)))),
        box(width=6,status='success',title = h1('GSEA plot',icon('chart-simple')),solidHeader = TRUE, collapsible=TRUE,
          fluidRow(column(width = 8, sliderInput("nbGseaPlot",label = "Number of pathways to display", min = 1, max =20, value = 10,step=1))),
          fluidRow(column(width=12,withSpinner(plotOutput("gseaPlot", height=1000,width="90%"), type = 8, color = "#CDCDE6", size = 1))))
      )
      ),
      fluidRow(
      column(width=12,
        box(width=12,status='success',title = h1('Tree Pathways',icon('square-poll-horizontal')),solidHeader = TRUE, 
          column(width = 8, sliderInput("k",label = "Number of different groups", min = 1, max =20, value = 5,step=1)),
          withSpinner(plotOutput("treePlot", inline=F, width = 1200, height=1500), type = 8, color = "#CDCDE6", size = 1), style = 'display:block;width:100%;overflow-y: scroll')
      )
      )



      ),tabPanel("MCPcounter",
      actionButton('gomcp',label="Run MCPcounter",icon('play')),br(),br(),br(),
      fluidRow(column(width=12,
      box(width=12,status='success',title = h1('All immune cell famillies',icon('chart-simple')),solidHeader = TRUE, 
      withSpinner(
      plotOutput("allboxMCP"#, inline=F, width = 1500
        , height=800
        ), 
      type = 8, color = "#CDCDE6", size = 1)))),

      br(),br(),
      fluidRow(column(width=12,
      box(width=12,status='info',title = h1('McpCounter projection',icon('table')),solidHeader = TRUE, 
        withSpinner(DT::dataTableOutput("mcptable"), type = 8, color = "#CDCDE6", size = 1)
      ))),
       br(),br(),
       fluidRow(column(width=12,
      box(
        width=6,status='success',title = h1('Choose one immune cell familly',icon('chart-simple')),solidHeader = TRUE, 
        uiOutput('mcpCond'),
        # selectInput("mcpPath", label = "Choose type to cell to plot",
        #           choices = c(
        #             "T cells"="Tcells",
        #             "CD8Tcells"="CD8Tcells",
        #             "Cytotox.lymph"="Cytotox.lymph",
        #             "NK" = "NK",
        #             "B.lineage"="B.lineage",
        #             "Mono.lineage"="Mono.lineage",    
        #           "Myeloid.dendritic" ="Myeloid.dendritic",
        #             "Neutrophils"="Neutrophils",
        #             "Endothelial"="Endothelial" ,
        #             "Fibroblasts"="Fibroblasts"
        #             ),selected=1),
        withSpinner(
          plotOutput('boxMCP'), 
          type = 8, color = "#CDCDE6", size = 1)
        ),
       box(
        width=6,status='success',title = h1('Choose one sample',icon('chart-simple')),solidHeader = TRUE, 
        uiOutput('sampleTarget'),
      withSpinner(plotOutput("densityPlot"), type = 8, color = "#CDCDE6", size = 1)

      )),


      )),




    id="tabBox",

    width = 12

    )
      
)),
  
  tabItem(tabName ="tab2",
  fluidRow(
            column(width = 12,
            box(width = 12, status = 'info', title = h1("Design your own count matrix", icon('table')), solidHeader = TRUE, collapsible=TRUE,
                fluidRow(
                    box(width = 4, title = h2('Normalization', icon('table')), collapsible = TRUE,
                      fluidRow(column(width=12,
                        radioButtons("normalize","Choose type of normalization",c("None"="no","Variance Stabilizing Transformation (VST)"="vst","DESeq2 normalization"="deseq"))

                      ))),
                    box(width = 4, title = h2('Genes custom', icon('magnifying-glass')), collapsible = TRUE,
                      fluidRow(column(width=12,
                        checkboxGroupInput("gene_custom","Choose customization",choiceNames=list("1000 most variant genes","Center by gene"), choiceValues = list("mostvar","center"))
                      ))),
                    box(width = 4, title = h2('Online tools', icon('globe')), collapsible = TRUE,
                      fluidRow())
                        

                          )))))

)))
