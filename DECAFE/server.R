#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# packages <- c("shiny", "DT","shinydashboard","shinycssloaders","BiocManager", "ggplot2", "plotly", "reshape2", "factoextra", "FactoMineR", "devtools", "ggupset", "fgsea","DESeq2","ggpubr")
# new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
# if(length(new_packages)) install.packages(new_packages)
# if (!require("BiocManager", quietly = TRUE) && "BiocManager" %in% new_packages)
#   install.packages("BiocManager")
# BiocManager::install(new_packages,update=FALSE)
packages <- c("shiny", "DT", "shinydashboard", "shinycssloaders", "BiocManager", "ggplot2", "plotly", "reshape2", "factoextra", "FactoMineR", "devtools", "ggupset", 
"fgsea", "DESeq2", "ggpubr", "stringr", "ggrepel", "UpSetR", "ggdendro", "dendextend","gplots","svglite", "shinyBS","grid","gridExtra")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(new_packages, update = FALSE)
}

options(shiny.maxRequestSize=100000*1024^2)
library(shiny)
library(shinyBS)
library(ggplot2)
library(plotly)
library(reshape2)
library(ggpubr)
library(DESeq2)
library(factoextra)
library(FactoMineR)
library(ggupset)
library(fgsea)
library(stringr)
library(ggpubr)
library(ggrepel)
library(UpSetR)
library(ggdendro)
library(dendextend)
library(gplots)
library(gridExtra)
library(grid)



# Define server logic required to draw a histogram
function(input, output, session) {
# HOME Page

output$pca_Image <- renderImage({

    filename = "pca.png"

    list(src = filename,
      width = 800,
      alt = "overview RNA-Seq")

  }, deleteFile = FALSE)


output$anaDiff_Image <- renderImage({

    filename = "analysediff.png"
 
    list(src = filename,
      width = 1300,
      alt = "overview RNA-Seq")

  }, deleteFile = FALSE)

output$rna_Image <- renderImage({

    filename = "RNA_seq.png"

    list(src = filename,
      width = 1200,
         
         alt = "overview RNA-Seq")

  }, deleteFile = FALSE)


output$gsea_Image <- renderImage({

    filename = "gsea.png"
 
    list(src = filename,
      width = 1200,
      alt = "overview RNA-Seq")

  }, deleteFile = FALSE)



output$annot_Image <- renderImage({
    filename = "annot.png"
 
    list(src = filename,
      width = 1400,
      alt = "overview RNA-Seq")

  }, deleteFile = FALSE)


  # Overview Panel
  ## Reactive function

  # Open annotation file
  annotFile <-reactive({
    req(input$'annot-file')
      annot = read.delim(input$'annot-file'$datapath, row.names = 1)
      count = countFile()$count
    if(length(intersect(rownames(annot), colnames(count)))==0)
    showModal(modalDialog(
        title = "Invalid input",
        "The intersection between colnames of count matrix and rownames of annotation is empty!",
        easyClose = TRUE
      ))
    return(annot)
  })
  
  # Format annotation matrix
  annotProcess <- reactive({
    annot = annotFile()
    annot$condshiny = apply(annot, 1, function(x) paste(colnames(annot),"_", x, collapse = ','))
    if(input$sex){
      req(input$sexAnnot)
      sex= read.delim(input$sexAnnot$datapath,row.names=1)
      annot$sex = sex
    }

    return(annot)
  })

    
  # Open count file
  countFile <-reactive({
    req(input$file)

    notif <<- showNotification("Opening count matrix in progress", duration = 0)
    count = read.delim(input$file$datapath, sep='\t', row.names = 1, header=T,as.is=T)
    removeNotification(notif)
    
    if(!all(as.matrix(count) == as.integer(as.matrix(count)))){
      showModal(modalDialog(
        title = "Invalid input",
        "The  count matrix must not be normalized, only integer are accepted!",
        easyClose = TRUE
      ))
    }

    genefile = switch(input$org, 
    'hs' = 'humanGeneannot.rds',
    'mm' = 'mouseGeneannot.rds',
    )
    geneannot = readRDS(genefile)
    geneannot =  geneannot[rownames(count),]

    if (input$coding) {
      geneannot = geneannot[which(geneannot$biotype == 'protein_coding'), ]
      count = count[geneannot$GeneID, ]
    }

    if( input$sex) {
          observeEvent(input$sex, {
          showModal(modalDialog(
        title = "Upload Sex Information",
        "Upload a TSV file with 1 column representing the sex of your samples in the same order as your metadata TSV file.",
        footer = tagList(
          actionButton("dismiss_modal", "Dismiss")
        ),
        easyClose = FALSE 
      ))
    })

    observeEvent(input$dismiss_modal, {
      removeModal()  
    })

      X = switch(input$org, 
        'hs' = 'X',
        'mm' = 'mmuX',
      )
      Y = switch(input$org, 
        'hs' = 'Y',
        'mm' = 'mmuY',
      )
      geneannot = geneannot[-which(geneannot$seqname == X |geneannot$seqname == Y) , ]
      count = count[geneannot$GeneID  , ]
    
    }

    # if (input$sex) {
    #   count <- count[!grepl("^(Y|X|MT)$", geneannot$seqname), ]
    #   geneannot <- geneannot[!grepl("^(Y|X|MT)$", geneannot$seqname), ]
    # }
    
    return(list(count=count, geneannot=geneannot))
  })
 

   # Create named vector will all annotation (condition)
  conditionVector <- reactive({
    annot = annotProcess()$condshiny
    name = unique(annot)
    num  = c(1:length(unique(name)))

    choice_table = data.frame(name, num)
    choix = setNames(as.numeric(choice_table$num), choice_table$name)
  })


  
## Output 

  output$tableAnnot <- DT::renderDT(server = FALSE, {

    DT::datatable(
      annotFile(),
      extensions = c("Buttons"),
      options = list(
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = list(
          list(extend = "copy", text = "Copy", filename = "annot",
               exportOptions = list(
                 modifier = list(page = "current")
               )
          ),
          list(extend = "csv", text = "CSV", filename = "annot",
               exportOptions = list(
                 modifier = list(page = "all")
               )
          ),
          list(extend = "excel", text = "Excel", filename = "annot",
               exportOptions = list(
                 modifier = list(page = "all")
               )
          )
        )
      )
    )
  })


  output$counthead <- DT::renderDT(server = FALSE, {
    data = countFile()$count
    data = data[1:10, ]

    DT::datatable(
      data,
      extensions = c("Buttons"),
      options = list(
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = list(
          list(extend = "copy", text = "Copy", filename = "annot",
               exportOptions = list(
                 modifier = list(page = "current")
               )
          ),
          list(extend = "csv", text = "CSV", filename = "annot",
               exportOptions = list(
                 modifier = list(page = "all")
               )
          ),
          list(extend = "excel", text = "Excel", filename = "annot",
               exportOptions = list(
                 modifier = list(page = "all")
               )
          )
        )
      )
    )
  })


# upset plot  
dataUpset <- reactive({
  annot2 <- annotFile()
  lab <- apply(annot2, 2, unique, simplify = FALSE)
  tt <- unlist(lapply(seq_along(lab), function(x) {
    lapply(seq_along(lab[[x]]), function(y) paste0(names(lab)[x], '_', lab[[x]][y]))
  }))
  nb_cond <- length(tt)

  conditions <- as.data.frame(matrix(0, nrow = nrow(annot2), ncol = nb_cond)) # Création d'une dataframe remplie de 0
  colnames(conditions) <- tt
  # Assignation de 1 à chaque élément correspondant dans la dataframe
  for (i in 1:nrow(annot2)) {
    for (j in 1:ncol(annot2)) {
      conditions[i, paste0(names(annot2)[j], '_', annot2[i, j])] <- 1
    }
  }
  return(conditions)
})


output$upsetPlot <- renderPlot({
  df = dataUpset()
  df$numbers <- rowSums(df)
  upset(df, main.bar.color = '#262686',
            sets.bar.color = '#262686', 
            keep.order = TRUE,
            text.scale = 2,
            shade.color = "#CDCDE6", 
            sets.x.label = "Samples per Annotations", 
            mainbar.y.label = "Samples per intersections",
            nsets=20 ) 
})

output$downloadUpsetPlot <- downloadHandler(
      filename = function() {
            plot_title <- 'overview'
            format = input$format
            format = if(format=='png') 'pdf' else format
            paste(gsub(" ", "_", plot_title), "_", Sys.Date(), ".", format, sep = "")
      },
      contentType = paste0("image/",input$format),
      content = function(file) { 
        df = dataUpset()

        if(input$format == 'png')
          pdf(file,onefile = F,width = 15,height = 8)
        else
          svg(file,onefile = F,width = 15,height = 8) 


        upset_plot = upset(df, main.bar.color = '#262686',
                              sets.bar.color = '#262686', 
                              keep.order = TRUE,
                              text.scale = 2, 
                              shade.color = "#CDCDE6", 
                              sets.x.label = "Samples per Annotations", 
                              mainbar.y.label = "Samples per intersections")   
        print(upset_plot)   
        dev.off()    
      })
  


  output$cond1 <-renderUI({

    selectInput("cond1", label = "Choose group 1",
                  choices = conditionVector(), selected=1)
  })

  output$cond2 <-renderUI({

    selectInput("cond2", label = "Choose group 2",
                  choices = conditionVector(), selected=2)
  })
  
 
# create annot filter by input group1 and 2
 annotationName <- reactive({

    annot = annotProcess()
    annot_name_cond1 = unique(annot$condshiny)[as.numeric(input$cond1)]
    annot_name_cond2 = unique(annot$condshiny)[as.numeric(input$cond2)]
    annot_gp1 = rownames(annot)[which(annot$condshiny == annot_name_cond1)]
    annot_gp2 = rownames(annot)[which(annot$condshiny == annot_name_cond2)]

    res = list(
      annotName1=annot_name_cond1,
      annotName2= annot_name_cond2,
      annotGP1 = annot_gp1,
      annotGP2 = annot_gp2
    )
    return(res)

  })

# filter count and annot by annotationName
  intersectCond <-reactive({
    
    annot_name = annotationName()
    count = countFile()$count
    annot = annotProcess()

    annot_gp1 = annot_name$annotGP1
    annot_gp2 = annot_name$annotGP2
    
    all_annot = c(annot_gp1,annot_gp2)
    annot = annot[all_annot,]
    
    count_intersect = count[, intersect(rownames(annot), colnames(count))]
    annot_intersect= as.data.frame(annot[intersect(rownames(annot), colnames(count)), ])

    # annot_intersect$condshiny = apply(annot_intersect,1,function(x) paste(colnames(annot_intersect), "_", x, collapse=','))

    res = list(
      annot=annot_intersect,
      count=count_intersect
    )

    return(res)
  })


  DDS_cond <- reactive({
    intersect = intersectCond()

    count_intersect = intersect$count
    annot_intersect = intersect$annot

    zero_threshold = as.numeric(input$zero_threshold)
    countfilt = count_intersect[rowMeans(count_intersect == 0) <= (zero_threshold ), ]

    annot_intersect$condshiny = as.factor(annot_intersect$condshiny)
    A = as.character(unique(annot_intersect$condshiny)[1])
    B = as.character(unique(annot_intersect$condshiny)[2])
    annot_intersect$condshiny <- factor(annot_intersect$condshiny, levels = c(A, B))



    if(input$sex){
      #annot_intersect$sex<- as.factor(annot_intersect$sex)
      dds = DESeqDataSetFromMatrix(countData = data.matrix(countfilt),
                                    colData = annot_intersect[,c('condshiny','sex')],
                                    design = ~condshiny + sex)
    }
    else{
      dds = DESeqDataSetFromMatrix(countData = data.matrix(countfilt),
                                    colData = annot_intersect,
                                    design = ~condshiny )
    }

    dds = dds[rowSums(counts(dds)) >= 10]
    removeNotification(notif)
    return(dds)
  })

  vstNormalization_cond <- reactive({
    dds = DDS_cond()
    notif <<- showNotification("VST", duration = 0)
    normalized_counts =  assay(vst(dds))

    

    removeNotification(notif)
    return(normalized_counts)

  })

  geneFiltered <- reactive({
    intersect = intersectCond()

    count_intersect = intersect$count
    annot_intersect = intersect$annot

    zero_threshold = as.numeric(input$zero_threshold)
    countfilt = count_intersect[rowMeans(count_intersect == 0) <= (zero_threshold ), ]

    return(list(filtered=nrow(countfilt),total=nrow(count_intersect)))



    })

  output$nbGene <- renderUI({
    return(
      p(icon('circle-info'),paste0(" You have ",geneFiltered()$filtered , " filtered  and ",geneFiltered()$total," total genes"))
    )
    })

  vstAll <-reactive({

    count = countFile()$count
    annot = annotProcess()


    count_intersect = count[, intersect(rownames(annot), colnames(count))]

    annot_intersect= as.data.frame(annot[intersect(rownames(annot), colnames(count)), ])

    zero_threshold = as.numeric(input$zero_threshold)
    countfilt = count_intersect[rowMeans(count_intersect == 0) <= (zero_threshold ), ]


    if(input$sex){
      dds = DESeqDataSetFromMatrix(countData = data.matrix(countfilt),
                                    colData = annot_intersect,
                                    design = ~condshiny + sex)
    }
    else{
      dds = DESeqDataSetFromMatrix(countData = data.matrix(countfilt),
                                    colData = annot_intersect,
                                    design = ~condshiny )
    }

    dds = dds[rowSums(counts(dds)) >= 10]

     notif <<- showNotification("VST", duration = 0)
    normalized_counts =  assay(vst(dds))
    
      
    #normalized_counts = t(scale(t(normalized_counts), scale=FALSE)) # Normalization per gene 
    removeNotification(notif)
    return(list(normalized_counts=normalized_counts, annot_intersect=annot_intersect, count_intersect=count_intersect))


    })



  desqNormalization_cond <- reactive({
    dds = DDS_cond()
    dds2 = estimateSizeFactors(dds)
    normalized_counts = counts(dds2, normalized = TRUE)
    return(normalized_counts)

  })

countNormGenePlot <-reactive({
  
    normalized_counts =  desqNormalization_cond()
    annotName = annotationName()

    annot_name_cond1 = annotName$annotName1
    annot_name_cond2 = annotName$annotName2
    annot_gp1 = annotName$annotGP1
    annot_gp2 = annotName$annotGP2

    norm1 = normalized_counts[as.numeric(input$geneTarget),intersect(annot_gp1, colnames(normalized_counts))]
    norm2 = normalized_counts[as.numeric(input$geneTarget),intersect(annot_gp2, colnames(normalized_counts))]
    
    norm1_rm = as.vector(norm1)
    norm2_rm = as.vector(norm2)

    res = data.frame(
      count = c(norm1_rm, norm2_rm),
      condition = c(rep(annot_name_cond1,length(norm1_rm)), rep(annot_name_cond2,length(norm2_rm)))
    )
    return(res)
  })

# Heatmap 
output$nbGene2 <- renderUI({
    return(
      p(icon('circle-info'),paste0(" You have ",geneFiltered()$filtered , " filtered  and ",geneFiltered()$total," total genes"))
    )
    })


heatmapData <- reactive({

  if(input$data_heat =="all"){
    normalized_counts =  vstAll()$normalized_counts
    condition = vstAll()$annot_intersect$condshiny
    normalized_counts = t(scale(t(normalized_counts), scale=FALSE)) # Normalization per gene 
    gvar = apply(normalized_counts, 1, sd) 




  }
  else{
    normalized_counts = vstNormalization_cond()
    intersect = intersectCond()
    count_intersect = intersect$count
    annot_intersect= intersect$annot
    condition = annot_intersect$condshiny
    normalized_counts = t(scale(t(normalized_counts), scale=FALSE)) # Normalization per gene 
    gvar = apply(normalized_counts, 1, sd) 
    

  }
  return(list(gvar=gvar, normalized_counts = normalized_counts,condition=condition))



})

output$heatMap <- renderPlot({

  data_heatmap = heatmapData()
  gvar = data_heatmap$gvar
  normalized_counts = data_heatmap$normalized_counts
  condition = data_heatmap$condition
  mostvargenes = order(gvar, decreasing=TRUE)[1:as.numeric(input$nb_gene_heat)]
  
  if (length(unique(condition)) == 2){
    palette.annot =c("#CDCDE6", '#262686')}
  else {
    set.seed(Sys.Date())
    palette.annot = sapply(1:length(unique(condition)),function(x) paste0('#',paste0(sample(c(0:9,LETTERS[1:6]),6,T),collapse='')))}
  condition.colors = palette.annot
  
  

  gplots::heatmap.2(x=normalized_counts[mostvargenes,], 
    dendrogram="column",
    srtCol=45,
    col = "bluered",
    scale="none",
    trace="none",
    ColSideColors=condition.colors[ as.factor(condition)],
    labRow=FALSE,
    #ylab="Genes",
    xlab=NULL,margins = c(15, 3),
    keysize = 0.5,
    key.par = list(cex=1),
  )

  legend("left",
    legend=paste0(sapply(strsplit(unique(condition),','),paste,collapse = '\n'),'\n'),
    fill=condition.colors[ unique(as.factor(condition))], 
    cex=0.7
  )



})
output$downloadHeatmap <- downloadHandler(
      filename = function() {
            plot_title <- 'heatmap'
            format = input$format
            format = if(format=='png') 'pdf' else format
            paste(gsub(" ", "_", plot_title), "_", Sys.Date(), ".", format, sep = "")
      },
      contentType = "application/pdf",
      content = function(file) { 
          
data_heatmap = heatmapData()
  gvar = data_heatmap$gvar
  normalized_counts = data_heatmap$normalized_counts
  condition = data_heatmap$condition
  mostvargenes = order(gvar, decreasing=TRUE)[1:as.numeric(input$nb_gene_heat)]
  
  if (length(unique(condition)) == 2){
    palette.annot =c("#CDCDE6", '#262686')}
  else {
    set.seed(Sys.Date())
    palette.annot = sapply(1:length(unique(condition)),function(x) paste0('#',paste0(sample(c(0:9,LETTERS[1:6]),6,T),collapse='')))}
  condition.colors = palette.annot

  if(input$format == 'png')
    pdf(file,onefile = F,width = 10,height = 8)
  else
    svg(file,onefile = F,width = 10,height = 8) 
  

  gplots::heatmap.2(x=normalized_counts[mostvargenes,], 
    dendrogram="column",
    srtCol=45,
    col = "bluered",
    scale="none",
    trace="none",
    ColSideColors=condition.colors[ as.factor(condition)],
    labRow=FALSE,
    #ylab="Genes",
    xlab=NULL,margins = c(15, 3),
  )

  legend("left",
    legend=paste0(sapply(strsplit(unique(condition),','),paste,collapse = '\n'),'\n'),
    fill=condition.colors[ unique(as.factor(condition))], 
    cex=0.7
  )
        dev.off()   
      })



#PCA
  pcaVST <- reactive({

    intersect = intersectCond()
    count_intersect = intersect$count
    annot_intersect= intersect$annot

    vst = vstNormalization_cond()
    vst = t(scale(t(vst), scale=FALSE)) # Normalization per gene 

   
    
    
    gvar = apply(vst, 1, sd) 
    mostvargenes = order(gvar, decreasing=TRUE)[1:as.numeric(input$nb_gene)] # Keep most variable gene

    # Run pca 
    notif <<- showNotification("PCA in progress", duration = 0)

    res_pca = prcomp(t(vst[mostvargenes,]), scale. = TRUE)
    removeNotification(notif)

    pca_sample = get_pca_ind(res_pca)
    pca_gene =  get_pca_var(res_pca)

    coordinates = pca_sample$coord[,1:5] # Keep just 5th first PC

    pca_df = as.data.frame(
      cbind(
        coordinates, 
        rownames(annot_intersect), 
        annot_intersect$condshiny,
        colSums(count_intersect>0)
      )
    )
    colnames(pca_df)=c(paste0("Dim",1:5),"Sample","Group","NumberGene")

    eigen_value = get_eig(res_pca)[1:5,2]

    res = list(
      pca_df = pca_df, # For ggplot
      pca_sample = pca_sample, # For contribution
      pca_gene = pca_gene, # For contribution
      eig = round(eigen_value,2) # For percentage variance explained
    )
    return(res)

  })

  pcaGraph <- reactive({

    res.pca = pcaVST()
    pca_df = res.pca$pca_df
    eigen_value =res.pca$eig

    pca_sub = pca_df[,c(as.numeric(input$dim1), as.numeric(input$dim2), 6, 7, 8)]
    colnames(pca_sub)= c("x","y","sample","Group","NumberGene")
    pca_sub$x = as.numeric(pca_sub$x)
    pca_sub$y = as.numeric(pca_sub$y)
    pca_sub$Sample= paste(pca_sub$sample,'\nNumber gene : ',pca_sub$NumberGene)
    
    plot=ggplot(data=pca_sub, aes(x=x, y=y, col=Group, tooltip=Sample)) +
      geom_point(size=1) + theme_minimal() +
      labs(x = paste0("Dimension ", input$dim1," (",eigen_value[as.numeric(input$dim1)],"%)"),y = paste0("Dimension ", input$dim2," (",eigen_value[as.numeric(input$dim2)],"%)")) + 
      scale_color_manual(values=c("lightcoral", '#4ab3d6')) + 
      geom_vline(xintercept=0, linetype="dashed", color = "grey") +
      geom_hline(yintercept=0, linetype="dashed", color = "grey") +
      theme(legend.position="bottom")
    return(plot)
  })

  output$pcaVST <- renderPlotly({ 
    ggplotly(pcaGraph()) %>% 
      layout(legend = list(orientation = 'h', x = 0.45, y = 1.1)) %>% config(toImageButtonOptions=list(format=input$format, filename=paste0("PCA", "_", Sys.Date())))
  })


  
pca_alldownload <- reactive({
    
    normalized_counts =  vstAll()$normalized_counts
    annot_intersect = vstAll()$annot_intersect
    count_intersect = vstAll()$count_intersect
    normalized_counts = t(scale(t(normalized_counts), scale=FALSE)) # Normalization per gene 
    gvar = apply(normalized_counts, 1, sd) 
    mostvargenes = order(gvar, decreasing=TRUE)[1:as.numeric(input$nb_gene)]

    res_pca = prcomp(t(normalized_counts[mostvargenes,]), scale. = TRUE)
    pca_sample = get_pca_ind(res_pca)
    coordinates = pca_sample$coord[,1:5] # Keep just 5th first PC

    pca_df = as.data.frame(
      cbind(
        coordinates, 
        rownames(annot_intersect), 
        annot_intersect$condshiny,
        colSums(count_intersect>0)
      )
    )
    colnames(pca_df)=c(paste0("Dim",1:5),"Sample","Group","NumberGene")

    eigen_value = get_eig(res_pca)[1:5,2]

    res = list(
      pca_df = pca_df, # For ggplot
      eig = round(eigen_value,2) # For percentage variance explained
    )
    return(res)

  })

  output$downloadPCAPlot <- downloadHandler(
      filename = function() {
            paste(gsub(" ", "_", 'PCA_allConditions'), "_", Sys.Date(), ".",input$format, sep = "")
      },
      content = function(file) {
        showModal(modalDialog("Downloading PCA with all samples, wait", footer=NULL))
        on.exit(removeModal())
        res.pca = pca_alldownload()
        pca_df = res.pca$pca_df
        eigen_value =res.pca$eig

        pca_sub = pca_df[,c(as.numeric(input$dim1), as.numeric(input$dim2), 6, 7, 8)]
        colnames(pca_sub)= c("x","y","sample","Group","NumberGene")
        pca_sub$x = as.numeric(pca_sub$x)
        pca_sub$y = as.numeric(pca_sub$y)
        pca_sub$Sample= paste(pca_sub$sample,'\nNumber gene : ',pca_sub$NumberGene)
        
        plot=ggplot(data=pca_sub, aes(x=x, y=y, col=Group, tooltip=Sample)) +
          geom_point(size=1) + theme_minimal() +
          labs(x = paste0("Dimension ", input$dim1," (",eigen_value[as.numeric(input$dim1)],"%)"),
               y = paste0("Dimension ", input$dim2," (",eigen_value[as.numeric(input$dim2)],"%)")) + 
          geom_vline(xintercept=0, linetype="dashed", color = "grey") +
          geom_hline(yintercept=0, linetype="dashed", color = "grey")

        ggsave(file, plot, width = 8, height = 6, units = "in", dpi = 300, device = input$format,bg='white')
      })


  sampleTopDim1 <- reactive({
    
    res.pca = pcaVST()
    pca_sample = res.pca$pca_sample
    contrib = pca_sample$contrib[,as.numeric(input$dim1)]
    intersect = intersectCond()
    count_intersect = intersect$count


    tab = as.data.frame(
      cbind(
        name = colnames(count_intersect),
        contrib = contrib
      )
    )
    tab = tab[order(abs(as.numeric(tab$contrib)),decreasing=TRUE),]
    rownames(tab)= NULL

    # return(tab[1:min(10,nrow(tab)),]) # In the case where number sample least than 10
    return(tab)
  })

   output$sampleTop1 <- DT::renderDT(server = FALSE, {
      DT::datatable(
        sampleTopDim1(),
        extensions = c("Buttons"),
        options = list(
          dom = 'Bfrtip',
          scrollX = TRUE,
          buttons = list(
            list(extend = "copy", text = "Copy", filename = "sample1_contrib",
                 exportOptions = list(
                   modifier = list(page = "current")
                 )
            ),
            list(extend = "csv", text = "CSV", filename = "sample1_contrib",
                 exportOptions = list(
                   modifier = list(page = "all")
                 )
            ),
            list(extend = "excel", text = "Excel", filename = "sample1_contrib",
                 exportOptions = list(
                   modifier = list(page = "all")
                 )
            )
          )
        )
      )
    })

  sampleTopDim2 <- reactive({
   
    
    res.pca = pcaVST()
    pca_sample = res.pca$pca_sample
    contrib = pca_sample$contrib[,as.numeric(input$dim2)]
    intersect = intersectCond()
    count_intersect = intersect$count


    tab = as.data.frame(
      cbind(
        name = colnames(count_intersect),
        contrib = contrib
      )
    )
    tab = tab[order(abs(as.numeric(tab$contrib)),decreasing=TRUE),]
    rownames(tab)= NULL

    # return(tab[1:min(10,nrow(tab)),]) # In the case where number sample least than 10
    return(tab)
  })

   output$sampleTop2 <- DT::renderDT(server = FALSE, {
      DT::datatable(
        sampleTopDim2(),
        extensions = c("Buttons"),
        options = list(
          dom = 'Bfrtip',
          scrollX = TRUE,
          buttons = list(
            list(extend = "copy", text = "Copy", filename = "sample2_contrib",
                 exportOptions = list(
                   modifier = list(page = "current")
                 )
            ),
            list(extend = "csv", text = "CSV", filename = "sample2_contrib",
                 exportOptions = list(
                   modifier = list(page = "all")
                 )
            ),
            list(extend = "excel", text = "Excel", filename = "sample2_contrib",
                 exportOptions = list(
                   modifier = list(page = "all")
                 )
            )
          )
        )
      )
    })

   geneTopDim1 <- reactive({

   res.pca = pcaVST()
    pca_gene = res.pca$pca_gene
    contrib = pca_gene$contrib[,as.numeric(input$dim1)]
    intersect = intersectCond()
    count_intersect = intersect$count
    geneannot = countFile()$geneannot

    tab = as.data.frame(
      cbind(
        name = as.vector(geneannot[rownames(count_intersect),'GeneName']),
        contrib = contrib
      )
    )
    tab = tab[order(abs(as.numeric(tab$contrib)),decreasing=TRUE),]
    rownames(tab)= NULL

    # return(tab[1:min(10,nrow(tab)),]) # In the case where number sample least than 10
    return(tab)
  })

   output$geneTop1 <- DT::renderDT(server = FALSE, {
      DT::datatable(
        geneTopDim1(),
        extensions = c("Buttons"),
        options = list(
          dom = 'Bfrtip',
          scrollX = TRUE,
          buttons = list(
            list(extend = "copy", text = "Copy", filename = "gene1_contrib",
                 exportOptions = list(
                   modifier = list(page = "current")
                 )
            ),
            list(extend = "csv", text = "CSV", filename = "gene1_contrib",
                 exportOptions = list(
                   modifier = list(page = "all")
                 )
            ),
            list(extend = "excel", text = "Excel", filename = "gene1_contrib",
                 exportOptions = list(
                   modifier = list(page = "all")
                 )
            )
          )
        )
      )
    })
 geneTopDim2 <- reactive({

    res.pca = pcaVST()
    pca_gene = res.pca$pca_gene
    contrib = pca_gene$contrib[,as.numeric(input$dim2)]
    intersect = intersectCond()
    count_intersect = intersect$count
    geneannot = countFile()$geneannot


    tab = as.data.frame(
      cbind(
        name = as.vector(geneannot[rownames(count_intersect),'GeneName']),
        contrib = contrib
      )
    )
    tab = tab[order(abs(as.numeric(tab$contrib)),decreasing=TRUE),]
    rownames(tab)= NULL

    # return(tab[1:min(10,nrow(tab)),]) # In the case where number sample least than 10
    return(tab)


  })
 output$geneTop2 <- DT::renderDT(server = FALSE, {
    DT::datatable(
      geneTopDim2(),
      extensions = c("Buttons"),
      options = list(
        dom = 'Bfrtip',
        scrollX = TRUE,
        buttons = list(
          list(extend = "copy", text = "Copy", filename = "gene2_contrib",
               exportOptions = list(
                 modifier = list(page = "current")
               )
          ),
          list(extend = "csv", text = "CSV", filename = "gene2_contrib",
               exportOptions = list(
                 modifier = list(page = "all")
               )
          ),
          list(extend = "excel", text = "Excel", filename = "gene2_contrib",
               exportOptions = list(
                 modifier = list(page = "all")
               )
          )
        )
      )
    )
  })



# Dispersion analysis panel 

  resDeseq <- reactive({

    dds = DDS_cond()

    notif <<- showNotification("Differential analysis in progress", duration = 0)
    nb_thread = as.numeric(input$nb_thread)
    if(nb_thread > 1){
      BiocParallel::register(BiocParallel::MulticoreParam())
      dds = DESeq(dds,parallel = TRUE)
    }
    else{
      dds = DESeq(dds,parallel = FALSE)
    }
    geneannot = countFile()$geneannot
    table = results(dds)
    removeNotification(notif)
    table$name = as.vector(geneannot[rownames(table), 'GeneName'])
    table = as.data.frame(table)
    # table = table[order(abs(table$stat),decreasing=TRUE),]
    table = table[, c(7, 1:6)]
     

    return(list(deseq = dds, res = table) )

  })

  volcano <- reactive({

    table = resDeseq()$res
    tsPadj = as.numeric(input$ts_padj)
    tsFC = 1

    table$diffexpressed = "NO"
    table$diffexpressed[table$log2FoldChange > tsFC & table$padj < tsPadj] = "UP"
    table$diffexpressed[table$log2FoldChange < -tsFC & table$padj < tsPadj] = "DOWN"

    cols=c("lightcoral", "lightgrey","#4ab3d6")

    if(nrow(table[table$diffexpressed == 'NO', ])   == 0) { cols = c('lightcoral' , '#4ab3d6')}
    if(nrow(table[table$diffexpressed == 'DOWN',]) == 0) { cols = c('lightgrey', '#4ab3d6')}
    if(nrow(table[table$diffexpressed == 'UP',])   == 0)  { cols = c('lightgrey', 'lightcoral')}

    plot = ggplot(data=as.data.frame(table), aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, tooltip=name)) +
      geom_point(size=1) + theme_minimal()+
      geom_vline(xintercept=c(-tsFC, tsFC), col="#919191") +
      geom_hline(yintercept=-log10(tsPadj), col="#919191")+
      scale_color_manual(values=cols)
     
    return(plot)
  })

  output$volcano <- renderPlotly({ ggplotly(volcano())%>% config(toImageButtonOptions=list(format=input$format, filename=paste0("Volcano", "_", Sys.Date()))) })

  output$downloadVolcanoPlot <- downloadHandler(
      filename = function() {
            group1 <- annotationName()$annotName1
            group2 <- annotationName()$annotName2
            plot_title <- paste(group1, "vs", group2)
            paste(gsub(" ", "_", plot_title), "_", Sys.Date(), ".",input$format, sep = "")
      },
      content = function(file) {
        group1 <- annotationName()$annotName1
        group2 <- annotationName()$annotName2

        table = resDeseq()$res
        table = table[order(abs(table$stat),decreasing=TRUE),]
        tsPadj = as.numeric(input$ts_padj)
        tsFC = 1
        table$diffexpressed = "NO"
        table$diffexpressed[table$log2FoldChange > tsFC & table$padj < tsPadj] = "UP"
        table$diffexpressed[table$log2FoldChange < -tsFC & table$padj < tsPadj] = "DOWN"
        volcanoplot = volcano()  + ggtitle(paste(group1, "vs", group2)) +
        geom_text_repel(
                  data = head(as.data.frame(table[which(table$diffexpressed != "NO"),]), 15),
                  aes(label = name),show.legend=FALSE,
                  box.padding = 0.5, point.padding = 0.1,
                  segment.color = 'grey', segment.size = 0.2,
                  nudge_y = 0.2
                )
        print(file)
        ggsave(file, volcanoplot,  width = 8, height = 6, units = "in", dpi = 300, device = input$format, bg='white')
      })


  output$degTable <- DT::renderDT(server = FALSE, {
    table = resDeseq()$res
    table = table[order(abs(table$stat),decreasing=TRUE),]
    DT::datatable(
      table,
      extensions = c("Buttons"),
      options = list(
        dom = 'Bfrtip', 
        scrollX = TRUE,
        buttons = list(
          list(extend = "copy", text = "Copy", filename = "res_deseq",
               exportOptions = list(
                 modifier = list(page = "current")
               )
          ),
          list(extend = "csv", text = "CSV", filename = "res_deseq",
               exportOptions = list(
                 modifier = list(page = "all")
               )
          ),
          list(extend = "excel", text = "Excel", filename = "res_deseq",
               exportOptions = list(
                 modifier = list(page = "all")
               )
          )
        )
      )
    )
  })

  # Boxplot Panel
  geneTargetChoice <- reactive({
    # count= countFile()$count
    # geneannot = countFile()$geneannot
    # genename_list = as.vector(geneannot[rownames(count),'GeneName'])
    genename_list = as.vector(resDeseq()$res$name)

    name = genename_list
    num  = c(1:length(name))

    choiceTable = data.frame(name, num)
    choix = setNames(as.numeric(choiceTable$num), choiceTable$name)

    return(list(choix_indice=choix, choix_name=choiceTable$name[choix]))

  })

  output$geneTarget <-renderUI({  
    return(selectInput("geneTarget", label = "Choose gene target to observe",
                  choices = geneTargetChoice()$choix_indice))
  })

# boxplot
 output$bpGeneTarget <- renderPlot({
    res = countNormGenePlot()
    res$count = as.numeric(res$count)
    
    res$genetarget = geneTargetChoice()$choix_name[as.numeric(input$geneTarget)]
    padj = resDeseq()$res$padj[which(resDeseq()$res$name == res$genetarget)]
   
    df_p_val <- data.frame(
      group1 = unique(res$condition)[1],
      group2 = unique(res$condition)[2],
      label = padj,
      y.position = max(res$count)*1.1
    )

    ggplot(res, aes(x=condition,y=count, color=condition)) + geom_boxplot() +
      rotate_x_text(45)+ labs(x = "dispersion", y = paste0('count normDESq2 of target gene: ',res$genetarget))+
      scale_color_manual(values=c("lightcoral", '#4ab3d6')) + 
      scale_x_discrete(labels=NULL) +  theme_minimal() + theme(legend.position="right", legend.text=element_text(size=10)) + 
      stat_pvalue_manual(df_p_val, xmin = "group1", xmax = "group2", label = "label", y.position = "y.position") + 
      stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="#262686", fill="#262686")
  })

  output$histGeneTarget <- renderPlot({
    res = countNormGenePlot()
    res$genetarget = geneTargetChoice()$choix_name[as.numeric(input$geneTarget)]

    ggplot(res, aes(x=count, color=condition)) +
      geom_histogram(fill="white", alpha=0.5, position="identity") + 
      theme_minimal() + theme(legend.position="top", legend.text=element_text(size=10)) +  
      labs(x = paste0('count normDESq2 of target gene: ',res$genetarget),y = "number of samples") + 
      scale_color_manual(values=c("lightcoral", '#4ab3d6'))
  })

  
output$downloadboxplot <- downloadHandler(
      filename = function() {
            plot_title <- 'boxplot'
            format = input$format
            format = if(format=='png') 'pdf' else format
            paste(gsub(" ", "_", plot_title), "_", Sys.Date(), ".", format, sep = "")
      },
      contentType = "application/pdf",
      content = function(file) { 
          res = countNormGenePlot()
          res$genetarget = geneTargetChoice()$choix_name[as.numeric(input$geneTarget)]
          res$count = as.numeric(res$count)
          padj = resDeseq()$res$padj[which(resDeseq()$res$name == res$genetarget)]

          df_p_val <- data.frame(
          group1 = unique(res$condition)[1],
          group2 = unique(res$condition)[2],
          label = padj,
          y.position = max(res$count)*1.1)


          boxplot = ggplot(res, aes(x=condition,y=count, color=condition)) + geom_boxplot() +
                    rotate_x_text(45)+ labs(x = "condition", y = paste0('count normDESq2 of target gene: ',res$genetarget)) +
                    scale_color_manual(values=c("lightcoral", '#4ab3d6')) + 
                    scale_x_discrete(labels=NULL) +  theme_minimal() + theme(legend.position="right", legend.text=element_text(size=10)) + 
                    stat_pvalue_manual(df_p_val, xmin = "group1", xmax = "group2", label = "label", y.position = "y.position") + 
                    stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="#262686", fill="#262686")
          
          histo = ggplot(res, aes(x=count, color=condition)) +
                  geom_histogram(fill="white", alpha=0.5, position="identity") + 
                  theme_minimal() + theme(legend.position="top", legend.text=element_text(size=10)) +  
                  labs(x = paste0('count normDESq2 of target gene: ',res$genetarget),y = "number of samples") + 
                  scale_color_manual(values=c("lightcoral", '#4ab3d6'))

        if(input$format == 'png'){
            pdf(file, width = 10, height = 8)
              print(histo)
              print(boxplot)
            dev.off()}   
        else{
          svg(file, width = 10, height = 8)
            print(ggarrange(histo, boxplot,
          labels = c("A", "B"),
          ncol = 2, nrow = 1))
            dev.off()}

      })



  # select collection pathway
  appendCollection <- reactive({
    pathwayfile = switch(input$org, 
    'hs' = 'human_pathays.rds',
    'mm' = 'mouse_pathays.rds',
    )
    pathways_list = readRDS(pathwayfile)
    return(list(namesp=names(pathways_list), pathways=pathways_list))
  })

  output$dbpath <- renderUI({
    selectInput('path_list', label = 'Select list of pathways', choices= appendCollection()$namesp, multiple = TRUE)
  })

  # photo pathways
  output$path_Image <- renderImage({
    imagepath = switch(input$org, 
    'hs' = 'human_pathway.png',
    'mm' = 'mouse_pathway.png',
    )
    list(src = imagepath,
      width = 800,
         alt = "image_pathways")
  }, deleteFile = FALSE)


  gsea <- eventReactive(input$gogsea,{

    pathways_list = appendCollection()$pathways
    pathwaySelect = input$path_list
    pathways_list = pathways_list[pathwaySelect]
    
    pathways = pathways_list[[1]]
    names(pathways) = paste0(pathwaySelect[1], '_:_', names(pathways))
    if (length(pathways_list) > 1){
    for (i in 2:length(pathways_list)){
      path_temp = pathways_list[[i]]
      names(path_temp) = paste0(pathwaySelect[i], '_:_', names(path_temp))
      pathways = append(pathways, path_temp)
    }}
    
    table = resDeseq()$res
    vec = table$stat
    names(vec) = table$name

    notif <<- showNotification("GSEA in progress", duration = 0)
    nb_thread = as.numeric(input$nb_thread)
    
    if(nb_thread > 1){
      res = fgseaSimple(pathways,vec,BPPARAM = BiocParallel::MulticoreParam(nb_thread), nperm = 1000)
    }
  
    else{
      res = fgseaSimple(pathways,vec, nperm = 1000)
    }
    res.original = res
    splitcolpath = strsplit(res$pathway, '_:_')
    res$collection = sapply(splitcolpath, function(x) x[1])
    res$pathway = sapply(splitcolpath, function(x) x[2])
    res$pathway = str_replace_all(res$pathway, '_', ' ')
    removeNotification(notif)

    lE_list = res$leadingEdge
    lE_vector = sapply(lE_list, paste, collapse=", ")
    res$leadingEdge = lE_vector
    res = as.data.frame(na.omit(res))
    res_sig = res[which(as.numeric(res$padj) < 0.05),]
    res_sort = res_sig[order(abs(as.numeric(res_sig$NES)),decreasing=TRUE),]
    res_sort = res_sort[, c(9, 1:5, 7:8)]

    return(list(sort=res_sort,orig = res.original, vec = vec, pathways=pathways ))
  })



  



  # GSEA Panel
  output$gsea <- DT::renderDT(server = FALSE, {
      DT::datatable(
        gsea()$sort,
        extensions = c("Buttons"),
        caption= "Pathways enrichment",
        selection = 'single',
        options = list(
          dom = 'Bfrtip',
          scrollX = TRUE,
          buttons = list(
            list(extend = "copy", text = "Copy", filename = "gsea",
                 exportOptions = list(
                   modifier = list(page = "current")
                 )
            ),
            list(extend = "csv", text = "CSV", filename = "gsea",
                 exportOptions = list(
                   modifier = list(page = "all")
                 )
            ),
            list(extend = "excel", text = "Excel", filename = "gsea",
                 exportOptions = list(
                   modifier = list(page = "all")
                 )
            )
          )
        )
      )
    })

  browse <- eventReactive(input$browsebutton, {
    species = switch(input$org, 
    'hs' = 'human',
    'mm' = 'mouse',
    )
    indice_clickpath = input$gsea_rows_selected
    if (length(indice_clickpath) == 1){
    clickpath = gsea()$sort[indice_clickpath, 'pathway']
    clickpath = str_replace_all(clickpath, ' ', '_')
    if(gsea()$sort[indice_clickpath, 'collection']=='sigGeNeHetX'){
      load('signatures.rda')
      annot = signatures$annotation
      pmid = strsplit(strsplit(annot[clickpath,'src'],';')[[1]][2],'[.]')[[1]][2]
      url = paste0('https://pubmed.ncbi.nlm.nih.gov/',pmid,'/')
    }else{
      url = paste0("https://www.gsea-msigdb.org/gsea/msigdb/", species, "/geneset/", clickpath ,".html")
    }
    browseURL(url)
    return(clickpath)
}})

  output$comment <- renderText({
    clickpath = browse()
    clickpath = str_replace_all(clickpath, '_', ' ')
    return(paste0('Opening url for : ', clickpath))
  })



  graph <- reactive({
    res <- gsea()$sort
    rownames(res) <- make.unique(res$pathway)
    
    res_filtered <- res[res$pval <= 0.01, ]
      if (nrow(res_filtered) < 100) {
        name_pathway <- rownames(res_filtered)
    } else {
        name_pathway <- rownames(res_filtered)[1:100]
    }   
    collection <- res$collection
    mat <- matrix(0, nrow = length(name_pathway), ncol = length(name_pathway))
    dimnames(mat) <- list(name_pathway, name_pathway)
    
    for (i in name_pathway) {
      for (j in name_pathway) {
        le_i <- strsplit(res[i, 'leadingEdge'], ',')[[1]]
        le_j <- strsplit(res[j, 'leadingEdge'], ',')[[1]]
        mat[i, j] <- 1 - (length(intersect(le_i, le_j)) / length(union(le_i, le_j)))
      }
    }
    
  hc <- hclust(as.dist(mat), method = "ward")
  return(list(hc = hc, collection = collection))

    })

  output$treePlot <- renderPlot({
      dend <- as.dendrogram(graph()$hc) 
      par(cex=0.6, mar=c(6, 1, 1, 60))
      dend %>% set("branches_k_color", k = min(as.numeric(input$k), length(graph()$hc$labels))) %>% plot(horiz = TRUE)

    }, width = 1200, height = 1300, res = 96)

  dotplotGSEA<- reactive({
    df <- gsea()$sort
    colnames(df)[which(colnames(df)=='size')]<-"Count"
    colnames(df)[which(colnames(df)=='p.adj')]<-"padj"

    df$sign = sign(df$NES)
    df$sign[df$sign == -1] <- "suppressed"
    df$sign[df$sign == 1]<- "activated"

    df$pathway = str_replace_all(df$pathway,'_',' ')
    df$pathway=str_replace_all(df$pathway, "(.{25})([^\\s])", "\\1\\2-\n")

    activated = df[which(df$NES > 0),]
    activated = activated[1:min(nrow(activated),as.numeric(input$nbDotplot)),]
    suppressed = df[which(df$NES < 0),]
    suppressed = suppressed[1:min(nrow(suppressed),as.numeric(input$nbDotplot)),]

    sub = as.data.frame(rbind(activated,suppressed))
    sub =  sub[order(abs(sub$NES)),]


    plot=ggplot(sub, aes(x=NES, y=pathway, colour=padj, size=Count)) +
    geom_point() + facet_grid(.~sign,scales="free_x")+scale_color_gradient(low="blue", high="red")
    return(plot)
    
  })

  output$dotplot<-renderPlotly({ ggplotly(dotplotGSEA())%>% config(toImageButtonOptions=list(format=input$format, filename=paste0("Dotplot", "_", Sys.Date()))) })



  qGseaTable <- function (
    pathways, stats, fgseaRes, gseaParam = 1,
    colwidths = c(4,3, 0.8, 1.2, 1.2),
    rename=NULL, colors= list(low = '#ca0020', mid="#f7f7f7", high = '#0571b0')
  ) {

  require(fgsea)
  rnk = rank(-stats)
  ord = order(rnk)
  statsAdj = stats[ord]
  statsAdj = sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj = statsAdj / max(abs(statsAdj))
  n = length(statsAdj)
  pathways = lapply(pathways, function(p) {
      unname(as.vector(na.omit(match(p, names(statsAdj)))))
  })

  new_names = str_replace_all(names(pathways),'_',' ')
  new_names=str_replace_all(new_names, "(.{25})([^\\s])", "\\1\\2-\n")
  new_names = sub(".*?:", "", new_names)


  rename=new_names
  names(rename)=names(pathways)



  ps = lapply(names(pathways), function(pn) {
      p = pathways[[pn]]
      annotation = fgseaRes[match(pn, fgseaRes$pathway), ]

      list(
        textGrob(rename[pn], just = "right", x = unit(0.95, "npc"), gp=gpar(fontsize=12)),
          ggplot() + 
            geom_segment(
              aes(
                x = p, 
                xend = p, 
                y = 0,
                        yend = statsAdj[p],
                        col=p
                    ), 
                    size = 0.2
                ) + 
                scale_x_continuous(
                  limits = c(0, length(statsAdj)), 
                  expand = c(0, 0)
                ) + 
                scale_y_continuous(
                  limits = c(-1, 1),
                  expand = c(0, 0)
                ) + 
                xlab(NULL) + 
                ylab(NULL) +
              scale_colour_gradient2(
                midpoint=n/2,
                low=colors$low,
                mid=colors$mid,
                high=colors$high
              ) +
              theme(
                panel.background = element_rect(fill = "transparent",
                                  colour = NA_character_),
                legend.position="none", 
                axis.line = element_blank(),
                  axis.text = element_blank(), 
                  axis.ticks = element_blank(),
                  panel.grid = element_blank(), 
                  axis.title = element_blank(),
                  plot.margin = rep(unit(0, "null"), 4), 
                  panel.spacing = rep(unit(0,"null"), 4)
                ), 
            textGrob(sprintf("%.2f", annotation$NES)),
          textGrob(sprintf("%.1e", annotation$pval)), 
          textGrob(sprintf("%.1e",annotation$padj))
        )
    })
    
    rankPlot = ggplot() + 
      geom_blank() + 
      scale_x_continuous(
        limits = c(0,length(statsAdj)), 
        expand = c(0, 0)
      ) + 
      scale_y_continuous(
        limits = c(-1,1), 
        expand = c(0, 0)
      ) + 
      xlab(NULL) + 
      ylab(NULL) + 
      theme(
        panel.background = element_rect(fill = "transparent",
                                  colour = NA_character_),
          axis.line = element_blank(), 
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), 
            panel.grid = element_blank(),
            axis.title = element_blank(), 
            plot.margin = unit(c(0,0, 0.5, 0), "npc"), 
            panel.spacing = unit(c(0, 0, 0, 0), "npc")
        )
    grobs = c(
      lapply(
        c("Pathway", "Gene ranks", "NES", "pval","padj"), 
        textGrob
      ), 
      unlist(ps, recursive = FALSE), 
      list(nullGrob(),rankPlot, nullGrob(), nullGrob(), nullGrob())
    )
    grobsToDraw = rep(colwidths != 0, length(grobs)/length(colwidths))
    gridExtra::grid.arrange(
      grobs = grobs[grobsToDraw], 
      ncol = sum(colwidths !=0),
      widths = colwidths[colwidths != 0]
    )


} 



  gseaTablePlot<-reactive({
    gsea = gsea()
    xgsea = gsea$orig
    xgsea=xgsea[which(as.numeric(xgsea$padj) < 0.05),]
    xgsea_sub=head(xgsea[order(abs(xgsea$NES),decreasing=TRUE),],input$nbGseaPlot)
    plot=qGseaTable(gsea$pathways[xgsea_sub$pathway],gsea$vec,xgsea)
    return(plot)


    })

  output$gseaPlot<-renderPlot({ gseaTablePlot() })



  }
