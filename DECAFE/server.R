#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
packages <- c("shiny", "DT","shinydashboard","shinycssloaders","BiocManager", "ggplot2", "plotly", "reshape2", "factoextra", "FactoMineR", "devtools", "ggupset", "fgsea","DESeq2","ggpubr")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
if (!require("BiocManager", quietly = TRUE) && "BiocManager" %in% new_packages)
  install.packages("BiocManager")
BiocManager::install(new_packages,update=FALSE)

options(shiny.maxRequestSize=100000*1024^2)
library(shiny)
library(ggplot2)
library(plotly)
library(reshape2)
library(ggpubr)
library(DESeq2)
library(factoextra)
library(FactoMineR)
library(ggupset)
library(fgsea)

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
    # When input$n is 3, filename is ./images/image3.jpeg
    filename = "annot.png"
 
    # Return a list containing the filename and alt text
    list(src = filename,
      width = 1500,
      alt = "overview RNA-Seq")

  }, deleteFile = FALSE)


  # Overview Panel
  ## Reactive function

  # Open annotation file
  annotFile <-reactive({
    req(input$'annot-file')

    annot = read.delim(input$'annot-file'$datapath, row.names = 1)

  })
  
  # Open count file
  countFile <-reactive({
    req(input$file)

    notif <<- showNotification("Opening count matrix in progress", duration = 0)
    count = read.delim(input$file$datapath, sep='\t', row.names = 1, header=T,as.is=T)
    removeNotification(notif)
    
    return(count)
  })

  # Format annotation matrix
  annotProcess <- reactive({

    annot = annotFile()
                                 
    annot$name = apply(annot, 1, function(x) paste(colnames(annot),"_", x, collapse = ','))
    annot_combn = annot %>%
      mutate(
          all_list = as.list(strsplit(name, ","))
      )
     
    return(annot_combn)
  })


  # Create named vector will all annotation (condition)
  conditionVector <- reactive({
    annot = annotProcess()$name
    name = unique(annot)
    num  = c(1:length(unique(name)))

    choice_table = data.frame(name, num)
    choix = setNames(as.numeric(choice_table$num), choice_table$name)

  
  })


  
  ## Output 
  output$upsetPlot <- renderPlot(

    ggplot(data = annotProcess(),mapping = aes(x = all_list)) +
      geom_bar(fill = '#262686') +
      scale_x_upset(order_by = "degree")+
      geom_text(stat='count', aes(label=..count..), vjust=2, size=5, color="white") +
      labs(x='Annotation combination', y='Distribution') +
      theme(axis.title = element_text(size=15))+
      theme_combmatrix(combmatrix.label.text = element_text(color ='#262686', size=12)) + 
      axis_combmatrix(levels = levels(as.factor(unlist(annotProcess()$all_list)))) 
    
  )
  
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
                 modifier = list(page = "all")
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

  output$cond1 <-renderUI({

    selectInput("cond1", label = "Choose group 1",
                  choices = conditionVector(), selected=1)
  })

  output$cond2 <-renderUI({

    selectInput("cond2", label = "Choose group 2",
                  choices = conditionVector(), selected=2)
  })

# Boxplot Panel


  geneTargetChoice <- reactive({
    count= countFile()
    gene = rownames(count)

    name = gene
    num  = c(1:length(name))

    choiceTable = data.frame(name, num)
    choix = setNames(as.numeric(choiceTable$num), choiceTable$name)

    return(choix)

  })

  output$geneTarget <-renderUI({
    
    
   
    return(selectInput("geneTarget", label = "Choose gene target to observe",
                  choices = geneTargetChoice()))

  })



  annotationName <- reactive({

    annot = annotProcess()


    annot_name_cond1 = unique(annot$name)[as.numeric(input$cond1)]
    annot_name_cond2 = unique(annot$name)[as.numeric(input$cond2)]
    annot_gp1 = rownames(annot)[which(annot$name == annot_name_cond1)]
    annot_gp2 = rownames(annot)[which(annot$name == annot_name_cond2)]

    res = list(
      annotName1=annot_name_cond1,
      annotName2= annot_name_cond2,
      annotGP1 = annot_gp1,
      annotGP2 = annot_gp2
    )

    return(res)

  })

  intersectCond <-reactive({
    
    annot_name = annotationName()
    count = countFile()
    annot = annotProcess()

    annot_gp1 = annot_name$annotGP1
    annot_gp2 = annot_name$annotGP2

    
    all_annot = c(annot_gp1,annot_gp2)
    annot = annot[all_annot,]

  
    
    count_intersect = count[, intersect(rownames(annot), colnames(count))]
    annot_intersect= as.data.frame(annot[intersect(rownames(annot), colnames(count)), ])

    annot_intersect$condshiny = apply(annot_intersect,1,function(x) paste(colnames(annot_intersect), "_", x, collapse=','))

  
    res = list(
      annot=annot_intersect,
      count=count_intersect
    )

    return(res)
  })



  DDS <- reactive({

    intersect = intersectCond()
    notif <<- showNotification("DESeq", duration = 0)
    count_intersect = intersect$count
    annot_intersect = intersect$annot
    
      dds = DESeqDataSetFromMatrix(
        countData = count_intersect,
        colData = annot_intersect,
        design = ~condshiny
      )
    removeNotification(notif)
      return(dds)
  })


  vstNormalization <- reactive({

    dds = DDS()
    notif <<- showNotification("VST", duration = 0)
    normalized_counts =  assay(vst(dds))

    removeNotification(notif)
    return(normalized_counts)

  })

countNormGenePlot <-reactive({
  
    normalized_counts =  vstNormalization()
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

  output$bpGeneTarget <- renderPlot({
    
    res = countNormGenePlot()
    
    ggboxplot(res, 'condition','count',col='condition',shape='condition',xlab=NULL)+
      rotate_x_text(45)+
      scale_x_discrete(labels=NULL)+
      stat_compare_means(label = "p.signif")
  })

  output$histGeneTarget <- renderPlot({
    res = countNormGenePlot()
    
    ggplot(res, aes(x=count, color=condition)) +
      geom_histogram(fill="white", alpha=0.5, position="identity") + 
      theme(legend.position="top")
  })

  

#PCA
  pcaVSD <- reactive({

    intersect = intersectCond()
    count_intersect = intersect$count
    annot_intersect= intersect$annot

    vst = vstNormalization()
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
        annot_intersect$name,
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

    res.pca = pcaVSD()
    pca_df = res.pca$pca_df
    print(dim(pca_df))
    eigen_value =res.pca$eig

    pca_sub = pca_df[,c(as.numeric(input$dim1), as.numeric(input$dim2), 6, 7,8)]
    colnames(pca_sub)= c("x","y","sample","Group","NumberGene")
    pca_sub$x = as.numeric(pca_sub$x)
    pca_sub$y = as.numeric(pca_sub$y)
    pca_sub$Sample= paste(pca_sub$sample,'\nNumber gene : ',pca_sub$NumberGene)
    
    plot=ggplot(data=pca_sub, aes(x=x, y=y, col=Group, tooltip=Sample)) +
      geom_point(size=1) + theme_minimal() +
      labs(x = paste0("Dimension ", input$dim1," (",eigen_value[as.numeric(input$dim1)],"%)"),y = paste0("Dimension ", input$dim2," (",eigen_value[as.numeric(input$dim2)],"%)")) + 
      geom_vline(xintercept=0, linetype="dashed", color = "grey") +
      geom_hline(yintercept=0, linetype="dashed", color = "grey") +
      theme(legend.position="bottom")
    return(plot)
  })

  output$pcaVST <- renderPlotly({ 
    ggplotly(pcaGraph()) %>% 
      layout(legend = list(orientation = 'h', x = 0.45, y = 1.1)) 
  })


  sampleTopDim1 <- reactive({
    
    res.pca = pcaVSD()
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

    return(tab[1:min(10,nrow(tab)),]) # In the case where number sample least than 10
  })

   output$sampleTop1 <- DT::renderDT(server = FALSE, {
      DT::datatable(
        sampleTopDim1(),
        extensions = c("Buttons"),
        options = list(
          dom = 'Bfrtip',
          buttons = list(
            list(extend = "copy", text = "Copy", filename = "sample1_contrib",
                 exportOptions = list(
                   modifier = list(page = "all")
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
   
    
    res.pca = pcaVSD()
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

    return(tab[1:min(10,nrow(tab)),]) # In the case where number sample least than 10

  })

   output$sampleTop2 <- DT::renderDT(server = FALSE, {
      DT::datatable(
        sampleTopDim2(),
        extensions = c("Buttons"),
        options = list(
          dom = 'Bfrtip',
          buttons = list(
            list(extend = "copy", text = "Copy", filename = "sample2_contrib",
                 exportOptions = list(
                   modifier = list(page = "all")
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

   res.pca = pcaVSD()
    pca_gene = res.pca$pca_gene
    contrib = pca_gene$contrib[,as.numeric(input$dim1)]
    intersect = intersectCond()
    count_intersect = intersect$count


    tab = as.data.frame(
      cbind(
        name = rownames(count_intersect),
        contrib = contrib
      )
    )
    tab = tab[order(abs(as.numeric(tab$contrib)),decreasing=TRUE),]
    rownames(tab)= NULL

    return(tab[1:min(10,nrow(tab)),]) # In the case where number sample least than 10

  })

   output$geneTop1 <- DT::renderDT(server = FALSE, {
      DT::datatable(
        geneTopDim1(),
        extensions = c("Buttons"),
        options = list(
          dom = 'Bfrtip',
          buttons = list(
            list(extend = "copy", text = "Copy", filename = "gene1_contrib",
                 exportOptions = list(
                   modifier = list(page = "all")
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

    res.pca = pcaVSD()
    pca_gene = res.pca$pca_gene
    contrib = pca_gene$contrib[,as.numeric(input$dim2)]
    intersect = intersectCond()
    count_intersect = intersect$count


    tab = as.data.frame(
      cbind(
        name = rownames(count_intersect),
        contrib = contrib
      )
    )
    tab = tab[order(abs(as.numeric(tab$contrib)),decreasing=TRUE),]
    rownames(tab)= NULL

    return(tab[1:min(10,nrow(tab)),]) # In the case where number sample least than 10


  })
 output$geneTop2 <- DT::renderDT(server = FALSE, {
    DT::datatable(
      geneTopDim2(),
      extensions = c("Buttons"),
      options = list(
        dom = 'Bfrtip',
        buttons = list(
          list(extend = "copy", text = "Copy", filename = "gene2_contrib",
               exportOptions = list(
                 modifier = list(page = "all")
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

    intersect = intersectCond()

    count_intersect = intersect$count
    annot_intersect = intersect$annot

    zero_threshold = as.numeric(input$zero_threshold)
    countfilt = count_intersect[rowMeans(count_intersect == 0) <= (zero_threshold ), ]

    dds = DESeqDataSetFromMatrix(countData = countfilt,
                                    colData = annot_intersect,
                                    design = ~name)

    dds = dds[rowSums(counts(dds)) >= 10]

    notif <<- showNotification("Differential analysis in progress", duration = 0)
    nb_thread = as.numeric(input$nb_thread)
    if(nb_thread > 1){
      BiocParallel::register(BiocParallel::MulticoreParam())
      dds = DESeq(dds,parallel = TRUE)
    }
    else{
      dds = DESeq(dds,parallel = FALSE)
    }

    table = results(dds)
    removeNotification(notif)
    table$name = rownames(table)
    table = as.data.frame(table)
    table = table[order(abs(table$stat),decreasing=TRUE),]

    return(list(deseq = dds, res = table) )

  })

  volcano <- reactive({

    table = resDeseq()$res
    tsPadj = as.numeric(input$ts_padj)
    tsFC = 2

          
    table$diffexpressed = "NO"
    table$diffexpressed[table$log2FoldChange > tsFC & table$padj < tsPadj] = "UP"
    table$diffexpressed[table$log2FoldChange < -tsFC & table$padj < tsPadj] = "DOWN"

    cols=c("gold", "black","blue")

    if(nrow(subset(table,diffexpressed = 'NO'))   == 0) { cols = c('gold' , 'blue')}
    if(nrow(subset(table,diffexpressed = 'DOWN')) == 0) { cols = c('black', 'blue')}
    if(nrow(subset(table,diffexpressed = 'UP'))   ==0)  { cols = c('black', 'gold')}

    plot = ggplot(data=as.data.frame(table), aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, tooltip=name)) +
      geom_point(size=1) + theme_minimal()+
      geom_vline(xintercept=c(-tsFC, tsFC), col="red") +
      geom_hline(yintercept=-log10(tsPadj), col="red")+
      scale_color_manual(values=cols)
     
    return(plot)
  })

  output$volcano <- renderPlotly({ ggplotly(volcano()) })

  output$degTable <- DT::renderDT(server = FALSE, {
    DT::datatable(
      resDeseq()$res,
      extensions = c("Buttons"),
      options = list(
        dom = 'Bfrtip',
        buttons = list(
          list(extend = "copy", text = "Copy", filename = "res_deseq",
               exportOptions = list(
                 modifier = list(page = "all")
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

  gsea <- reactive({

    pathways = readRDS('pathway_human.rds')
    table = resDeseq()$res
    vec = table$stat
    names(vec) = table$name

    notif <<- showNotification("GSEA in progress", duration = 0)
    nb_thread = as.numeric(input$nb_thread)
    
    if(nb_thread > 1){
      res = fgsea(pathways,vec,BPPARAM = BiocParallel::MulticoreParam(nb_thread))
    }
    else{
      res = fgsea(pathways,vec)
    }
    
    removeNotification(notif)

    lE_list = res$leadingEdge

    
    lE_vector = sapply(lE_list, paste, collapse=",")
    res$leadingEdge = lE_vector
    res = as.data.frame(na.omit(res))
    res_sig = res[which(as.numeric(res$padj) < 0.05),]
    res_sort = res_sig[order(abs(as.numeric(res_sig$NES)),decreasing=TRUE),]
    return(res_sort)
  })

  output$barplotGsea <- renderPlot({

    df = gsea()[1:15,]
    df$NES = as.numeric(df$NES)
# En dehors des barres
  ggplot(data=df, aes(x=NES, y=pathway, fill = padj)) +
  geom_bar(stat="identity")+
  theme_minimal()


    })

  # GSEA Panel
  output$gsea <- DT::renderDT(server = FALSE, {
      DT::datatable(
        gsea(),
        extensions = c("Buttons"),
        caption= "Pathways enrichment",
        options = list(
          dom = 'Bfrtip',
          scrollX = TRUE,
          buttons = list(
            list(extend = "copy", text = "Copy", filename = "gsea",
                 exportOptions = list(
                   modifier = list(page = "all")
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


  }



