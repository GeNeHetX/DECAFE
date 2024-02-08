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
#if(length(new_packages)) install.packages(new_packages)
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
#HOME
output$pca_Image <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- "pca.png"
 
    # Return a list containing the filename and alt text
    list(src = filename,
      width = 800,
         
         alt = "overview RNA-Seq")

  }, deleteFile = FALSE)


output$anaDiff_Image <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- "analysediff.png"
 
    # Return a list containing the filename and alt text
    list(src = filename,
      width = 1000,
         
         alt = "overview RNA-Seq")

  }, deleteFile = FALSE)

output$rna_Image <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- "RNA_seq.png"
 
    # Return a list containing the filename and alt text
    list(src = filename,
      width = 800,
         
         alt = "overview RNA-Seq")

  }, deleteFile = FALSE)


output$gsea_Image <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- "gsea.png"
 
    # Return a list containing the filename and alt text
    list(src = filename,
      width = 1000,
         
         alt = "overview RNA-Seq")

  }, deleteFile = FALSE)



output$annot_Image <- renderImage({
    # When input$n is 3, filename is ./images/image3.jpeg
    filename <- "annot.png"
 
    # Return a list containing the filename and alt text
    list(src = filename,
      width = 1000,
         
         alt = "overview RNA-Seq")

  }, deleteFile = FALSE)


# Overview Panel
  annotFile <-reactive({
    req(input$'annot-file')
    annot = read.delim(input$'annot-file'$datapath,row.names=1)

  })
  
  countFile <-reactive({
    req(input$file)
    count= read.delim(input$file$datapath, sep='\t', row.names = 1, header=T,as.is=T)
  })
  annot_process <- reactive({

    annot = annotFile()
                                 
    annot$name = apply(annot,1,function(x) paste(colnames(annot),"_",x,collapse=','))
    annot2 = annot %>%
      mutate(
          all_list = as.list(strsplit(name, ","))
      )
     
    return(annot2)
  })
  
  output$barplot <- renderPlot(

    print(ggplot(data = annot_process(),mapping = aes(x = all_list)) +
      geom_bar(fill = '#262686') +
      scale_x_upset(order_by = "degree")+
      geom_text(stat='count', aes(label=..count..), vjust=2, size=5, color="white") +
      labs(x='Annotation combination', y='Distribution') +
      theme(axis.title = element_text(size=15))+
      theme_combmatrix(combmatrix.label.text = element_text(color ='#262686', size=12))) + 
      axis_combmatrix(levels = levels(as.factor(unlist(annot_process()$all_list))) 
      
    )  
    
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

  output$cond_1 <-renderUI({

    annot = annot_process()$name
    name = unique(annot)
    num  = c(1:length(unique(name)))

    choiceTable = data.frame(name, num)
    choix = setNames(as.numeric(choiceTable$num), choiceTable$name)

    
    return(selectInput("cond1", label = "Choose group 1",
                  choices = choix))
  })

  output$cond_2 <-renderUI({

    annot = annot_process()$name


    name = unique(annot)
    num  = c(1:length(unique(name)))

    choiceTable = data.frame(name, num)
    choix = setNames(as.numeric(choiceTable$num), choiceTable$name)

    
    return(selectInput("cond2", label = "Choose group 2",
                  choices = choix,selected=2))

  })

# Boxplot Panel
  output$geneTarget <-renderUI({
    req(input$file)
    count= countFile()
  
    gene = rownames(count)


    name = unique(gene)
    num  = c(1:length(unique(name)))

    choiceTable = data.frame(name, num)
    choix = setNames(as.numeric(choiceTable$num), choiceTable$name)

    
    return(selectInput("geneTarget", label = "Choose gene target to observe",
                  choices = choix))

  })



gene_target_plot <-reactive({
    req(input$'file')
    count = countFile()
    annot = annotFile()
    annot_name_cond1 = unique(annot_process()$name)[as.numeric(input$cond1)]
    annot_name_cond2 = unique(annot_process()$name)[as.numeric(input$cond2)]


    count_intersect = count[, intersect(rownames(annot), colnames(count))]


    annot_intersect= as.data.frame(annot[intersect(rownames(annot), colnames(count)), ])

     annot_intersect$condshiny = apply(annot_intersect,1,function(x) paste(colnames(annot_intersect),"_",x,collapse=','))

    dds <- DESeqDataSetFromMatrix(countData = count_intersect,
                                    colData = annot_intersect,
                                    design = ~condshiny)

    normalized_counts =  assay(vst(dds))

    annot_gp1 = rownames(annot_process())[which(annot_process()$name == annot_name_cond1)]
    annot_gp2 = rownames(annot_process())[which(annot_process()$name == annot_name_cond2)]


    norm_1 = normalized_counts[as.numeric(input$geneTarget),intersect(annot_gp1,colnames(normalized_counts))]
    norm_2 = normalized_counts[as.numeric(input$geneTarget),intersect(annot_gp2,colnames(normalized_counts))]

    norm1_rm = as.vector(norm_1)
    norm2_rm = as.vector(norm_2)

    res = data.frame( 
      count = c(norm1_rm,norm2_rm),
      condition = c(rep(annot_name_cond1,length(norm1_rm)), rep(annot_name_cond2,length(norm2_rm)))
    )
  return(res)
    })

output$bp_geneTarget <- renderPlot({
  res = gene_target_plot()
ggboxplot(res, 'condition','count',col='condition',shape='condition',xlab=NULL)+rotate_x_text(45)+scale_x_discrete(labels=NULL)
  

})

output$hist_geneTarget <- renderPlot({
  res = gene_target_plot()
  ggplot(res, aes(x=count, color=condition)) +
    geom_histogram(fill="white", alpha=0.5, position="identity") + 
    theme(legend.position="top")
  

  
})

#PCA
pcaVSD <- reactive({

  count = countFile()
annot = annot_process()

annot_name_cond1 = unique(annot_process()$name)[as.numeric(input$cond1)]
annot_name_cond2 = unique(annot_process()$name)[as.numeric(input$cond2)]
annot_gp1 = rownames(annot)[which(annot$name == annot_name_cond1)]
annot_gp2 = rownames(annot)[which(annot$name == annot_name_cond2)]

all_row = c(annot_gp1,annot_gp2)
annot = annot[all_row,]

count_intersect = count[, intersect(rownames(annot), colnames(count))]
annot_intersect= as.data.frame(annot[intersect(rownames(annot), colnames(count)), ])

dds <- DESeqDataSetFromMatrix(countData = count_intersect,
                                colData = annot_intersect,
                                design = ~name)



        dds =  vst(dds)
vst = assay(dds)

# center by gene 
vst = t(scale(t(vst), scale=F))

# Select most variable genes
gvar <- apply(vst, 1, sd)

mostvargenes <- order(gvar, decreasing=TRUE)[1:as.numeric(input$nb_gene)]

## Run PCA
# res_pca <- PCA(t(vst[mostvargenes,]), graph = F)
res_pca <- prcomp(t(vst[mostvargenes,]), scale. = TRUE)
coordinates= get_pca_ind(res_pca)$coord[,1:5]
res_pca$sample = rownames(annot_intersect)
pca_df = as.data.frame(cbind(coordinates,rownames(annot_intersect),annot_intersect$name))
print(head(pca_df))
colnames(pca_df)=c(paste0("Dim",1:5),"Sample","Group")
eig = get_eig(res_pca)[1:5,2]
print(eig)

  return(list(pca_ind=pca_df,pca=res_pca,eig =eig ,name=annot_intersect$name, name_samp = rownames(annot_intersect),name_gene= rownames(count_intersect)))

})

pcaGraph <- reactive({
  data= pcaVSD()$pca_ind
  eig = round(pcaVSD()$eig,2)
  data_sub = data[c(as.numeric(input$dim1),as.numeric(input$dim2),6,7)]
  colnames(data_sub)= c("x","y","Sample","Group")
  data_sub$x = as.numeric(data_sub$x)
  data_sub$y = as.numeric(data_sub$y)
  plot=ggplot(data=data_sub, aes(x=x, y=y, col=Group, tooltip=Sample)) +
    geom_point(size=1) + theme_minimal() +
    labs(x = paste0("Dimension ", input$dim1," (",eig[as.numeric(input$dim1)],"%)"),y = paste0("Dimension ", input$dim2," (",eig[as.numeric(input$dim2)],"%)")) + 
   #scale_x_continuous(limits = c(-max(data_sub$x),max(data_sub$x))) +
    #scale_y_continuous(limits = c(-max(data_sub$y),max(data_sub$y))) +
    geom_vline(xintercept=0, linetype="dashed", color = "grey") +
    geom_hline(yintercept=0, linetype="dashed", color = "grey") +
    theme(legend.position="bottom")
  return(plot)
})

output$pcaVsd <- renderPlotly({ ggplotly(pcaGraph()) %>% 
    layout(legend = list(orientation = 'h', x = 0.45, y = 1.1)) })


sample1 <- reactive({
  pca = pcaVSD()$pca
  
  contrib=get_pca_ind(pca)$contrib[,as.numeric(input$dim1)]
  tab = as.data.frame(cbind(name=pcaVSD()$name_samp,contrib=contrib))

  tab= tab[order(abs(as.numeric(tab$contrib)),decreasing=TRUE),]
  rownames(tab)= NULL

  return(tab[1:min(10,nrow(tab)),])

})

 output$sampleTop1 <- DT::renderDT(server = FALSE, {
    DT::datatable(
      sample1(),
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

sample2 <- reactive({
  pca = pcaVSD()$pca
  contrib=get_pca_ind(pca)$contrib[,as.numeric(input$dim2)]
  tab = as.data.frame(cbind(name=pcaVSD()$name_samp,contrib=contrib))

  tab= tab[order(abs(as.numeric(tab$contrib)),decreasing=TRUE),]
  rownames(tab)= NULL

  return(tab[1:min(10,nrow(tab)),])

})

 output$sampleTop2 <- DT::renderDT(server = FALSE, {
    DT::datatable(
      sample2(),
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
 gene1 <- reactive({
  pca = pcaVSD()$pca
  contrib=get_pca_var(pca)$contrib[,as.numeric(input$dim1)]
  tab = as.data.frame(cbind(name=pcaVSD()$name_gene,contrib=contrib))

  tab= tab[order(abs(as.numeric(tab$contrib)),decreasing=TRUE),]
  rownames(tab)= NULL

  return(tab[1:min(10,nrow(tab)),])

})

 output$geneTop1 <- DT::renderDT(server = FALSE, {
    DT::datatable(
      gene1(),
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
gene2 <- reactive({
  pca = pcaVSD()$pca
  contrib=get_pca_var(pca)$contrib[,as.numeric(input$dim2)]
  tab = as.data.frame(cbind(name=pcaVSD()$name_gen,contrib=contrib))

  tab= tab[order(abs(as.numeric(tab$contrib)),decreasing=TRUE),]
  rownames(tab)= NULL

  return(tab[1:min(10,nrow(tab)),])

})

 output$geneTop2 <- DT::renderDT(server = FALSE, {
    DT::datatable(
      gene2(),
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

res_deseq <- reactive({
  req(input$file)
count = countFile()
annot = annot_process()

annot_name_cond1 = unique(annot_process()$name)[as.numeric(input$cond1)]
annot_name_cond2 = unique(annot_process()$name)[as.numeric(input$cond2)]
annot_gp1 = rownames(annot)[which(annot$name == annot_name_cond1)]
annot_gp2 = rownames(annot)[which(annot$name == annot_name_cond2)]

all_row = c(annot_gp1,annot_gp2)
annot = annot[all_row,]

count_intersect = count[, intersect(rownames(annot), colnames(count))]
annot_intersect= as.data.frame(annot[intersect(rownames(annot), colnames(count)), ])
zero_threshold = as.numeric(input$zero_threshold)
countfilt = count_intersect[rowMeans(count_intersect == 0) <= (zero_threshold ), ]

dds = DESeqDataSetFromMatrix(countData = countfilt,
                                colData = annot_intersect,
                                design = ~name)
dds = dds[rowSums(counts(dds)) >= 10]
    dds = DESeq(dds)
    table = results(dds)
    table$name = rownames(table)
    table = as.data.frame(table)
    table = table[order(abs(table$stat),decreasing=TRUE),]

   return(list(deseq = dds, res = table) )

   })

volcano <- reactive({

table = res_deseq()$res
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
      res_deseq()$res,
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
  table = res_deseq()$res
  vec = table$stat
  names(vec) = table$name

  res = fgsea(pathways,vec)


  lE_list = res$leadingEdge

  
  lE_vector = sapply(lE_list, paste, collapse=",")
  res$leadingEdge = lE_vector
  res = as.data.frame(na.omit(res))
  res_sig = res[which(as.numeric(res$padj) < 0.05),]
  res_sort = res_sig[order(abs(as.numeric(res_sig$NES)),decreasing=TRUE),]
  return(res_sort)
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



