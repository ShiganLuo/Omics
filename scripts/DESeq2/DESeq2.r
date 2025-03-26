suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(argparse))
DESeq2Analysis = function(dat,resfile="",normMethods="cpm",sampleNumber){
    cat("\nDESeq2Analysis function\n")
    print(paste0("outfile path : ",resfile))
    print(paste0("the methods of normalizing raw counts : ",normMethods))

    "DESeq2Analysis function : use the DESeq2 to perform differential gene analysis
     return results and dds
    "
    ####### DESeq2 analysis
    #prepare Step : Specify the order of grouping factors
    #Note that the sample order in the expression matrix and the grouping order here must correspond one to one
    coldata <- data.frame(condition = factor(rep(c('control', 'experiment'), each = sampleNumber), levels = c('control', 'experiment')))
    
    #Step 1 : build a DESeqDataSet object
    dds <- DESeqDataSetFromMatrix(countData = dat, colData = coldata, design= ~condition)

    #Step 2 : Calculate the difference fold and obtain the p value
    #Note: parallel = TRUE can run in multiple threads. It is recommended to enable it when the amount of data is large.
    dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE)

    # Step 3 : get the result
    #Note : that treat should be placed first and control should be placed second, 
    #which means which genes are upregulated/downregulated in treat compared to control

    # Step 4 : normalize raw counts
    if ( normMethods == "DESeq2"){
      #normalized_counts is matrix
      normalized_counts <- round(counts(dds1, normalized = TRUE),3)
    } else if ( normMethods == "cpm" ){
      raw_counts <- counts(dds1, normalized = F)
      # y muste be a DGEList object
      y <- edgeR::DGEList(counts = raw_counts)
      y <- edgeR::calcNormFactors(y, method="TMM")
      #prior.count prevents a gene from being expressed in all samples and causing a division by 0 error
      #normalized_counts is matrix
      normalized_counts <- edgeR::cpm(y,normalized.lib.sizes = TRUE,log = FALSE,prior.count=2) 
    } else {
       stop("Invalid normMethods (should be DESeq2 or cpm)")
    }
      # print(head(normalized_counts))
      # print(class(normalized_counts))
    normalized_counts = as.data.frame(normalized_counts)
    res <- results(dds1, contrast = c('condition', 'experiment', 'control'))
    res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
    res2 <- cbind(res1,normalized_counts)
      # print(head(res2))
    if(resfile != ""){
        write.table(res2, outfile, col.names = NA, sep = '\t', quote = FALSE)
    }
    return(list(res2,dds1))
}

TEFilter = function(res){
    "
    TEFilter function: filter the TEcount data
    return datatframe only contain TE row
    "
  df_TE <- res[grepl("^[^:]+:[^:]+:[^:]+$", rownames(res)), ]
  return(df_TE)
}
ScreenFeature = function(res1,upfile="",downfile="",updownfile = "",padj_cutoff = 0.05){
    "
    ScreenFeature function: get Upregulated and Downregulated genes
    retrun list(res1_up,res1_down)
    "
    print(paste0("outpath of Upregulated genes: ",upfile))
    print(paste0("outpath of Downregulated genes: ",downfile))
    print(paste0("outpath of disordered genes: ",updownfile))
    ##Screen differentially expressed genes
    #First sort the table, sort by padj value in ascending order, and continue to sort by log2FC in descending order for the same padj value
    res1 <- res1[order(res1$padj, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

    #log2FC≥1 & padj<0.01 indicates up, representing significantly upregulated genes
    #log2FC≤-1 & padj<0.01 indicates down, representing significantly downregulated genes
    #Others indicate none, representing non-differential genes
    
    res1[which(res1$log2FoldChange >= 0.58 & res1$padj < padj_cutoff),'sig'] <- 'up'
    res1[which(res1$log2FoldChange <= -0.58 & res1$padj < padj_cutoff),'sig'] <- 'down'
    res1[which(abs(res1$log2FoldChange) <= 0.58 | res1$padj >= padj_cutoff),'sig'] <- 'none'
    col_order <- colnames(res1)
    # Remove "sig" first, then insert "sig" after "padj"
    col_order <- append(setdiff(col_order, "sig"), "sig", after = match("padj", col_order))
    res1 <- res1[, col_order]

   #Separate output according to up and down
    res1_up <- subset(res1, sig == 'up')
    res1_down <- subset(res1, sig == 'down')
    if ( upfile != "" ){
      write.table(res1_up, file = upfile, sep = '\t', col.names = NA, quote = FALSE)
    }
    if ( downfile != ""){
      write.table(res1_down, file = downfile, sep = '\t', col.names = NA, quote = FALSE)
    }
    if (updownfile != ""){
      res1_updown = subset(res1, sig == 'up' | sig == 'down')
      write.table(res1_updown, file = updownfile, sep = '\t', col.names = NA, quote = FALSE)
    }
    return(list(res1_up,res1_down))
}

updown_heatmap <- function(res, outjpeg,coldata=NULL) {
    cat("\nupdown_heatmap function\n")
    # generate the color palette
    brewer_palette <- "RdBu"
    ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer_palette))
    mr <- ramp(256)[256:1]
    # obtain the significant genes and order by log2FoldChange
    heatmap_values <- as.matrix(res[,-c(1:7)])
    # print(heatmap_values)
    # plot the heatmap using pheatmap
    ngenes = length(heatmap_values)
    if (is.null(coldata)) {
      print("coldata 为空，将不添加注释信息")
      jpeg(width = 150, height = 200, units = "mm", res = 300, file = outjpeg)
        p = pheatmap::pheatmap(heatmap_values, color = mr, scale = "row", fontsize_col = 10, 
                          fontsize_row = 200/ngenes, fontsize = 5, border_color = NA)
      dev.off()
    } else {
      jpeg(width = 150, height = 200, units = "mm", res = 300, file = outjpeg)
        p = pheatmap::pheatmap(heatmap_values, color = mr, scale = "row", annotation_col = coldata,
                          fontsize_col = 10, fontsize_row = 200/ngenes, fontsize = 5, border_color = NA)
      dev.off()
    }
    
    
}

volcano <- function (res, outjpeg,mode = "gene",geneAnnotation = "", nlabel = 10, label.by = "padj"){
  # add gene_name column to the results table
  if ( mode == 'TE'){
    res = res %>% 
    rownames_to_column(var = "index") %>%
    mutate(gene_name = word(index, 1, sep = ":"))
  } else if ( mode == 'gene' && geneAnnotation != ""){
    df <- read.csv(geneAnnotation, sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
    # it is gene mode as long as it contains ensembl_id ,this can resolve that both str1:str2:…… ensembl_id
    res <- res %>%
      rownames_to_column(var = "index") %>%
      mutate(gene_name = ifelse(
      str_count(index, ":") >= 2, 
      str_extract(index, "^[^:]+"),  # else extract gene_name from df
      df$gene_name[match(index, rownames(df))]))
  } else if ( mode == 'gene' && geneAnnotation == ""){
    stop("Invalid geneAnnotation argument. Please provide a gene annotation file, because you choose mode gene.")
  } else {
    stop("Invalid mode argument. Choose either gene or TE.")
  }

    # print(head(res))
  # assign significance to results based on padj
  res <- res %>%
  mutate(
    gene_status = case_when(
      log2FoldChange > 0.58 & -log10(padj) > 1.3 ~ "Upregulated genes",
      log2FoldChange < -0.58 & -log10(padj) > 1.3 ~ "Downregulated genes",
      TRUE ~ "Other genes"
    )
  )
  res = res[!is.na(res$padj),]
  # Calculate the number of each gene_status
  status_counts <- res %>%
    group_by(gene_status) %>%
    summarise(count = n())

  # Create new legend labels containing quantities
  new_labels <- setNames(
    sapply(status_counts$gene_status, function(status) {
      count <- status_counts$count[status_counts$gene_status == status]
      paste0(status, " (n = ", count, ")")
    }),
    status_counts$gene_status
  )
  # print(head(res))
  
  # get labels for the highest or lowest genes according to either padj or log2FoldChange
  if (label.by == "padj") {
    up_genes <- res %>%
      filter(gene_status == "Upregulated genes") %>%
      arrange(padj) %>%
      head(nlabel)
    down_genes <- res %>%
      filter(gene_status == "Downregulated genes") %>%
      arrange(padj) %>%
      head(nlabel)
  } else if (label.by == "log2FoldChange") {
    up_genes <- head(arrange(res[res[,"gene_status"] == "Upregulated genes" ,], desc(log2FoldChange)),nlabel)
    down_genes <- head(arrange(res[res[,"gene_status"] == "Upregulated genes" ,], log2FoldChange),nlabel)
  } else
    stop ("Invalid label.by argument. Choose either padj or log2FoldChange.")
  
    p = ggplot(res, aes(log2FoldChange, -log10(padj))) +
      geom_point(aes(col=gene_status)) + 
      # scale_fill_manual(values=c("Upregulated genes" ="red","Downregulated genes" = "blue","Other genes" = "black"), 
      #                       name = "gene classification", 
      #                       labels = c("Upregulated genes" = uplegendLabel, "Downregulated genes" = downlegendLabel,"Other genes" = otherlegendLabel))+
      scale_color_manual(values=c("Upregulated genes" ="red","Downregulated genes" = "blue","Other genes" = "black"),
        labels = new_labels) + 
      ggrepel::geom_text_repel(data=up_genes, aes(label=head(gene_name,nlabel)), color = "#F59494", size = 3)+
      ggrepel::geom_text_repel(data=down_genes, aes(label=head(gene_name,nlabel)), color = "#93ACF6", size = 3)+
      labs ( x = expression(log[2]("FoldChange")), y = expression(-log[10]("adjusted p-value")))+
      geom_vline(xintercept = 0.58, linetype = "dotted")+
      geom_vline(xintercept = -0.58, linetype = "dotted")+
      geom_hline(yintercept = 1.3, linetype = "dotted")+
      theme_classic() +  
      theme(
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA,linewidth = 0.2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()   
      )
  jpeg(width = 200, height = 150, units = "mm", res = 300, file = outjpeg)
    print(p)
  dev.off()
}
MAplot = function (res, outjpeg,mode = "gene",geneAnnotation = "", nlabel = 10, label.by = "padj"){
  # add gene_name column to the results table
  if ( mode == 'TE'){
    res = res %>% 
    rownames_to_column(var = "index") %>%
    mutate(gene_name = word(index, 1, sep = ":"))
  } else if ( mode == 'gene' && geneAnnotation != ""){
    df <- read.csv(geneAnnotation, sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
    # it is gene mode as long as it contains ensembl_id ,this can resolve that both str1:str2:…… ensembl_id
    res <- res %>%
      rownames_to_column(var = "index") %>%
      mutate(gene_name = ifelse(
      str_count(index, ":") >= 2, 
      str_extract(index, "^[^:]+"),  # else extract gene_name from df
      df$gene_name[match(index, rownames(df))]))
  } else if ( mode == 'gene' && geneAnnotation == ""){
    stop("Invalid geneAnnotation argument. Please provide a gene annotation file, because you choose mode gene.")
  } else {
    stop("Invalid mode argument. Choose either gene or TE.")
  }

    # print(head(res))
  # assign significance to results based on padj
  res <- res %>%
  mutate(
    gene_status = case_when(
      log2FoldChange > 0.58 & -log10(padj) > 1.3 ~ "Upregulated genes",
      log2FoldChange < -0.58 & -log10(padj) > 1.3 ~ "Downregulated genes",
      TRUE ~ "Other genes"
    )
  )
  res = res[!is.na(res$padj),]
  # Calculate the number of each gene_status
  status_counts <- res %>%
    group_by(gene_status) %>%
    summarise(count = n())

  # Create new legend labels containing quantities
  new_labels <- setNames(
    sapply(status_counts$gene_status, function(status) {
      count <- status_counts$count[status_counts$gene_status == status]
      paste0(status, " (n = ", count, ")")
    }),
    status_counts$gene_status
  )
  # print(head(res))
  
  # get labels for the highest or lowest genes according to either padj or log2FoldChange
  if (label.by == "padj") {
    up_genes <- res %>%
      filter(gene_status == "Upregulated genes") %>%
      arrange(padj) %>%
      head(nlabel)
    down_genes <- res %>%
      filter(gene_status == "Downregulated genes") %>%
      arrange(padj) %>%
      head(nlabel)
  } else if (label.by == "log2FoldChange") {
    up_genes <- head(arrange(res[res[,"gene_status"] == "Upregulated genes" ,], desc(log2FoldChange)),nlabel)
    down_genes <- head(arrange(res[res[,"gene_status"] == "Upregulated genes" ,], log2FoldChange),nlabel)
  } else
    stop ("Invalid label.by argument. Choose either padj or log2FoldChange.")
  
    p = ggplot(res, aes(log2(baseMean), log2FoldChange)) +
      geom_point(aes(col=gene_status)) + 
      # scale_fill_manual(values=c("Upregulated genes" ="red","Downregulated genes" = "blue","Other genes" = "black"), 
      #                       name = "gene classification", 
      #                       labels = c("Upregulated genes" = uplegendLabel, "Downregulated genes" = downlegendLabel,"Other genes" = otherlegendLabel))+
      scale_color_manual(values=c("Upregulated genes" ="red","Downregulated genes" = "blue","Other genes" = "grey"),
        labels = new_labels) + 
      ggrepel::geom_text_repel(data=up_genes, aes(label=head(gene_name,nlabel)), color = "#F59494", size = 3)+
      ggrepel::geom_text_repel(data=down_genes, aes(label=head(gene_name,nlabel)), color = "#93ACF6", size = 3)+
      labs ( x = expression(log[2]("Expression in baseMean")), y = expression(log[2]("FoldChange")))+
      geom_hline(yintercept = 0.58, linetype = "dotted")+
      geom_hline(yintercept = -0.58, linetype = "dotted")+
      theme_classic() +  
      theme(
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA,linewidth = 0.2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()   
      )
  jpeg(width = 200, height = 150, units = "mm", res = 300, file = outjpeg)
    print(p)
  dev.off()
}

plot_PCANorm = function(res,outjpeg,sampleNumber,controName,experimentName){
  cpm = res[,-c(1:6)]
  # print(head(cpm))
  pca_result <- prcomp(t(cpm))
  # print(pca_result)
  # 步骤3：提取主成分信息
  pca_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Sample = colnames(cpm)
  )
  # 假设样本有分组信息，这里简单模拟
  pca_data$Group <- factor(rep(c(controlName, experimentName), each = sampleNumber))

  # 步骤4：绘制PCA图
  p = ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
    geom_point(size = 3) +
    geom_text(hjust = 0.5, vjust = -0.5,show.legend = FALSE,size = 0.5) +
    labs(
      x = paste0("PC1: ", round(pca_result$sdev[1]^2 / sum(pca_result$sdev^2) * 100, 2), "% variance"),
      y = paste0("PC2: ", round(pca_result$sdev[2]^2 / sum(pca_result$sdev^2) * 100, 2), "% variance"),
      title = "PCA Plot of CPM Values"
    ) + 
  theme_minimal() +
  theme(
    panel.border = element_blank(),  # 去除整个面板边框
    axis.line = element_line(colour = "black"),  # 重新添加坐标轴线条
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank(),  # 去除次要网格线
    axis.ticks = element_line(colour = "black"),  # 设置坐标轴刻度线颜色
    axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
    plot.title = element_text(hjust = 0.5),  # 标题居中
    axis.title = element_text(colour = "black")  # 设置坐标轴标题颜色
  )
  jpeg(width = 200, height = 150, units = "mm", res = 300, file = outjpeg)
    print(p)
  dev.off()
}

TEcount_read = function(infile,control,experiment){
    ###### read format
    # TEcount_read function: read the output of TEcount
    # control and experiment: c() ->
    # decide which column is control,which column is experiment,and sort the dataframe according to it

    ###Differential analysis of TEs with different genetic backgrounds
    TEcount = read.csv(infile,sep="\t",row.names=1,header=TRUE)
    TEcount = TEcount[,c(control,experiment)] # sampe with DESeq2 condition

    ### only TE
    df_TE <- TEcount[grepl("^[^:]+:[^:]+:[^:]+$", rownames(TEcount)), ]
    df_TE =  df_TE[rowMeans(df_TE) > 5, ] # the expression of feature must > specific number,make plot beautiful
    ### only Gene
    df_Gene <- TEcount[!grepl("^[^:]+:[^:]+:[^:]+$", rownames(TEcount)), ]
    df_Gene = df_Gene[rowMeans(df_TE) > 5, ] # the expression of feature must > specific number,make plot beautiful

    return(list(Gene_TE = TEcount,TE = df_TE, Gene = df_Gene))
}
TElocal_read = function(infile,control,experiment){
    ###### read format
    # TElocal_read function: read the output of TEcount
    # control and experiment: c() ->
    # decide which column is control,which column is experiment,and sort the dataframe according to it
    TElocal = read.csv(infile,sep="\t",row.names=1,header=TRUE)
    TElocal = TElocal[,c(control,experiment)] # sampe with DESeq2 condition
    ### only TE
    df_TE <- TElocal[grepl("^[^:]+:[^:]+:[^:]+:[^:]+$", rownames(TElocal)), ]
    df_TE = df_TE[rowMeans(df_TE) > 5,]
    return(df_TE)
}
#############################################################################################################

parser <- ArgumentParser(description='DESeq2 analysis for TEcount and TElocal')
parser$add_argument('-m', '--mode', type='character', required=TRUE,choices = c("TEcount", "TElocal"),
                    help='mode: TEcount,TElocal,All')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='input count matrix. Rows are genes. Columns are samples.')
parser$add_argument('-g', '--group', type='character', required=TRUE,
                    help='path to group file,sep with table, two columns: sample,group')
parser$add_argument('-p', '--pattern', type='character', required=TRUE,nargs="+",
                    help='DESeq2 pattern, --pattern control experiment')           
parser$add_argument('-f', '--figure', type='character', required=TRUE,nargs="+",choices = c("heatmap","volcano","pca"),
                    help='normalize method: "TMM","TMMwsp","RLE","upperquartile","none", default: TMM')
parser$add_argument('-o', '--outdir', type='character', required=TRUE,
                    help='output directory')
parser$add_argument('-a', '--annotation', type='character', required=TRUE,
                    help='annotation file, three column gene_id,gene_type,gene_name')
parser$add_argument('-Tcm', '--TEcountMode', type='character',nargs="+",choices = c("all","Gene_TE", "TE", "Gene"),
                    help='TEcount mode: Gene_TE,TE,Gene or all')
args <- parser$parse_args()

mode = args$mode
infile = args$matrix
group = args$group
outdir = args$outdir
pattern = args$pattern
figure = args$figure
annotationG = args$annotation
TEcountMode = args$TEcountMode
print("--------------Start DESeq2 analysis--------------")
print(paste("mode:",mode))
print(paste("infile:",infile))
print(paste("group:",group))
print(paste("outdir:",outdir))
cat("pattern:",pattern,"\n")
cat("figure:",figure,"\n")
print(paste("annotation:",annotationG))
cat("TEcountMode:",TEcountMode,"\n")

paraserPattern = function(pattern,group){
  df_group = read.csv(group,sep="\t",header=TRUE)
  controlStr = pattern[1]
  experimentStr = pattern[2]
  control = as.vector(df_group$sample[df_group$group == controlStr])
  experiment = as.vector(df_group$sample[df_group$group == experimentStr])
  if (length(control) != length(experiment)){
    stop("control and experiment samples are not equal")
  }
  ### for heatmap
  sample_anno = subset(df_group, group %in% c(controlStr, experimentStr))
  rownames(sample_anno) = sample_anno$sample
  sample_anno$sample = NULL
  sample_anno$group = factor(sample_anno$group)
  return(list(control,experiment,sample_anno))
}
ce = paraserPattern(pattern,group)
control = ce[[1]]
experiment = ce[[2]]
sample_annoQ = ce[[3]]
sampleNumberQ = length(control)
cat("control:",control,"\n")
cat("experiment:",experiment,"\n")
cat("sampleNumber:",sampleNumberQ,"\n")
controlName = pattern[1]
experimentName = pattern[2]

plot_heatmap = function(res,mode,sample_anno = sample_annoQ,outdir = args$outdir,dataType=""){
    print(sample_anno)
    if (!dir.exists(paste(outdir,"DESeq2/upDown/",sep=""))){
        dir.create(paste(outdir,"DESeq2/upDown/",sep=""),recursive = TRUE,showWarnings = FALSE)
    }
    if (!dir.exists(paste(outdir,"DESeq2/heatmap/",sep=""))){
        dir.create(paste(outdir,"DESeq2/heatmap/",sep=""),recursive = TRUE,showWarnings = FALSE)
    }
    if ( mode == "TEcount"){
      upfile = paste(outdir,"DESeq2/upDown/TEcount_",dataType,"_up.csv",sep="")
      downfile = paste(outdir,"DESeq2/upDown/TEcount_",dataType,"_down.csv",sep="")
      updownfile = paste(outdir,"DESeq2/upDown/TEcount_",dataType,"_updown.csv",sep="")
      updown = ScreenFeature(res,upfile,downfile,updownfile)
      df_up = updown[[1]]
      df_down = updown[[2]]
      df = rbind(df_up,df_down)
      if ( nrow(df_up) >= 2) {
        updown_heatmap(df_up,outjpeg=paste(outdir,"DESeq2/heatmap/TEcount_",dataType,"_up.jpeg",sep=""),coldata = sample_anno)
      } else {
        print("上调基因数量小于2,无法绘制热图")
      }
      if ( nrow(df_down) >= 2) {
        updown_heatmap(df_down,outjpeg=paste(outdir,"DESeq2/heatmap/TEcount_",dataType,"_down.jpeg",sep=""),coldata = sample_anno)
        
      } else {
        print("下调基因数量小于2,无法绘制热图")
      }
      if ( nrow(df) >= 2) {
        updown_heatmap(df,outjpeg=paste(outdir,"DESeq2/heatmap/TEcount_",dataType,"_updown.jpeg",sep=""),coldata = sample_anno)
        
      } else {
        print("失调基因数量小于2,无法绘制热图")
      }
    } else if ( mode == "TElocal" ){
        upfile = paste(outdir,"DESeq2/upDown/TElocal_TE_up.csv",sep = "")
        downfile = paste(outdir,"DESeq2/upDown/TElocal_TE_down.csv",sep = "")
        updownfile = paste(outdir,"DESeq2/upDown/TElocal_TE_updown.csv",sep = "")
        updown = ScreenFeature(res,upfile,downfile,updownfile)
        df_up = updown[[1]]
        df_down = updown[[2]]
        df = rbind(df_up,df_down)
        if ( nrow(df_up) >= 2) {
          updown_heatmap(df_up,outjpeg=paste(outdir,"DESeq2/heatmap/TElocal_up.jpeg",sep=""),coldata = sample_anno)
        } else {
          print("上调基因数量小于2,无法绘制热图")
        }
        if ( nrow(df_down) >= 2) {
          updown_heatmap(df_down,outjpeg=paste(outdir,"DESeq2/heatmap/TElocal_down.jpeg",sep=""),coldata = sample_anno) 
        } else {
          print("下调基因数量小于2,无法绘制热图")
        }
        if ( nrow(df) >= 2) {
          updown_heatmap(df,outjpeg=paste(outdir,"DESeq2/heatmap/TElocal_updown.jpeg",sep=""),coldata = sample_anno)
        } else {
          print("失调基因数量小于2,无法绘制热图")
        }
    } else {
      stop("mode error: During heatmap creation, please enter the correct mode, TEcount or TElocal")
    }

}
plot_volcano = function(mode,res,dataType = "",annotation = annotationG){
    if ( !dir.exists(paste(outdir,"DESeq2/volcano/",sep=""))){
        dir.create(paste(outdir,"DESeq2/volcano/",sep=""),recursive = TRUE,showWarnings = FALSE)
    }
    if ( mode == "TEcount"){
      outjpegV = paste(outdir, "DESeq2/volcano/TEcount_",dataType,"_volcano.jpeg",sep="")
      outjpegM = paste(outdir, "DESeq2/volcano/TEcount_",dataType,"_MAplot.jpeg",sep="")
      if (grepl("Gene", dataType) ){
          volcano(res,outjpegV,mode = 'gene',geneAnnotation = annotation,nlabel = 10, label.by = "padj")
          MAplot(res,outjpegM,mode = 'gene',geneAnnotation = annotation,nlabel = 10, label.by = "padj")
      } else {
          volcano(res,outjpegV,mode = 'TE',nlabel = 10, label.by = "padj")
          MAplot(res,outjpegM,mode = 'TE',nlabel = 10, label.by = "padj")
      }
    } else if ( mode == "TElocal") {
        outjpegV = paste(outdir ,"DESeq2/volcano/TElocal_volcano.jpeg",sep="")
        volcano(res,outjpegV,mode = 'TE',nlabel = 10, label.by = "padj")
        outjpegM = paste(outdir , "DESeq2/volcano/TElocal_MAplot.jpeg",sep="")
        MAplot(res,outjpegM,mode = 'TE',nlabel = 10, label.by = "padj")
    } else {
      stop("mode error: During heatmap creation, please enter the correct mode, TEcount or TElocal")
    }


}

if (args$mode == "TEcount"){
  dfList = TEcount_read(infile,control,experiment)
  if ( "pca" %in% figure){
    df = dfList[["Gene"]]
      res = DESeq2Analysis(df,normMethods = "cpm",sampleNumber = sampleNumberQ)[[1]]
      if(!dir.exists(paste(outdir,"DESeq2/plot/",sep=""))){
          dir.create(paste(outdir,"DESeq2/plot/",sep=""),recursive = TRUE,showWarnings = FALSE)
      }
      outjpeg = paste(outdir , "DESeq2/plot/cpmPCA.jpeg",sep = "")
      plot_PCANorm(res,outjpeg,sampleNumberQ,controlName,experimentName)
  } 
  if ( "all" %in% TEcountMode){
  ## TEcount_read generate three dataframe: Gene_TE,TE,Gene
    for ( dataType in names(dfList)){
      df = dfList[[dataType]]
      ## DESeq2Analysis generate a list, first element is res(dataframe combine counts), second element is dds(DESeqDataSet object)
      outfile = paste(outdir,"DESeq2/TEcount_",dataType,".csv",sep="")
      resDds = DESeq2Analysis(df,resfile = outfile,normMethods = "DESeq2", sampleNumber = sampleNumberQ)
      if ( dataType == "Gene_TE" ){
        res1 = resDds[[1]]
        dds1 = resDds[[2]]
        res1 = TEFilter(res1)
          if ( "heatmap" %in% figure){
            plot_heatmap(res = res1,mode = "TEcount",dataType = dataType)
          }
          if ( "volcano" %in% figure){
            plot_volcano(mode = "TEcount",res = res1,dataType = dataType)
          }
      } else {
      res1 = resDds[[1]]
      dds1 = resDds[[2]]
      if ( "heatmap" %in% figure){
        plot_heatmap(res = res1,mode = "TEcount",dataType = dataType)
      }
      if ( "volcano" %in% figure){
        plot_volcano(mode = "TEcount",res = res1,dataType = dataType)
      }
      }
    }
  # !("all" %in% TEcountMode)避免重复执行，浪费资源
  } else if ( "Gene_TE" %in% TEcountMode && !("all" %in% TEcountMode) ){
    df = dfList[["Gene_TE"]]
    ## DESeq2Analysis generate a list, first element is res(dataframe combine counts), second element is dds(DESeqDataSet object)
    outfile = paste(outdir,"DESeq2/TEcount_Gene_TE.csv",sep="")
    resDds = DESeq2Analysis(df,resfile = outfile,normMethods = "DESeq2", sampleNumber = sampleNumberQ)
    res1 = resDds[[1]]
    dds1 = resDds[[2]]
    res1 = TEFilter(res1)
    if ( "heatmap" %in% figure){
      plot_heatmap(res = res1,mode = "TEcount",dataType = "Gene_TE")
    }
    if ( "volcano" %in% figure){
      plot_volcano(mode = "TEcount",res = res1,dataType = "Gene_TE")
    }
  #!("all" %in% TEcountMode)避免重复执行，浪费资源
  } else if ( "TE" %in% TEcountMode && !("all" %in% TEcountMode) ){
    df = dfList[["TE"]]
    ## DESeq2Analysis generate a list, first element is res(dataframe combine counts), second element is dds(DESeqDataSet object)
    outfile = paste(outdir,"DESeq2/TEcount_TE.csv",sep="")
    resDds = DESeq2Analysis(df,resfile = outfile,normMethods = "DESeq2", sampleNumber = sampleNumberQ)
    res1 = resDds[[1]]
    dds1 = resDds[[2]]
    if ( "heatmap" %in% figure){
      plot_heatmap(res = res1,mode = "TEcount",dataType = "TE")
    }
    if ( "volcano" %in% figure){
      plot_volcano(mode = "TEcount",res = res1,dataType = "TE")
    }
  #!("all" %in% TEcountMode)避免重复执行，浪费资源
  } else if ( "Gene" %in% TEcountMode && !("all" %in% TEcountMode) ){
    df = dfList[["Gene"]]
    ## DESeq2Analysis generate a list, first element is res(dataframe combine counts), second element is dds(DESeqDataSet object)
    outfile = paste(outdir,"DESeq2/TEcount_Gene.csv",sep="")
    resDds = DESeq2Analysis(df,resfile = outfile,normMethods = "DESeq2", sampleNumber = sampleNumberQ)
    res1 = resDds[[1]]
    dds = resDds[[2]]
    if ( "heatmap" %in% figure){
      plot_heatmap(res = res1,mode = "TEcount",dataType = "Gene")
    }
    if ( "volcano" %in% figure){
      plot_volcano(mode = "TEcount",res = res1,dataType = "Gene")
    }
  } else {
    stop("请输入正确的TEcountMode,Gene_TE,TE,Gene or all")
  }

} else if ( args$mode == "TElocal"){
  ## TElocal_read generate a dataframe: TEsite
  df = TElocal_read(infile,control,experiment)
  ## DESeq2Analysis generate a list, first element is res(dataframe combine counts), second element is dds(DESeqDataSet object)
  outfile = paste(outdir,"DESeq2/TElocal_TE.csv",sep="")
  resDds = DESeq2Analysis(df,resfile = outfile,sampleNumber = sampleNumberQ)
  res1 = resDds[[1]]
  dds1 = resDds[[2]]
  if ( "heatmap" %in% figure){
    plot_heatmap(res = res1,mode = "TElocal")
  }
  if ( "volcano" %in% figure){
    plot_volcano(mode = "TElocal",res = res1)
  }

}





