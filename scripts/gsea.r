suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))
prepare_ranked_list <- function(ranked_list) { 
  # if duplicate gene names present, average the values
  if( sum(duplicated(ranked_list$gene_name)) > 0) {
    ranked_list <- aggregate(.~gene_name, FUN = mean, data = ranked_list)
    ranked_list <- ranked_list[order(ranked_list$log2FoldChange, decreasing = T),]
    # print("--------")
  }
    # print(ranked_list)
  # omit rows with NA values

  ranked_list <- na.omit(ranked_list)
  # This will remove duplicate values ​​without changing the overall sorting.
  set.seed(123)
  ranked_list$log2FoldChange <- ranked_list$log2FoldChange + runif(nrow(ranked_list), -1e-5, 1e-5)
  # turn the dataframe into a named vector
  ranked_list <- tibble::deframe(ranked_list)
    # print(ranked_list)
}

GSEA_prepare = function(res,mode = "gene",AnnotationFile="",outfile=""){
  "GSEA_prepare function: prepare the data for GSEA analysis
  AnnotatilFile : gene_id,gene_type,gene_name. gmt file is symbols based" 
  if ( mode == "Gene"){
    if ( AnnotationFile != ""){
      df <- read.csv(AnnotationFile, sep = "\t", header = T, stringsAsFactors = F)
      res <- res %>% rownames_to_column(var = "gene_id") %>% left_join(df, by = "gene_id")
      res_prot <- res[which(res$gene_type == "protein_coding"),] # stanard gmt only contains protein coding genes
      res_prot_ranked <- res_prot[order(res_prot$log2FoldChange, decreasing = T),c("gene_name", "log2FoldChange")]
      res_prot_ranked <- na.omit(res_prot_ranked)
    } else {
      stop("Invalid AnnotationFile argument. Please provide a gene annotation file, because you choose mode gene.")
    }
  } else if ( mode == "TE"){
    res <- res %>% rownames_to_column(var = "index") %>% mutate(gene_name = word(index, 1, sep = ":"))
    res_prot_ranked <- res[order(res$log2FoldChange, decreasing = T),c("gene_name", "log2FoldChange")]
    res_prot_ranked <- na.omit(res_prot_ranked)
  } else {
    stop("Invalid mode argument. Choose either gene or TE.")
  }
  res_prot_ranked = prepare_ranked_list(res_prot_ranked)
  if ( outfile != "" ){
    write.table(res_prot_ranked, file = outfile, sep = "\t", row.names = T, quote = F,col.names = FALSE)
  }
  
  return(res_prot_ranked)
}

waterfall_plot <- function (rnk,geneset,outjpeg,graph_title) {
  # read in file containing lists of genes for each pathway
  hallmark_pathway <- gmtPathways(geneset)
  fgsea_results <- fgsea(pathways = hallmark_pathway,
                  stats = rnk,
                  minSize = 15,
                  maxSize = 500
                  )
  fgsea_results %>% arrange (desc(NES)) %>% select (pathway, padj, NES) %>% head()

  fgsea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>% # removes 'HALLMARK_' from the pathway title 
    ggplot( aes(reorder(short_name,NES), NES)) +
      geom_bar(stat= "identity", aes(fill = padj<0.05))+
      coord_flip()+
      labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
      theme(axis.text.y = element_text(size = 7), 
            plot.title = element_text(hjust = 0.5))
  ggsave(outjpeg, plot = last_plot(), device = "jpeg", width = 8, height = 6, dpi = 300)

}
# wrapper for fgsea::plotEnrichment()
plot_enrichment <- function (geneset, pathway, ranked_list) {
  p = plotEnrichment(geneset[[pathway]], ranked_list)+labs (title = pathway) +
      theme_classic() +  
      theme(
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA,linewidth = 0.2), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()   
      )
  outjpeg = paste("output/DESeq2/GSEA/",pathway,".jpeg",sep="")
  jpeg(width = 200, height = 150, units = "mm", res = 300, file = outjpeg)
    print(p)
  dev.off()
}
parser <- ArgumentParser(description='GSEA analysis for TE and Gene')
parser$add_argument('-m', '--mode', type='character',default = "Gene",choices = c("Gene","TE"),
                    help='mode: "TE", "Gene"')
parser$add_argument('-g', '--gmt', type='character', required=TRUE,
                    help='path to gmt file')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
                    help='path to expression file')
parser$add_argument('-o', '--outdir', type='character', required=TRUE,
                    help='path to output dir')
parser$add_argument('-a', '--annotation', type='character', required=TRUE,
            help='annotation file, three column gene_id,gene_type,gene_name')
parser$add_argument('-t', '--graphTitle', type='character', default="",
            help='graph title')
        
args <- parser$parse_args()
mode = args$mode
gmt = args$gmt
matrix = args$matrix
outdir = args$outdir
annotation = args$annotation
graphTitle = args$graphTitle
print(paste("mode:",mode,sep=""))
print(paste("gmt path:",gmt,sep=""))
print(paste("input matrix path:",matrix,sep=""))
print(paste("outfile:",outdir,sep=""))
print(paste("annotation file:",annotation,sep=""))
print(paste("graph title: ",graphTitle,sep=""))

if ( !dir.exists(paste(outdir,"DESeq2/GSEA/",sep=""))){
    dir.create(paste(outdir,"DESeq2/GSEA/",sep=""),recursive = TRUE,showWarnings = FALSE)
}
if ( mode == "Gene"){
    df = read.csv(matrix,sep="\t",row.names=1,header=T);
    outRnk = paste(outdir,"DESeq2/GSEA/TEcount_Gene_GSEA.rnk",sep="")

    rnk = GSEA_prepare(df,mode = "Gene",outfile = outRnk,AnnotationFile=annotation)
    outEnrichJpeg = paste(outdir,"DESeq2/GSEA/TEcount_Gene_GSEA.jpeg",sep="")
    outEnrich = waterfall_plot(rnk,geneset = gmt,outjpeg = outEnrichJpeg, graph_title= graphTitle)
} else if ( mode == "TE"){
    stop("haven't develop method for this mode")
} else {
    stop("please using correct mode,Gene or TE")
}
