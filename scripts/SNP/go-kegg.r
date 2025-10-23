suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(org.Hs.eg.db))
go_plot = function(gene_list,gofile,gojpeg,top=10,mainTitle = "GO Enrichment Analysis",species = "mouse"){
    if (!file.exists(gofile)) {
        if ( species == "mouse") {
            gene_info <- bitr(gene_list, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db,drop = TRUE)
            # Get a unique Entrez ID
            entrez_ids <- unique(gene_info$ENTREZID)
            print(paste("entrez_ids:",length(entrez_ids)))
            ego <- enrichGO(gene = entrez_ids,
                        OrgDb = org.Mm.eg.db,
                        ont = "ALL",  # 选择 "ALL" 表示同时分析 BP、CC 和 MF
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05,
                        pvalueCutoff = 0.05,
                        readable = TRUE)
            df = ego@result
        } else if ( species == "human") {
            gene_info <- bitr(gene_list, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db,drop = TRUE)
            # Get a unique Entrez ID
            entrez_ids <- unique(gene_info$ENTREZID)
            print(paste("entrez_ids:",length(entrez_ids)))
            ego <- enrichGO(gene = entrez_ids,
            OrgDb = org.Hs.eg.db,
            ont = "ALL",  # 选择 "ALL" 表示同时分析 BP、CC 和 MF
            pAdjustMethod = "BH",
            qvalueCutoff = 0.05,
            pvalueCutoff = 0.05,
            readable = TRUE)
            df = ego@result
        } else {
            stop("species not supported")
        }
        write.table(ego@result, gofile, sep="\t", quote=FALSE, row.names=FALSE)
    } else {
        df = read.csv(gofile,sep = "\t",header = T)
    }

    
    df = df %>% arrange(pvalue, Count) %>% slice_head(n = top) 
    p = ggplot(df, aes(x = reorder(df[,"Description"],  -log10(df[,"pvalue"])), y = -log10(df[,"pvalue"]))) +
        geom_point(aes(color = -log10(df[,"pvalue"]), size = Count )) +
        scale_color_gradient(low = "#E2A3A3", high = "#F54E4E")+
        coord_flip() +  # 翻转坐标轴，使GO_Term在纵轴上显示
        labs(title = mainTitle,
            x = "",
            y = expression(-log[10]("pvalue")),
            color = expression(-log[10]("pvalue")),
            size = "Count") +
        guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))+
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
    jpeg(width = 200, height = 150, units = "mm", res = 300, file = gojpeg)
        print(p)
    dev.off()
}
# if ( !dir.exists(paste(outdir,"/DESeq2/go/",sep=""))){
#     dir.create(paste(outdir,"/DESeq2/go/",sep=""),recursive = TRUE,showWarnings = FALSE)
# }
# df = read.csv("output/DESeq2/upDown/TEcount_Gene_up.csv", header = TRUE, sep = "\t",row.names =1)
# gene_list <- sapply(strsplit(rownames(df), "\\."), `[`, 1)
# outfile = "output/DESeq2/go/TEcount_Gene_up_GO.txt"
# go_plot(gene_list,gofile = outfile,gojpeg = "output/DESeq2/go/TEcount_Gene_up_GO.jpeg",mainTitle = "Upregulated genes after GCN2 knock-out")

# df = read.csv("output/DESeq2/upDown/TEcount_Gene_down.csv",header = TRUE, sep = "\t",row.names =1)
# gene_list <- sapply(strsplit(rownames(df), "\\."), `[`, 1)
# outfile = "output/DESeq2/go/TEcount_Gene_down_GO.txt"
# go_plot(gene_list,gofile = outfile,gojpeg = "output/DESeq2/go/TEcount_Gene_down_GO.jpeg",mainTitle = "Downregulated genes after GCN2 knock-out")
kegg_plot = function(gene_list,keggfile,keggjpeg,top=10,mainTitle = "kegg Enrichment Analysis",species = "mouse"){
    if (!file.exists(keggfile)){
        if ( species == "mouse") {
            gene_info <- bitr(gene_list, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
            entrez_ids <- unique(gene_info$ENTREZID)
            print(paste("entrez_ids:",length(entrez_ids)))
            ekegg <- enrichKEGG(
                gene = entrez_ids,     
                organism = "mmu",     
                keyType = "kegg",      
                pvalueCutoff = 0.05,  
                qvalueCutoff = 0.05,   
                pAdjustMethod = "BH"
            )
            df = ekegg@result
        } else if ( species == "human") {
            gene_info <- bitr(gene_list, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
            entrez_ids <- unique(gene_info$ENTREZID)
            print(paste("entrez_ids:",length(entrez_ids)))
            ekegg <- enrichKEGG(
                gene = entrez_ids,     
                organism = "hsa",     
                keyType = "kegg",      
                pvalueCutoff = 0.05,  
                qvalueCutoff = 0.05,   
                pAdjustMethod = "BH"
            )
            df = ekegg@result
        } else {
            stop("species not supported")
        }
        write.table(ekegg@result, keggfile, sep="\t", quote=FALSE, row.names=FALSE)
    } else {
        df = read.csv(keggfile,sep = "\t",header = T)
    }

    
    df = df %>% arrange(pvalue, Count) %>% slice_head(n = top) 
    p = ggplot(df, aes(x = reorder(df[,"Description"],  -log10(df[,"pvalue"])), y = -log10(df[,"pvalue"]))) +
        geom_point(aes(color = -log10(df[,"pvalue"]), size = Count )) +
        scale_color_gradient(low = "#E2A3A3", high = "#F54E4E")+
        coord_flip() +  # 翻转坐标轴，使GO_Term在纵轴上显示
        labs(title = mainTitle,
            x = "",
            y = expression(-log[10]("pvalue")),
            color = expression(-log[10]("pvalue")),
            size = "Count") +
        guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))+
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
    jpeg(width = 200, height = 150, units = "mm", res = 300, file = keggjpeg)
        print(p)
    dev.off()

}

parser <- ArgumentParser(description='go and kegg enrichment analysis for the output of alteredGene')
parser$add_argument('-m', '--mode', type='character', required=TRUE,nargs="+",choices=c("go","kegg"),
                    help='mode: go or kegg ')
parser$add_argument('-i', '--infile', type='character', required=TRUE,
                    help='up: the outfile of alteredGene')
parser$add_argument('-o', '--outdir', type='character', required=TRUE,
                    help='outdir: item output dir,such as output/GSE166216/')
parser$add_argument('-t', '--title', type='character',default="",
                    help='title: the main title suffix of go jpeg;such Upregulated/Downregulated genes after [title]')
parser$add_argument('-s', '--species', type='character', default="mouse",choices=c("human","mouse"),
                    help='species: the species of the input file; default is mouse')
args <- parser$parse_args()
mode = args$mode
infile = args$infile
outdir = args$outdir
title = args$title
species = args$species

# cat("mode:",mode,"\n")
print(paste("the output file of alteredGene: ",infile,sep=""))
# print(paste("title suffix: ",title,sep=""))


df = read.csv(infile,header = TRUE,sep="\t")

df = df %>% arrange(desc(difference))
write.csv(df,"c.txt")
df_up = df[which(df$difference > 0),] %>%
    filter(difference >= quantile(difference, probs = 0.95)) %>% 
    arrange(desc(difference))

df_down = df[which(df$difference < 0),] %>% 
    filter(difference <= quantile(difference, probs = 0.05)) %>% 
    arrange(difference) 

gene_list_up <- sapply(strsplit(df_up[,'gene_id'], "\\."), `[`, 1)
gene_list_down <- sapply(strsplit(df_down[,'gene_id'], "\\."), `[`, 1)
if ( "go" %in% mode) {
    if ( !dir.exists(paste(outdir,"SNP/go/",sep=""))){
        dir.create(paste(outdir,"SNP/go/",sep=""),recursive = TRUE,showWarnings = FALSE)
    }
    # print(upBaseName)
    if (title != "") {
        uptitle = paste("SNV upregulated genes after ",title,sep="")
        downtitle = paste("SNV downregulated genes after ",title,sep="")
    } else {
        uptitle = "SNV upregulated genes"
        downtitle = "SNV downregulated genes"
    }
    outfile = paste(outdir,"SNP/go/up_GO.csv",sep="")
    outjpeg = paste(outdir,"SNP/go/up_GO.jpeg",sep="")
    go_plot(gene_list_up,gofile = outfile,gojpeg = outjpeg,mainTitle = uptitle,species = species)

    outfile = paste(outdir,"SNP/go/down_GO.csv",sep="")
    outjpeg = paste(outdir,"SNP/go/down_GO.jpeg",sep="")
    go_plot(gene_list_down,gofile = outfile,gojpeg = outjpeg,mainTitle = downtitle,species = species)
}

if ( "kegg" %in% mode) {
    if ( !dir.exists(paste(outdir,"SNP/kegg/",sep=""))){
        dir.create(paste(outdir,"SNP/kegg/",sep=""),recursive = TRUE,showWarnings = FALSE)
    }
    if (title != "") {
        uptitle = paste("SNV upregulated genes after ",title,sep="")
        downtitle = paste("SNV downregulated genes after ",title,sep="")
    } else {
        uptitle = "SNV upregulated genes"
        downtitle = "SNV downregulated genes"
    }
    outfile = paste(outdir,"SNP/kegg/up_kegg.csv",sep="")
    outjpeg = paste(outdir,"SNP/kegg/up_kegg.jpeg",sep="")
    kegg_plot(gene_list_up,keggfile = outfile,keggjpeg = outjpeg,mainTitle = uptitle,species = species)

    outfile = paste(outdir,"SNP/kegg/down_kegg.csv",sep="")
    outjpeg = paste(outdir,"SNP/kegg/down_kegg.jpeg",sep="")
    kegg_plot(gene_list_down,keggfile = outfile,keggjpeg = outjpeg,mainTitle = downtitle,species = species)
}







