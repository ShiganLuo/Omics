suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
kegg_plot = function(gene_list,keggfile,keggjpeg,top=10,mainTitle = "GO Enrichment Analysis"){
    gene_info <- bitr(gene_list, fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
    entrez_ids <- unique(gene_info$ENTREZID)
    ekegg <- enrichKEGG(
        gene = entrez_ids,     
        organism = "mmu",     
        keyType = "kegg",      
        pvalueCutoff = 0.05,  
        qvalueCutoff = 0.05,   
        pAdjustMethod = "BH"
    )
    df = ekegg@result
    write.table(ekegg@result, keggfile, sep="\t", quote=FALSE, row.names=FALSE)
    # df = read.csv(keggfile,sep = "\t",header = T)
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
if ( !dir.exists(paste(outdir,"/DESeq2/kegg/",sep=""))){
    dir.create(paste(outdir,"/DESeq2/kegg/",sep=""),recursive = TRUE,showWarnings = FALSE)
}
df = read.csv("output/DESeq2/upDown/TEcount_Gene_up.csv", header = TRUE, sep = "\t",row.names =1)
gene_list <- sapply(strsplit(rownames(df), "\\."), `[`, 1)
outfile = "output/DESeq2/kegg/TEcount_Gene_up_KEGG.txt"
outjpeg = "output/DESeq2/kegg/TEcount_Gene_up_KEGG.jpeg"
kegg_plot(gene_list,keggfile = outfile,keggjpeg = outjpeg,mainTitle = "Upregulated genes after GCN2 knock-out")

df = read.csv("output/DESeq2/upDown/TEcount_Gene_down.csv",header = TRUE, sep = "\t",row.names =1)
gene_list <- sapply(strsplit(rownames(df), "\\."), `[`, 1)
outfile = "output/DESeq2/kegg/TEcount_Gene_down_KEGG.txt"
outjpeg = "output/DESeq2/kegg/TEcount_Gene_down_KEGG.jpeg"
kegg_plot(gene_list,keggfile = outfile,keggjpeg = outjpeg,mainTitle = "Downregulated genes after GCN2 knock-out")
