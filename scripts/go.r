suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
go_plot = function(gene_list,gofile,gojpeg,top=10,mainTitle = "GO Enrichment Analysis"){
    gene_info <- bitr(gene_list, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db,drop = TRUE)
    # Get a unique Entrez ID
    entrez_ids <- unique(gene_info$ENTREZID)

    ego <- enrichGO(gene = entrez_ids,
                    OrgDb = org.Mm.eg.db,
                    ont = "ALL",  # 选择 "ALL" 表示同时分析 BP、CC 和 MF
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    pvalueCutoff = 0.05,
                    readable = TRUE)
    df = ego@result
    write.table(ego@result, gofile, sep="\t", quote=FALSE, row.names=FALSE)
    # df = read.csv(gofile,sep = "\t",header = T)
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
if ( !dir.exists(paste(outdir,"/DESeq2/go/",sep=""))){
    dir.create(paste(outdir,"/DESeq2/go/",sep=""),recursive = TRUE,showWarnings = FALSE)
}
df = read.csv("output/DESeq2/upDown/TEcount_Gene_up.csv", header = TRUE, sep = "\t",row.names =1)
gene_list <- sapply(strsplit(rownames(df), "\\."), `[`, 1)
outfile = "output/DESeq2/go/TEcount_Gene_up_GO.txt"
go_plot(gene_list,gofile = outfile,gojpeg = "output/DESeq2/go/TEcount_Gene_up_GO.jpeg",mainTitle = "Upregulated genes after GCN2 knock-out")

df = read.csv("output/DESeq2/upDown/TEcount_Gene_down.csv",header = TRUE, sep = "\t",row.names =1)
gene_list <- sapply(strsplit(rownames(df), "\\."), `[`, 1)
outfile = "output/DESeq2/go/TEcount_Gene_down_GO.txt"
go_plot(gene_list,gofile = outfile,gojpeg = "output/DESeq2/go/TEcount_Gene_down_GO.jpeg",mainTitle = "Downregulated genes after GCN2 knock-out")



