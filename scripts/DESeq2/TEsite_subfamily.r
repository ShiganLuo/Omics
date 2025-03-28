suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(argparse))
TEsiteSta = function(df,outjpeg,top = 10,mainTitle="TE subfamily",outfile = ""){
    df <- df %>% mutate(index = rownames(df)) %>%
    separate(index, into = c("site", "subfamily", "family", "class"), sep = ":")
    df_count <- df %>%
        count(subfamily, family)  %>% arrange(desc(n))
    if ( outfile != "") {
        write.csv(df_count,outfile)
    }
    df_count = df_count %>% head(top) # 统计 subfamily 数量，并保留 family 信息
    # 绘制柱状图
    p = ggplot(df_count, aes(x = reorder(subfamily, -n), y = n, fill = family)) +
    geom_bar(stat = "identity",width = 0.6) +  # 使用柱状图
    theme_minimal() +  
    labs(title = mainTitle, x = "Subfamily", y = "Number of loci", fill = "Family") +
       theme_minimal() +
    theme(
        panel.border = element_blank(),  # 去除整个面板边框
        axis.line = element_line(colour = "black"),  # 重新添加坐标轴线条
        panel.grid.major = element_blank(),  # 去除主要网格线
        panel.grid.minor = element_blank(),  # 去除次要网格线
        axis.ticks = element_line(colour = "black"),  # 设置坐标轴刻度线颜色
        axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
        plot.title = element_text(hjust = 0.5),  # 标题居中
        axis.title = element_text(colour = "black"),  # 设置坐标轴标题颜色
        axis.text.x = element_text(angle = 90, hjust = 1)
    )
    jpeg(width = 150, height = 150, units = "mm", res = 300, file = outjpeg)
        print(p)
    dev.off()

}
parser <- ArgumentParser(description='TE subfamily statistics for TElocal using DESeq2 results')
parser$add_argument('-u', '--up', type='character',required = TRUE,
                    help='path to upregulated TElocal DESeq2 result file')
parser$add_argument('-d', '--down', type='character',required = TRUE,
                    help='path to downregulated TElocal DESeq2 result file')
parser$add_argument('-o', '--outdir', type='character', required=TRUE,
                    help='path to output dir')
parser$add_argument('-t', '--graphTitle', type='character', default="",
                    help='graph title suffix')
args <- parser$parse_args()
up = args$up
down = args$down
outdir = args$outdir
graphTitle = args$graphTitle
print(paste("upregulated TElocal DESeq2 result file:",up,sep=""))
print(paste("downregulated TElocal DESeq2 result file:",down,sep=""))
print(paste("output dir:",outdir,sep=""))
print(paste("graph title: ",graphTitle,sep=""))
if ( !dir.exists(paste(outdir,"DESeq2/TE/",sep=""))){
    dir.create(paste(outdir,"DESeq2/TE/",sep=""),recursive = TRUE,showWarnings = FALSE)
}

outjpeg = paste(outdir,"DESeq2/TE/TElocal_TE_up_subfamily.jpeg",sep="")
title = paste("Upregulated TE subfamily after",graphTitle,sep="")
outfile = paste(outdir,"DESeq2/TE/TElocal_TE_up_subfamily.csv",sep="")
df = read.csv(up,sep = "\t",header = T,row.names=1)
TEsiteSta(df,outjpeg,mainTitle = title,outfile = outfile)

outjpeg = paste(outdir,"DESeq2/TE/TElocal_TE_down_subfamily.jpeg",sep="")
title = paste("Downregulated TE subfamily after",graphTitle,sep="")
outfile = paste(outdir,"DESeq2/TE/TElocal_TE_down_subfamily.csv",sep="")
df = read.csv(down,sep = "\t",header = T,row.names=1)
TEsiteSta(df,outjpeg,mainTitle = title,outfile = outfile)