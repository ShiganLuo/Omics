# 安装依赖包（比较耗时，确保网络状态良好）；
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# 安装所需富集分析R包和人基因注释包；
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
install.packages(c("aPEAR", "ggplot2"))
# 加载R包；
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(aPEAR)
library(ggplot2)

# 读取基因列表；
df <- read.csv("OS_diffgenes.csv")
# 提取目的基因列表；
gene_list <- df$GENE_SYMBOL
# 基因id转换；
id_map <- bitr(gene_list, 
               fromType = "SYMBOL", 
               toType = "ENTREZID", 
               OrgDb = org.Hs.eg.db)
entrez_ids <- id_map$ENTREZID

# 执行GO富集分析；
go_res <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# 导出GO结果到本地；
write.csv(go_res@result, "GO_enrichment_results.csv", row.names = FALSE)

# 执行KEGG富集分析；
kegg_res <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)
# 导出KEGG结果到本地；
write.csv(kegg_res@result, "KEGG_enrichment_results.csv", row.names = FALSE)

# 绘制GO富集条形图（Top15）；
barplot(go_res, 
        showCategory = 15, 
        x = "Count",
        color = "pvalue",
        label_format = 50,
        font.size = 11,
        title = "GO Pathway Enrichment") +
  theme(text = element_text(size = 12))

# 绘制GO富集气泡图（Top15）；
dotplot(go_res, 
        showCategory = 15, 
        color = "pvalue",
        label_format = 50,
        font.size = 11,
        title = "GO Pathway Enrichment") +
  theme(text = element_text(size = 12))


# 绘制KEGG富集条形图（Top15）；
barplot(kegg_res, 
        showCategory = 15, 
        x = "Count",
        color = "p.adjust",
        label_format = 45,
        font.size = 11,
        title = "KEGG Pathway Enrichment") +
  theme(text = element_text(size = 12))

# 绘制KEGG富集气泡图（Top15），尝试通过pvalue进行排序；
dotplot(kegg_res, 
        showCategory = 15, 
        color = "p.adjust",
        label_format = 50,
        font.size = 11,
        orderBy = "pvalue",
        decreasing = F,
        title = "KEGG Pathway Enrichment") +
  theme(text = element_text(size = 8))