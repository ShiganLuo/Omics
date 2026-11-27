#### 导入包 ----
library(ggplot2)

df <- readRDS("./results/function/GO/tc.rds")
head(df)
pathways <- c(
  "translational initiation",
  "mRNA catabolic process",
  "RNA catabolic process",
  "muscle system process",
  "platelet activation"
)
df <- df[df$Description %in% pathways,]

p <- ggplot(data = df, aes(
  x = -log10(pvalue),
  y = reorder(Description,-log10(pvalue))
)) +
  geom_bar(stat = 'identity', fill = "#E64B35FF") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 8),
    axis.title.x = element_text(color = "black", size = 8),
    axis.title.y = element_blank()
  )
p
ggsave(filename = "./results/function/GO/tc_up.png", plot = p,
       height = 3, width = 3)
