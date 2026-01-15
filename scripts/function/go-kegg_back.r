#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
})

# ============================================================
# 1. 从 enrichment 表中定义 Up / Down 基因
# ============================================================
define_up_down_genes <- function(
  infile,
  gene_col,
  value_col,
  p_col,
  up_cutoff = 1,
  p_cutoff = 0.05
) {
  df <- read.csv(infile, sep = "\t", header = TRUE)

  stopifnot(all(c(gene_col, value_col, p_col) %in% colnames(df)))

  up <- df %>%
    filter(.data[[value_col]] >= up_cutoff,
           .data[[p_col]] <= p_cutoff) %>%
    pull(.data[[gene_col]]) %>%
    unique()

  down <- df %>%
    filter(.data[[value_col]] <= 1 / up_cutoff,
           .data[[p_col]] <= p_cutoff) %>%
    pull(.data[[gene_col]]) %>%
    unique()

  list(up = up, down = down)
}

# ============================================================
# 2. GO / KEGG enrichment（真正的核心模块）
# ============================================================
run_go_kegg <- function(
  genes,
  species = c("human", "mouse"),
  type = c("go", "kegg")
) {
  species <- match.arg(species)
  type <- match.arg(type)
  genes_clean <- sub("\\..*$", "", genes)
  if (length(genes) == 0) return(data.frame())

  if (species == "human") {
    OrgDb <- org.Hs.eg.db
    kegg_org <- "hsa"
  } else {
    OrgDb <- org.Mm.eg.db
    kegg_org <- "mmu"
  }

  gene_df <- bitr(
    genes_clean,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = OrgDb,
    drop = TRUE
  )

  entrez <- unique(gene_df$ENTREZID)
  if (length(entrez) == 0) return(data.frame())

  if (type == "go") {
    res <- enrichGO(
      gene = entrez,
      OrgDb = OrgDb,
      ont = "BP",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      readable = TRUE
    )
  } else {
    res <- enrichKEGG(
      gene = entrez,
      organism = kegg_org,
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  }

  as.data.frame(res@result)
}

# ============================================================
# 3. 背靠背 GO / KEGG 柱状图
# ============================================================
plot_back_to_back <- function(
  up_df,
  down_df,
  top = 10,
  title = "",
  outfile
) {
  up <- up_df %>%
    arrange(pvalue) %>%
    slice_head(n = top) %>%
    mutate(
      Group = "Up",
      value = -log10(pvalue)
    )

  down <- down_df %>%
    arrange(pvalue) %>%
    slice_head(n = top) %>%
    mutate(
      Group = "Down",
      value = -(-log10(pvalue))
    )

  df <- bind_rows(up, down)

  p <- ggplot(
    df,
    aes(x = reorder(Description, value),
        y = value,
        fill = Group)
  ) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c(Up = "#D73027", Down = "#4575B4")) +
    labs(
      title = title,
      x = NULL,
      y = expression(-log[10](pvalue))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line.x = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5)
    )

  ggsave(outfile, p, width = 8, height = 6)
}

# ============================================================
# 4. 总控流程（围绕 GO / KEGG）
# ============================================================
run_pipeline <- function(
  infile,
  outdir,
  species = "human",
  gene_col = "gene_id",
  value_col = "total_enrichment",
  p_col = "total_pvalue",
  up_cutoff = 1,
  p_cutoff = 0.05,
  top = 10
) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  message("Defining Up / Down genes ...")
  genes <- define_up_down_genes(
    infile,
    gene_col,
    value_col,
    p_col,
    up_cutoff,
    p_cutoff
  )

  writeLines(genes$up, file.path(outdir, "up_genes.txt"))
  writeLines(genes$down, file.path(outdir, "down_genes.txt"))

  for (type in c("go", "kegg")) {
    message("Running ", toupper(type), " enrichment ...")

    up_res <- run_go_kegg(genes$up, species, type)
    down_res <- run_go_kegg(genes$down, species, type)

    write.csv(up_res, file.path(outdir, paste0(type, "_up.csv")), row.names = FALSE)
    write.csv(down_res, file.path(outdir, paste0(type, "_down.csv")), row.names = FALSE)

    plot_back_to_back(
      up_res,
      down_res,
      top = top,
      title = paste("LINE1", toupper(type), "enrichment"),
      outfile = file.path(outdir, paste0(type, "_back_to_back.png"))
    )
  }
}

# ============================================================
# 5. 示例：你的 LINE1 enrichment
# ============================================================
run_pipeline(
  infile = "/disk5/luosg/DIPseq20251215/output/LINE1/LINE1_enrichment_per_gene.tsv",
  outdir = "/disk5/luosg/DIPseq20251215/output/LINE1/GO_KEGG",
  species = "human",
  gene_col = "gene_id",
  value_col = "total_enrichment",
  p_col = "total_pvalue",
  up_cutoff = 1,
  p_cutoff = 0.05,
  top = 10
)
