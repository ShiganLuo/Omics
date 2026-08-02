#!/usr/bin/env Rscript

## =========================
## GO and KEGG enrichment analysis
## Accepts DESeq2 result TSV with a gene_name column.
## =========================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(argparse)
})

## =========================
## Unified logger
## =========================
log_msg <- function(level = c("INFO","WARN","ERROR"), ..., quit = FALSE) {
  level <- match.arg(level)
  msg <- paste(...)
  prefix <- paste0(
    "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
    "[", level, "] "
  )
  message(prefix, msg)
  if (quit) stop(msg, call. = FALSE)
}

## =========================
## Define up/down regulated genes
## =========================
define_up_down_genes <- function(
  infile,
  gene_col,
  value_col,
  p_col,
  lfc_cut = 1,
  p_cut = 0.05
) {
  df <- read.csv(infile, sep = "\t", header = TRUE, check.names = FALSE)

  if (!all(c(gene_col, value_col, p_col) %in% colnames(df))) {
    log_msg("ERROR",
            "Input file must contain columns:",
            paste(c(gene_col, value_col, p_col), collapse = ", "),
            quit = TRUE)
  }

  up <- df %>%
    filter(!is.na(.data[[p_col]]),
           .data[[value_col]] >= lfc_cut,
           .data[[p_col]] <= p_cut) %>%
    pull(.data[[gene_col]]) %>%
    unique()

  down <- df %>%
    filter(!is.na(.data[[p_col]]),
           .data[[value_col]] <= -lfc_cut,
           .data[[p_col]] <= p_cut) %>%
    pull(.data[[gene_col]]) %>%
    unique()

  list(up = up, down = down)
}

## =========================
## Run GO or KEGG enrichment
## =========================
run_go_kegg <- function(
  genes,
  species = c("human", "mouse"),
  type = c("go", "kegg")
) {
  species <- match.arg(species)
  type <- match.arg(type)
  if (length(genes) == 0) return(data.frame())

  if (species == "human") {
    OrgDb <- org.Hs.eg.db
    kegg_org <- "hsa"
  } else {
    OrgDb <- org.Mm.eg.db
    kegg_org <- "mmu"
  }

  gene_df <- bitr(
    genes,
    fromType = "SYMBOL",
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

## =========================
## Back-to-back bar plot
## =========================
clean_pathway_description <- function(description) {
  description <- as.character(description)
  description[is.na(description) | !nzchar(trimws(description))] <- "Unnamed pathway"
  sub(
    " - (Mus musculus \\(house mouse\\)|Homo sapiens \\(human\\))$",
    "",
    trimws(description),
    perl = TRUE
  )
}

wrap_plot_labels <- function(labels, width = 42) {
  newline <- intToUtf8(10)
  vapply(
    labels,
    function(label) {
      # Remove literal backslash-n sequences before adding real line breaks.
      label <- gsub("\\\\n", " ", label, fixed = FALSE)
      paste(strwrap(label, width = width), collapse = newline)
    },
    character(1)
  )
}

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

  if (nrow(df) == 0) {
    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No enriched terms", size = 6) +
      labs(title = title, x = NULL, y = NULL) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(outfile, p, width = 10, height = 4, dpi = 300, bg = "white")
    return(invisible(NULL))
  }

  # KEGG appends the organism to every Description. Keep it out of each label.
  df$Description <- clean_pathway_description(df$Description)
  df$plot_label <- wrap_plot_labels(df$Description)
  n_labels <- length(unique(df$plot_label))
  plot_height <- max(6, min(18, 2.5 + 0.32 * n_labels))

  p <- ggplot(
    df,
    aes(x = reorder(plot_label, value),
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
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.text.y = element_text(size = 10, lineheight = 0.9),
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(10, 24, 10, 10)
    )

  # Keep enough width for the numeric panel after allocating space to labels.
  ggsave(outfile, p, width = 12, height = plot_height, dpi = 300, bg = "white")
}

## =========================
## Pipeline
## =========================
run_pipeline <- function(
  infile,
  outdir,
  species = "human",
  gene_col = "gene_name",
  value_col = "log2FoldChange",
  p_col = "padj",
  lfc_cut = 1,
  p_cut = 0.05,
  top = 10
) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  log_msg("INFO", "Defining Up / Down genes ...")
  genes <- define_up_down_genes(
    infile,
    gene_col,
    value_col,
    p_col,
    lfc_cut,
    p_cut
  )

  log_msg("INFO", "Up genes:", length(genes$up),
          "  Down genes:", length(genes$down))

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
      title = paste(toupper(type), "enrichment"),
      outfile = file.path(outdir, paste0(type, "_back_to_back.png"))
    )
  }
}

## =========================
## CLI
## =========================
parser <- ArgumentParser(description = "GO and KEGG enrichment analysis for DESeq2 results")

parser$add_argument("-i", "--input", required = TRUE, type = "character",
                    help = "Input DESeq2 result TSV with gene_name column")
parser$add_argument("-o", "--outdir", required = TRUE, type = "character",
                    help = "Output directory")
parser$add_argument("-s", "--species", default = "mouse",
                    choices = c("human", "mouse"), type = "character",
                    help = "Species: human or mouse (default: mouse)")
parser$add_argument("--gene-col", default = "gene_name", type = "character",
                    help = "Column name for gene names (default: gene_name)")
parser$add_argument("--value-col", default = "log2FoldChange", type = "character",
                    help = "Column name for log2 fold change (default: log2FoldChange)")
parser$add_argument("--p-col", default = "padj", type = "character",
                    help = "Column name for adjusted p-value (default: padj)")
parser$add_argument("--lfc-cut", type = "double", default = 1,
                    help = "Absolute log2 fold change cutoff (default: 1)")
parser$add_argument("--p-cut", type = "double", default = 0.05,
                    help = "Adjusted p-value cutoff (default: 0.05)")
parser$add_argument("--top", type = "integer", default = 10,
                    help = "Top N pathways to plot (default: 10)")
args <- parser$parse_args()

## =========================
## Main
## =========================
log_msg("INFO", "Input:", args$input)
log_msg("INFO", "Outdir:", args$outdir)
log_msg("INFO", "Species:", args$species)

run_pipeline(
  infile = args$input,
  outdir = args$outdir,
  species = args$species,
  gene_col = args$gene_col,
  value_col = args$value_col,
  p_col = args$p_col,
  lfc_cut = args$lfc_cut,
  p_cut = args$p_cut,
  top = args$top
)

log_msg("INFO", "GO/KEGG analysis completed.")
