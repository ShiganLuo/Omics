#!/usr/bin/env Rscript

## =========================
## Libraries
## =========================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(fgsea)
  library(argparse)
  library(tibble)
  library(stringr)
  library(data.table)
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
## Utility: adaptive plot height
## =========================
calc_plot_height <- function(n,
                             base = 4,
                             per_item = 0.18,
                             min_h = 4,
                             max_h = 20) {
  h <- base + n * per_item
  h <- max(min_h, min(h, max_h))
  return(h)
}

## =========================
## Prepare ranked list
## =========================
prepare_ranked_list <- function(ranked_df) {
  ### ranked_df: data.frame with gene_name and log2FoldChange columns 
  if (!all(c("gene_name", "log2FoldChange") %in% colnames(ranked_df))) {
    log_msg("ERROR",
            "ranked_df must contain gene_name and log2FoldChange",
            quit = TRUE)
  }

  if (any(duplicated(ranked_df$gene_name))) {
    log_msg("WARN",
            "Duplicated gene_name detected, averaging log2FoldChange")
    ranked_df <- ranked_df %>%
      group_by(gene_name) %>%
      summarise(log2FoldChange = mean(log2FoldChange, na.rm = TRUE),
                .groups = "drop")
  }

  ranked_df <- ranked_df %>%
    filter(!is.na(log2FoldChange)) %>%
    arrange(desc(log2FoldChange))

  ## break ties
  set.seed(123)
  ranked_df$log2FoldChange <- ranked_df$log2FoldChange +
    runif(nrow(ranked_df), -1e-5, 1e-5)

  tibble::deframe(ranked_df)
}

## =========================
## GSEA prepare
## =========================
GSEA_prepare <- function(res,
                         mode = "Gene",
                         AnnotationFile = "",
                         outfile = "") {

  log_msg("INFO", "Preparing ranked list, mode =", mode)

  if (mode == "Gene") {

    if (AnnotationFile == "") {
      log_msg("ERROR",
              "AnnotationFile is required for Gene mode",
              quit = TRUE)
    }

    anno <- fread(AnnotationFile)
    required_cols <- c("gene_id", "gene_type", "gene_name")
    if (!all(required_cols %in% colnames(anno))) {
      log_msg("ERROR",
              "Annotation file must contain:",
              paste(required_cols, collapse = ", "),
              quit = TRUE)
    }

    res <- res %>%
      rownames_to_column("gene_id") %>%
      left_join(anno, by = "gene_id") %>%
      filter(gene_type == "protein_coding")

    ranked_df <- res %>%
      select(gene_name, log2FoldChange)

  } else if (mode == "TE") {

    res <- res %>%
      rownames_to_column("index") %>%
      mutate(gene_name = word(index, 1, sep = ":"))

    ranked_df <- res %>%
      select(gene_name, log2FoldChange)

  } else {
    log_msg("ERROR", "Invalid mode:", mode, quit = TRUE)
  }

  ranked_vec <- prepare_ranked_list(ranked_df)

  if (outfile != "") {
    fwrite(
      data.table(
        gene = names(ranked_vec),
        score = ranked_vec
      ),
      file = outfile,
      sep = "\t"
    )
    log_msg("INFO", "Ranked list written to", outfile)
  }

  ranked_vec
}

gradient_bar_panel <- function(n) {

  df_grad <- data.frame(
    x = seq_len(n),
    y = 1,
    value = seq(1, -1, length.out = n)
  )

  ggplot(df_grad, aes(x, y, fill = value)) +
    geom_raster() +
    scale_fill_gradient2(
      low = "#4575B4",   # Negative
      mid = "white",
      high = "#D73027", # Positive
      midpoint = 0,
      guide = "none"
    ) +
    # annotate(
    #   "text",
    #   x = 1,
    #   y = 1.3,
    #   label = "Negative",
    #   hjust = 0,
    #   size = 3
    # ) +
    # annotate(
    #   "text",
    #   x = n,
    #   y = 1.3,
    #   label = "Positive",
    #   hjust = 1,
    #   size = 3
    # ) +
    theme_void() +
    theme(
      plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
    )
}

## =========================
## plot enrichment in gseapy style
## =========================
plot_enrichment_gseapy_style <- function(pathway,
                                         geneset,
                                         ranked_list,
                                         fgsea_row,
                                         outjpeg,
                                         width = 7,
                                         height = 5) {

  library(ggplot2)
  library(patchwork)

  ## ---- geneset ----
  pathway_genes <- geneset[[pathway]]
  if (is.null(pathway_genes)) {
    log_msg("WARN", "Pathway not found:", pathway)
    return(invisible(NULL))
  }

  genes <- names(ranked_list)

  hit <- genes %in% pathway_genes
  hit_idx <- which(hit)

  if (sum(hit) < 2) {
    log_msg("WARN", "Too few overlapping genes for:", pathway)
    return(invisible(NULL))
  }

  ## ---- running ES (gseapy definition) ----
  stats <- ranked_list
  Nh <- sum(hit)
  Nm <- length(stats) - Nh

  Phit  <- sum(abs(stats[hit]))
  Pmiss <- Nm

  running_es <- cumsum(
    ifelse(
      hit,
      abs(stats) / Phit,
      -1 / Pmiss
    )
  )

  df_es <- data.frame(
    position = seq_along(running_es),
    ES = running_es
  )

  ## ES peak (for annotation)
  peak_pos <- which.max(abs(running_es))
  peak_es  <- running_es[peak_pos]

  ## ---- hit ticks ----
  df_hits <- data.frame(
    position = hit_idx,
    y = 1
  )

  ## ---- rank metric ----
  df_rank <- data.frame(
    position = seq_along(stats),
    score = stats
  )

  nes  <- round(fgsea_row$NES, 3)
  padj <- signif(fgsea_row$padj, 3)

  ## ================= Panel 1 =================
  p1 <- ggplot(df_es, aes(position, ES)) +
    geom_line(color = "#D62728", linewidth = 1) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               linewidth = 0.4,
               color = "grey40") +
    geom_vline(xintercept = peak_pos,
               linetype = "dotted",
               linewidth = 0.4,
               color = "grey50") +
    labs(
      title = pathway,
      subtitle = paste0("NES = ", nes, ", adj.P = ", padj),
      y = "Enrichment score",
      x = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11,hjust = 0.5),
      plot.subtitle = element_text(size = 9,hjust = 0),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  ## ================= Panel 2 =================
  p2 <- ggplot(df_hits, aes(position, y)) +
    geom_linerange(aes(ymin = 0, ymax = 1),
                   linewidth = 0.4,
                   color = "black") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void()

  p_grad <- gradient_bar_panel(length(ranked_list))
  ## ================= Panel 3 =================
  # score represent log2FC in there,decided by ranked_list
  # scores <- as.numeric(ranked_list)
  # zero_cross_idx <- which(scores < 0)[1]
  # if (is.na(zero_cross_idx)) zero_cross_idx <- length(scores)
  p3 <- ggplot(df_rank, aes(position, score, fill = score)) +
    geom_area() +
    scale_fill_gradient2(
      low = "#4575B4",   # Negative（右侧）
      mid = "grey85",
      high = "#D73027", # Positive（左侧）
      midpoint = 0,
      guide = "none"
    ) +
    geom_hline(yintercept = 0,
              linewidth = 0.3,
              color = "grey40") +
    labs(
      y = "Rank metric",
      x = "Rank in ordered dataset"
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    annotate(
      "text",
      x = 1,
      y = max(df_rank$score),
      label = "Positive",
      hjust = 0,
      vjust = 0.5,
      size = 3,
      color = "red"
    ) +
    annotate(
      "text",
      x = nrow(df_rank),
      y = min(df_rank$score),
      label = "Negative",
      hjust = 1,
      vjust = 0.2,
      size = 3,
      color = "blue"
    )

  ## ---- combine ----
  p <- p1 / p2 / p_grad / p3 +
    plot_layout(
      heights = c(2.2, 0.4, 0.25, 1)
    )

  ggsave(
    outjpeg,
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )

  invisible(p)
}




## =========================
## Waterfall plot (adaptive)
## =========================
waterfall_plot <- function(rnk,
                           geneset,
                           outjpeg,
                           graph_title,
                           outfile) {

  log_msg("INFO", "Running fgsea")

  if (!file.exists(outfile)) {

    pathways <- gmtPathways(geneset)

    fgsea_res <- fgsea(
      pathways = pathways,
      stats = rnk,
      minSize = 15,
      maxSize = 500,
      nPermSimple = 10000
    )

    fwrite(fgsea_res, outfile)

  } else {
    fgsea_res <- fread(outfile)
  }

  fgsea_res <- fgsea_res %>%
    filter(!is.na(pval)) %>%
    arrange(desc(NES)) %>%
    select(pathway, padj, NES) %>%
    mutate(short_name = str_remove(pathway, "^HALLMARK_"))

  n_pathway <- nrow(fgsea_res)
  height <- calc_plot_height(n_pathway)

  log_msg("INFO",
          "Plotting", n_pathway, "pathways;",
          "height =", height)

  p <- ggplot(fgsea_res,
              aes(reorder(short_name, NES), NES)) +
    geom_col(aes(fill = padj < 0.05)) +
    coord_flip() +
    labs(
      x = "Hallmark Pathway",
      y = "Normalized Enrichment Score",
      title = graph_title
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 7),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )

  ggsave(
    outjpeg,
    plot = p,
    width = 8,
    height = height,
    dpi = 300
  )
  ## ===== gseapy-style enrichment plots =====
  top_pathways <- fgsea_res %>%
    arrange(desc(NES)) %>%
    slice_head(n = 5)

  geneset_list <- gmtPathways(geneset)

  enrich_dir <- file.path(dirname(outjpeg), "enrichment_plots")
  dir.create(enrich_dir, showWarnings = FALSE)

  for (i in seq_len(nrow(top_pathways))) {
    pw <- top_pathways$pathway[i]

    outp <- file.path(
      enrich_dir,
      paste0(pw, "_enrichment.jpeg")
    )

    plot_enrichment_gseapy_style(
      pathway = pw,
      geneset = geneset_list,
      ranked_list = rnk,
      fgsea_row = top_pathways[i, ],
      outjpeg = outp
    )
  }
  invisible(fgsea_res)
}

## =========================
## plotEnrichment wrapper
## =========================
plot_enrichment <- function(geneset,
                            pathway,
                            ranked_list,
                            outdir) {

  p <- plotEnrichment(geneset[[pathway]], ranked_list) +
    labs(title = pathway) +
    theme_classic() +
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.2
      )
    )

  outjpeg <- file.path(outdir, paste0(pathway, ".jpeg"))

  ggsave(
    outjpeg,
    plot = p,
    width = 8,
    height = 4,
    dpi = 300
  )
}

## =========================
## CLI
## =========================
parser <- ArgumentParser(description = "GSEA analysis for TE and Gene")

parser$add_argument("-m", "--mode",default = "Gene",choices = c("Gene", "TE"),type="character", help = "GSEA mode: Gene or TE")
parser$add_argument("-g", "--gmt", required = TRUE, type="character", help = "Input GMT file for gene sets")
parser$add_argument("-i", "--matrix",required = TRUE,type="character", help = "Input matrix file (DESeq2 result)")
parser$add_argument("-o", "--outdir",required = TRUE, type="character", help = "Output directory")
parser$add_argument("-a", "--annotation", required = TRUE, type="character",
                    help = "Gene annotation file for Gene mode  (including gene_id, gene_type, gene_name columns)")
parser$add_argument("-t", "--graphTitle",type = "character", default = "", help = "Title for the GSEA graph")
parser$add_argument("-r", "--rewrite",action = "store_true", default = FALSE, help = "Whether to rewrite existing files")
args <- parser$parse_args()

## =========================
## Main
## =========================
log_msg("INFO", "Mode:", args$mode)
log_msg("INFO", "GMT:", args$gmt)
log_msg("INFO", "Matrix:", args$matrix)
log_msg("INFO", "Outdir:", args$outdir)

gsea_dir <- file.path(args$outdir, "GSEA")
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)

if (args$mode == "Gene") {

  df <- fread(args$matrix, data.table = FALSE)
  rownames(df) <- df[[1]]
  df[[1]] <- NULL

  out_rnk <- file.path(gsea_dir, "TEcount_Gene_GSEA.rnk")
  out_fgsea <- file.path(gsea_dir, "TEcount_Gene_GSEA.csv")
  out_jpeg <- file.path(gsea_dir, "TEcount_Gene_GSEA.jpeg")

  if (!file.exists(out_rnk) || args$rewrite) {
    rnk <- GSEA_prepare(
      df,
      mode = "Gene",
      AnnotationFile = args$annotation,
      outfile = out_rnk
    )
  } else {
    rnk_df <- fread(out_rnk)
    rnk <- setNames(rnk_df[[2]], rnk_df[[1]])
  }

  waterfall_plot(
    rnk = rnk,
    geneset = args$gmt,
    outjpeg = out_jpeg,
    graph_title = args$graphTitle,
    outfile = out_fgsea
  )

} else {
  log_msg("ERROR",
          "TE mode not implemented yet",
          quit = TRUE)
}
