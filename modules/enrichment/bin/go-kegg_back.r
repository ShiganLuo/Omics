#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
})

## =========================
## Multi-format figure output (modules.md #9)
## =========================
## ggsave derives the device from the file extension, so writing every
## requested format is a plain loop over "{out_base}.{fmt}". pgf/raw/rgba
## are matplotlib-only and fail fast here (R ggsave device set).
VALID_FORMATS <- c("png", "pdf", "svg", "ps", "eps", "tif", "tiff",
                   "jpg", "jpeg")

## Fail fast on unsupported formats, listing what ggsave can write.
validate_formats <- function(formats) {
  bad <- setdiff(formats, VALID_FORMATS)
  if (length(bad) > 0) {
    stop(sprintf("Unsupported image format(s): %s; valid: %s",
                 paste(bad, collapse = ", "),
                 paste(VALID_FORMATS, collapse = ", ")))
  }
  invisible(formats)
}

## Save one ggplot per requested format as "{out_base}.{fmt}".
save_ggsave_formats <- function(plot, out_base, formats, width, height,
                                dpi = 300, ...) {
  validate_formats(formats)
  for (fmt in formats) {
    out <- sprintf("%s.%s", out_base, fmt)
    ggsave(out, plot = plot, width = width, height = height, dpi = dpi, ...)
  }
  invisible(out_base)
}

define_up_down_genes <- function(
  infile,
  gene_col,
  value_col,
  p_col,
  up_cutoff = 1,
  p_cutoff = 0.05
) {
  df <- read.csv(infile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  missing_cols <- setdiff(c(gene_col, value_col, p_col), colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "missing columns: ", paste(missing_cols, collapse = ", "),
      "; available: ", paste(colnames(df), collapse = ", ")
    )
  }

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


  # 判断输入是否已经是ENTREZID
  if (all(grepl("^[0-9]+$", genes_clean))) {
    entrez <- unique(genes_clean)
    message("Input genes look like ENTREZID, skipping bitr mapping.")
  } else {
    gene_df <- bitr(
      genes_clean,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = OrgDb,
      drop = TRUE
    ) # from clusterProfiler, need organism-specific OrgDb loaded
    entrez <- unique(gene_df$ENTREZID)
    if (length(entrez) == 0) return(data.frame())
  }

  if (type == "go") {
    res <- enrichGO(
      gene = entrez,
      OrgDb = OrgDb,
      ont = "BP",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      readable = TRUE
    ) # don't rely on internet for analysis, use OrgDb for GO annotation
  } else {
   res <- enrichKEGG(
      gene = entrez,
      organism = kegg_org,
      keyType = "kegg",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      use_internal_data = FALSE
    ) # need internet for KEGG analysis, but can specify organism to avoid mapping issues
    if (!is.null(res)) {
      message("enrichKEGG result row count: ", nrow(as.data.frame(res@result)))
    }
  }
  if (is.null(res)) {
    message("No enrichment results returned for ", toupper(type), ".")
    return(data.frame())
  }
  as.data.frame(res@result)
}


plot_back_to_back <- function(
  up_df,
  down_df,
  top = 10,
  title = "",
  out_base,
  formats = "png"
) {
  required_cols <- c("pvalue", "Description")

  prepare_df <- function(df, group_label) {
    if (nrow(df) == 0) {
      return(data.frame())
    }
    if (!all(required_cols %in% colnames(df))) {
      message(
        group_label, " results missing columns (",
        paste(required_cols, collapse = ", "), ") for: ", title,
        "; available: ", paste(colnames(df), collapse = ", ")
      )
      return(data.frame())
    }
    df %>%
      arrange(pvalue) %>%
      slice_head(n = top) %>%
      mutate(
        Group = group_label,
        value = if (group_label == "Up") -log10(pvalue) else -(-log10(pvalue))
      )
  }

  up <- prepare_df(up_df, "Up")
  down <- prepare_df(down_df, "Down")

  df <- bind_rows(up, down)
  if (nrow(df) == 0) {
    message("No enrichment results to plot for: ", title)
    return(invisible(NULL))
  }

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

  save_ggsave_formats(p, out_base, formats, width = 8, height = 6)
}


run_pipeline <- function(
  infile,
  outdir,
  species = "mouse",
  gene_col = "gene_name",
  value_col = "fold_enrchment",
  p_col = "lg.p",
  up_cutoff = 1,
  p_cutoff = 0.05,
  top = 10,
  formats = "png"
) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  # fail fast on bad formats before any enrichment work starts
  validate_formats(formats)

  message("Defining Up / Down genes ...")
  genes <- define_up_down_genes(
    infile,
    gene_col,
    value_col,
    p_col,
    up_cutoff,
    p_cutoff
  )
  print(genes)
  writeLines(as.character(genes$up), file.path(outdir, "up_genes.txt"))
  writeLines(as.character(genes$down), file.path(outdir, "down_genes.txt"))

  for (type in c("go", "kegg")) {
    if (type == "go") {
      next
    }
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
      out_base = file.path(outdir, paste0(type, "_back_to_back")),
      formats = formats
    )
  }
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  res <- list(
    infile = NULL,
    outdir = NULL,
    species = "mouse",
    gene_col = "gene_name",
    value_col = "fold_enrchment",
    p_col = "lg.p",
    up_cutoff = 1,
    p_cutoff = 0.05,
    top = 10,
    formats = NULL
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "-i" || key == "--infile") {
      i <- i + 1
      res$infile <- args[[i]]
    } else if (key == "-o" || key == "--outdir") {
      i <- i + 1
      res$outdir <- args[[i]]
    } else if (key == "--species") {
      i <- i + 1
      res$species <- args[[i]]
    } else if (key == "--gene_col") {
      i <- i + 1
      res$gene_col <- args[[i]]
    } else if (key == "--value_col") {
      i <- i + 1
      res$value_col <- args[[i]]
    } else if (key == "--p_col") {
      i <- i + 1
      res$p_col <- args[[i]]
    } else if (key == "--up_cutoff") {
      i <- i + 1
      res$up_cutoff <- as.numeric(args[[i]])
    } else if (key == "--p_cutoff") {
      i <- i + 1
      res$p_cutoff <- as.numeric(args[[i]])
    } else if (key == "--top") {
      i <- i + 1
      res$top <- as.integer(args[[i]])
    } else if (key == "--format") {
      # repeatable: one output image per --format flag
      i <- i + 1
      res$formats <- c(res$formats, args[[i]])
    }
    i <- i + 1
  }
  res
}

main <- function() {
  args <- parse_args()
  if (is.null(args$infile) || is.null(args$outdir)) {
    stop("usage: go-kegg_back.r -i <infile> -o <outdir> [--species <human|mouse>] [--gene_col <col>] [--value_col <col>] [--p_col <col>] [--up_cutoff <num>] [--p_cutoff <num>] [--top <num>] [--format <png|pdf|svg|ps|eps|tif|tiff|jpg|jpeg>] (repeatable, default: png)")
  }
  # default keeps the historical png output; fail fast on unsupported formats
  formats <- if (is.null(args$formats)) "png" else args$formats
  validate_formats(formats)
  run_pipeline(
    infile = args$infile,
    outdir = args$outdir,
    species = args$species,
    gene_col = args$gene_col,
    value_col = args$value_col,
    p_col = args$p_col,
    up_cutoff = args$up_cutoff,
    p_cutoff = args$p_cutoff,
    top = args$top,
    formats = formats
  )
}

if (sys.nframe() == 0) {
  main()
}
