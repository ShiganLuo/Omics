#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  pkgs <- c("DESeq2","edgeR","pheatmap","ggplot2","ggrepel","RColorBrewer",
            "argparse","dplyr","stringr","tibble")
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing_pkgs)){
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Install them before running the script.")
  }
  library(DESeq2)
  library(edgeR)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(argparse)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ensure_dir <- function(dir){
  if(!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

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

VALID_FORMATS <- c("png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg", "jpeg")

#' Validate image output formats, fail fast on anything unsupported
#'
#' @param formats Character vector of requested formats.
#' @return The validated formats (invisibly usable as a vector).
check_formats <- function(formats) {
  formats <- as.character(formats)
  bad <- setdiff(formats, VALID_FORMATS)
  if (length(bad)) {
    stop(sprintf("Unsupported image format(s): %s; valid: %s",
                 paste(bad, collapse = ", "), paste(VALID_FORMATS, collapse = " ")))
  }
  formats
}

#' Open a graphics device for one output format (base graphics / pheatmap)
#'
#' Device-per-extension helper; sizing is in inches for every format so the
#' geometry is identical across devices.
#'
#' @param path Output file path (with extension).
#' @param fmt One of VALID_FORMATS.
#' @param width_in Device width in inches.
#' @param height_in Device height in inches.
#' @param res Resolution for raster devices. Default 300.
open_image_device <- function(path, fmt, width_in, height_in, res = 300) {
  switch(fmt,
    png = grDevices::png(path, width = width_in, height = height_in,
                         units = "in", res = res),
    pdf = grDevices::pdf(path, width = width_in, height = height_in),
    svg = grDevices::svg(path, width = width_in, height = height_in),
    ps  = grDevices::postscript(path, width = width_in, height = height_in),
    eps = grDevices::postscript(path, width = width_in, height = height_in,
                                onefile = FALSE, horizontal = FALSE),
    tif =, tiff = grDevices::tiff(path, width = width_in, height = height_in,
                                  units = "in", res = res),
    jpg =, jpeg = grDevices::jpeg(path, width = width_in, height = height_in,
                                  units = "in", res = res),
    stop(sprintf("Unsupported image format(s): %s; valid: %s", fmt,
                 paste(VALID_FORMATS, collapse = " "))))
}

#' Perform DESeq2 differential expression and return results + dds
#'
#' @param counts_df Data frame or matrix. Raw counts where rows are features (genes) 
#'   and columns are samples.
#' @param colData Data frame. Metadata with a 'condition' factor column. 
#'   Rows must match the order of columns in counts_df.
#' @param resfile Character. Optional. File path to save the combined results 
#'   (normalized counts + statistics).
#' @param normMethods Character. Normalization method to return: "DESeq2" or "cpm". 
#'   Default is "cpm".
#' @param fitType Character. The type of fitting used for dispersion. Default is "mean".
#' @param minReplicatesForReplace Numeric. Threshold for Outlier replacement in DESeq2. 
#'   Default is 7.
#' @param parallel Logical. Whether to use parallel processing for DESeq2. 
#'   Default is FALSE.
#'
#' @return A list containing three elements:
#'   - `results`: Combined data frame of DESeq2 statistics and normalized counts.
#'   - `dds`: The processed DESeqDataSet object.
#'   - `normalized`: The data frame of normalized counts.
#' @export
DESeq2Analysis <- function(counts_df, colData, resfile = NULL,
                           normMethods = c("cpm","DESeq2"),
                           fitType = "mean", minReplicatesForReplace = 7,
                           parallel = FALSE){
  normMethods <- match.arg(normMethods)
  if (!(is.data.frame(counts_df) || is.matrix(counts_df))) {
    log_msg("ERROR","'counts_df' must be a data.frame or a matrix.", quit = TRUE)
  }
  if (!is.data.frame(colData)) {
    log_msg("ERROR","'colData' must be a data.frame.", quit = TRUE)
  }

  if(ncol(counts_df) != nrow(colData)){
    log_msg("ERROR","Number of samples in count matrix (", ncol(counts_df),
         ") does not match number of rows in colData (", nrow(colData), ").", quit = TRUE)
  }

  log_msg("INFO","Building DESeqDataSet...")
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts_df),
                                colData = colData,
                                design = ~ condition)

  log_msg("INFO","Running DESeq() ...")
  dds <- tryCatch({
    DESeq(dds, fitType = fitType, minReplicatesForReplace = minReplicatesForReplace,
          parallel = parallel)
  }, error = function(e){
    log_msg("ERROR","DESeq() failed: ", e$message, quit = TRUE)
  })

  if(normMethods == "DESeq2"){
    normalized_counts <- counts(dds, normalized = TRUE)
  } else if(normMethods == "cpm"){
    raw_counts <- counts(dds, normalized = FALSE)
    y <- edgeR::DGEList(counts = raw_counts)
    y <- edgeR::calcNormFactors(y, method = "TMM")
    normalized_counts <- edgeR::cpm(y, normalized.lib.sizes = TRUE, prior.count = 2)
  }

  normalized_counts <- as.data.frame(normalized_counts, check.names = FALSE)
  res <- results(dds, contrast = c("condition", "experiment", "control"))
  res_df <- as.data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

  combined <- cbind(res_df, normalized_counts)
  rownames(combined) <- rownames(res_df)

  if(!is.null(resfile) && nzchar(resfile)){
    log_msg("INFO","Writing results to: ", resfile)
    ensure_dir(dirname(resfile))
    write.table(combined, file = resfile, sep = "\t", quote = FALSE, col.names = NA)
  }

  return(list(results = combined, dds = dds, normalized = normalized_counts))
}

#' Filter Transposable Element (TE) entries from a data frame
#'
#' This function filters rows based on a specific naming convention in the rownames.
#' It expects rownames to follow a triple-colon format (e.g., "string1:string2:string3"),
#' which is a common convention for representing TE genomic coordinates or classifications.
#'
#' @param df A data frame or matrix. The rownames must be strings in the format "str1:str2:str3".
#'
#' @return A data frame containing only the rows that match the "str1:str2:str3" pattern.
#' @export
TEFilter <- function(df, pattern = "^[^:]+:[^:]+:[^:]+$"){
  if(is.null(rownames(df))) {
    log_msg("ERROR","the df rownames can't be null, must be colon-separated fields")
  }
  df_TE <- df[grepl(pattern, rownames(df)), , drop = FALSE]
  
  df_TE
}

#' Screen and Categorize Differential Features
#'
#' This function filters a differential expression result data frame based on 
#' log2FoldChange and adjusted p-value (padj) thresholds. It labels features 
#' as "up", "down", or "none" and optionally exports these lists to files.
#'
#' @param res_df Data frame. The output from a differential expression analysis 
#'   (e.g., DESeq2). Must contain 'log2FoldChange' and 'padj' columns.
#' @param upfile Character. Optional path to save the up-regulated features (TSV).
#' @param downfile Character. Optional path to save the down-regulated features (TSV).
#' @param updownfile Character. Optional path to save both up and down-regulated features (TSV).
#' @param lfc_cut Numeric. The absolute log2FoldChange threshold. Default is 0.58 
#'   (approx. 1.5-fold change).
#' @param padj_cut Numeric. The adjusted p-value (FDR) threshold. Default is 0.05.
#'
#' @return A list containing three data frames:
#'   \itemize{
#'     \item \code{up}: Features significantly up-regulated.
#'     \item \code{down}: Features significantly down-regulated.
#'     \item \code{all}: The input data frame with an added 'sig' column, sorted by padj.
#'   }
#' @export
ScreenFeature <- function(res_df, upfile = NULL, downfile = NULL, updownfile = NULL,
                          lfc_cut = 0.58, padj_cut = 0.05){
  if (!is.data.frame(res_df)) {
    log_msg("ERROR","The res_df must be dataframe",quit = TRUE)
  }
  if(!("padj" %in% colnames(res_df)) || !("log2FoldChange" %in% colnames(res_df))){
    stop("res_df must contain 'padj' and 'log2FoldChange' columns.")
  }

  res_df <- res_df[order(res_df$padj, -res_df$log2FoldChange, na.last = TRUE), , drop = FALSE]

  res_df$sig <- "none"
  res_df$sig[ which(!is.na(res_df$padj) & res_df$log2FoldChange >= lfc_cut & res_df$padj < padj_cut) ] <- "up"
  res_df$sig[ which(!is.na(res_df$padj) & res_df$log2FoldChange <= -lfc_cut & res_df$padj < padj_cut) ] <- "down"

  if("padj" %in% colnames(res_df)){
    cols <- colnames(res_df)
    cols <- append(setdiff(cols, "sig"), "sig", after = match("padj", cols))
    res_df <- res_df[, cols, drop = FALSE]
  }

  # write files if requested
  ensure_out <- function(path, df){
    if(!is.null(path) && nzchar(path)){
      ensure_dir(dirname(path))
      write.table(df, file = path, sep = "\t", quote = FALSE, col.names = NA)
      log_msg("INFO","Wrote file:", path)
    }
  }

  res_up <- subset(res_df, sig == "up")
  res_down <- subset(res_df, sig == "down")
  res_updown <- subset(res_df, sig %in% c("up","down"))

  ensure_out(upfile, res_up)
  ensure_out(downfile, res_down)
  ensure_out(updownfile, res_updown)

  return(list(up = res_up, down = res_down, all = res_df))
}

#' Generate Heatmap for Up/Down Regulated Features
#'
#' This function creates a row-scaled (Z-score) heatmap from differential expression 
#' results. It automatically handles color palettes, dynamic font sizing based on 
#' the number of genes, and optional sample annotations.
#'
#' @param res_df Data frame. The input data where the first 7 columns are typically 
#'   statistics (log2FC, padj, etc.) and subsequent columns are normalized counts.
#' @param out_base Character. Output path without extension; one file is
#'   written per requested format as \code{<out_base>.<fmt>} (e.g.
#'   \code{"heatmap"} -> \code{"heatmap.png"}).
#' @param coldata Data frame. Optional metadata for columns (samples) to display 
#'   annotation bars at the top of the heatmap. Default is NULL.
#' @param brewer_palette Character. A color palette name from RColorBrewer (e.g., "RdBu", "PiYG"). 
#'   Default is "RdBu".
#' @param formats Character vector of image formats (png/pdf/svg/ps/eps/tif/
#'   tiff/jpg/jpeg). Default \code{"png"}.
#'
#' @details 
#' The function assumes that normalized counts start from the 8th column. 
#' It applies row-scaling (\code{scale = "row"}) to highlight relative expression 
#' changes across samples. Row names (gene IDs) are automatically hidden if the 
#' number of genes exceeds 50 to maintain legibility.
#'
#' @return None. Saves one image per requested format to \code{<out_base>.<fmt>}.
#' @export
updown_heatmap <- function(res_df, out_base, coldata = NULL, brewer_palette = "RdBu", formats = "png"){
  
  if (!is.data.frame(res_df)) {
    log_msg("ERROR","The res_df must be dataframe", quit = TRUE)
  }
  heatmat <- as.matrix(res_df[, -c(1:7), drop = FALSE])
  if (nrow(heatmat) < 1) stop("Empty matrix for heatmap.")
  ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer_palette))(256)
  ramp <- rev(ramp)

  formats <- check_formats(formats)
  ensure_dir(dirname(out_base))
  n_genes <- nrow(heatmat)
  row_font <- ifelse(n_genes <= 50, 8, max(3, 200 / n_genes))
  show_rows <- n_genes <=50
  for (fmt in formats) {
    outfile <- sprintf("%s.%s", out_base, fmt)
    open_image_device(outfile, fmt, width_in = 10, height_in = 7.5, res = 200)
    pheatmap::pheatmap(heatmat, color = ramp, scale = "row",
                       annotation_col = if(!is.null(coldata)) coldata else NA,
                       fontsize_col = 10,
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_rownames = show_rows,
                       fontsize_row = row_font,
                       border_color = NA)
    dev.off()
    log_msg("INFO", "Heatmap saved to:", outfile)
  }
}

#' Add Gene Names with Fallback Logic
#'
#' This function appends a `gene_name` column to the results. It extracts names 
#' from TE-style identifiers directly or looks them up in an annotation file 
#' for standard genes. If a lookup fails, it preserves the original ID.
#'
#' @param res_df Data frame. The results table with feature IDs as rownames.
#' @param mode Character. Either "gene" or "TE".
#' @param geneAnnotation Character. Path to a TSV annotation file.
#'
#' @return A data frame with an additional `gene_name` column.
#' @export
.add_gene_name <- function(res_df, mode = c("gene", "TE"), geneAnnotation = NULL) {
  mode <- match.arg(mode)

  res <- res_df %>% tibble::rownames_to_column(var = "index")

  if (mode == "TE") {
    res <- res %>% 
      dplyr::mutate(gene_name = stringr::word(index, 1, sep = ":"))

  } else {
    if (!is.null(geneAnnotation)) {
      ann <- read.csv(geneAnnotation, sep = "\t", header = TRUE, 
                      row.names = 1, stringsAsFactors = FALSE)

      res <- res %>% dplyr::mutate(
        gene_name = ifelse(
          stringr::str_count(index, ":") >= 2,
          stringr::str_extract(index, "^[^:]+"),
          ann$gene_name[match(index, rownames(ann))]
        ),
        gene_name = dplyr::coalesce(gene_name, index)
      )
    } else {
      log_msg("WARN", "Mode 'gene' requires geneAnnotation. Using original IDs as names.")
      res <- res %>% dplyr::mutate(gene_name = index)
    }
  }

  res
}

#' Generate a Volcano Plot for Differential Expression Analysis
#'
#' This function creates a volcano plot visualizing the relationship between 
#' fold change and statistical significance. It labels top differentially expressed 
#' features and colors points based on user-defined thresholds.
#'
#' @param res_df Data frame. The results from differential expression analysis 
#'   (containing at least `log2FoldChange` and `padj` columns).
#' @param out_base Character. Output path without extension; one file is
#'   written per requested format as \code{<out_base>.<fmt>} (e.g.
#'   \code{"volcano"} -> \code{"volcano.png"}).
#' @param mode Character. Either "gene" or "TE", used for gene name mapping via `.add_gene_name`.
#' @param geneAnnotation Character. Path to the annotation file required if `mode` is "gene".
#' @param lfc_cut Numeric. Absolute log2 Fold Change threshold for significance. Default is 0.58.
#' @param padj_cut Numeric. Adjusted p-value threshold for significance. Default is 0.05.
#' @param nlabel Integer. Number of top up-regulated and down-regulated features to label. Default is 10.
#' @param label.by Character. Criterion to select features for labeling: "padj" (most significant) 
#'   或 "log2FoldChange" (highest magnitude of change).
#'
#' @details 
#' The plot uses a $-\log_{10}$ transformation for the y-axis (adjusted p-value). 
#' Points are categorized into "Upregulated", "Downregulated", or "Other" based on 
#' the `lfc_cut` and `padj_cut`. The legend automatically includes the count (n) 
#' for each category. Text labeling is handled by `ggrepel` to avoid overlapping.
#'
#' @param formats Character vector of image formats (png/pdf/svg/ps/eps/tif/
#'   tiff/jpg/jpeg). Default \code{"png"}.
#'
#' @return None. Saves one high-resolution (300 DPI) image per requested
#'   format to \code{<out_base>.<fmt>}.
#' @export
plot_volcano <- function(res_df, out_base, mode = c("gene","TE"), geneAnnotation = NULL,
                         lfc_cut = 0.58, padj_cut = 0.05, nlabel = 10, label.by = c("padj","log2FoldChange"),
                         formats = "png"){
  mode <- match.arg(mode)
  label.by <- match.arg(label.by)
  res <- .add_gene_name(res_df, mode = mode, geneAnnotation = geneAnnotation)
  res <- res %>% mutate(padj = as.numeric(padj), log2FoldChange = as.numeric(log2FoldChange))
  res <- res %>% filter(!is.na(padj))
  res <- res %>% mutate(
    gene_status = case_when(
      log2FoldChange > lfc_cut & -log10(padj) > -log10(padj_cut) ~ "Upregulated genes",
      log2FoldChange < -lfc_cut & -log10(padj) > -log10(padj_cut) ~ "Downregulated genes",
      TRUE ~ "Other genes"
    )
  )

  if(label.by == "padj"){
    up_genes <- res %>% filter(gene_status == "Upregulated genes") %>% arrange(padj) %>% head(nlabel)
    down_genes <- res %>% filter(gene_status == "Downregulated genes") %>% arrange(padj) %>% head(nlabel)
  } else {
    up_genes <- res %>% filter(gene_status == "Upregulated genes") %>% arrange(desc(log2FoldChange)) %>% head(nlabel)
    down_genes <- res %>% filter(gene_status == "Downregulated genes") %>% arrange(log2FoldChange) %>% head(nlabel)
  }

  status_counts <- res %>% group_by(gene_status) %>% summarise(n = n())
  new_labels <- setNames(paste0(status_counts$gene_status, " (n = ", status_counts$n, ")"),
                         status_counts$gene_status)

  p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = gene_status), alpha = 0.7) +
    scale_color_manual(values = c("Upregulated genes" = "red", "Downregulated genes" = "blue", "Other genes" = "grey"),
                       labels = new_labels) +
    ggrepel::geom_text_repel(data = up_genes, aes(label = gene_name), color = "#F59494", size = 3) +
    ggrepel::geom_text_repel(data = down_genes, aes(label = gene_name), color = "#93ACF6", size = 3) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dotted") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dotted") +
    labs(x = expression(log[2]("FoldChange")), y = expression(-log[10]("adjusted p-value"))) +
    theme_classic() +
    theme(legend.title = element_blank())

  formats <- check_formats(formats)
  ensure_dir(dirname(out_base))
  for (fmt in formats) {
    outfile <- sprintf("%s.%s", out_base, fmt)
    ggsave(outfile, plot = p, width = 8, height = 6, dpi = 300)
    log_msg("INFO","Volcano plot saved to:", outfile)
  }
}

#' Generate an MA Plot for Differential Expression Analysis
#'
#' This function creates an MA plot (Minus-Average plot) visualizing the relationship between 
#' the average expression (baseMean) and the log2 Fold Change. It helps in identifying 
#' expression-dependent biases and highlighting significant differentially expressed features.
#'
#' @param res_df Data frame. The differential expression results (e.g., from DESeq2). 
#'   Must contain `baseMean`, `log2FoldChange`, and `padj` columns.
#' @param out_base Character. Output path without extension; one file is
#'   written per requested format as \code{<out_base>.<fmt>} (e.g.
#'   \code{"MA_plot"} -> \code{"MA_plot.png"}).
#' @param mode Character. Either "gene" or "TE", used for identifier mapping via `.add_gene_name`.
#' @param geneAnnotation Character. Optional path to the gene annotation file for "gene" mode.
#' @param lfc_cut Numeric. The log2 Fold Change threshold for significance. Default is 0.58.
#' @param padj_cut Numeric. The adjusted p-value threshold for significance. Default is 0.05.
#' @param nlabel Integer. Number of top features to label in both up and down directions. Default is 10.
#' @param label.by Character. Criterion to prioritize features for labeling: "padj" (significance) 
#'   or "log2FoldChange" (effect size).
#'
#' @details 
#' The x-axis represents \eqn{\log_2(baseMean + 1)}, providing a measure of the average abundance 
#' of a feature across all samples. The y-axis represents the \eqn{\log_2(Fold Change)}. 
#' Significant features are colored (Red for Up, Blue for Down) based on both the 
#' `lfc_cut` and `padj_cut` thresholds.
#'
#' @param formats Character vector of image formats (png/pdf/svg/ps/eps/tif/
#'   tiff/jpg/jpeg). Default \code{"png"}.
#'
#' @return None. Saves one high-resolution image per requested format to
#'   \code{<out_base>.<fmt>}.
#' @export
plot_MA <- function(res_df, out_base, mode = c("gene","TE"), geneAnnotation = NULL,
                    lfc_cut = 0.58, padj_cut = 0.05, nlabel = 10, label.by = c("padj","log2FoldChange"),
                    formats = "png"){
  mode <- match.arg(mode)
  label.by <- match.arg(label.by)
  res <- .add_gene_name(res_df, mode = mode, geneAnnotation = geneAnnotation)
  res <- res %>% mutate(padj = as.numeric(padj), log2FoldChange = as.numeric(log2FoldChange), baseMean = as.numeric(baseMean))
  res <- res %>% filter(!is.na(padj))
  res <- res %>% mutate(
    gene_status = case_when(
      log2FoldChange > lfc_cut & -log10(padj) > -log10(padj_cut) ~ "Upregulated genes",
      log2FoldChange < -lfc_cut & -log10(padj) > -log10(padj_cut) ~ "Downregulated genes",
      TRUE ~ "Other genes"
    )
  )

  if(label.by == "padj"){
    up_genes <- res %>% filter(gene_status == "Upregulated genes") %>% arrange(padj) %>% head(nlabel)
    down_genes <- res %>% filter(gene_status == "Downregulated genes") %>% arrange(padj) %>% head(nlabel)
  } else {
    up_genes <- res %>% filter(gene_status == "Upregulated genes") %>% arrange(desc(log2FoldChange)) %>% head(nlabel)
    down_genes <- res %>% filter(gene_status == "Downregulated genes") %>% arrange(log2FoldChange) %>% head(nlabel)
  }

  status_counts <- res %>% group_by(gene_status) %>% summarise(n = n())
  new_labels <- setNames(paste0(status_counts$gene_status, " (n = ", status_counts$n, ")"),
                         status_counts$gene_status)

  p <- ggplot(res, aes(x = log2(baseMean + 1), y = log2FoldChange)) +
    geom_point(aes(color = gene_status), alpha = 0.7) +
    scale_color_manual(values = c("Upregulated genes" = "red", "Downregulated genes" = "blue", "Other genes" = "grey"),
                       labels = new_labels) +
    ggrepel::geom_text_repel(data = up_genes, aes(label = gene_name), color = "#F59494", size = 3) +
    ggrepel::geom_text_repel(data = down_genes, aes(label = gene_name), color = "#93ACF6", size = 3) +
    geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dotted") +
    labs(x = "log2(Expression in baseMean)", y = "log2(FoldChange)") +
    theme_classic() +
    theme(legend.title = element_blank())

  formats <- check_formats(formats)
  ensure_dir(dirname(out_base))
  for (fmt in formats) {
    outfile <- sprintf("%s.%s", out_base, fmt)
    ggsave(outfile, plot = p, width = 8, height = 6, dpi = 300)
    log_msg("INFO","MA plot saved to:", outfile)
  }
}

#' Read and Classify TEcount Expression Matrix
#'
#' This function reads a raw count matrix from a file, subsets the samples based on 
#' experimental design, separates Transposable Elements (TEs) from standard Genes 
#' using regular expressions, and filters out low-expression features.
#'
#' @param infile Character. Path to the input count matrix file (TSV format).
#' @param control Character vector. Names of control samples (must match column names in the matrix).
#' @param experiment Character vector. Names of experiment samples (must match column names in the matrix).
#' @param min_mean Numeric. The threshold for row mean counts. Rows with an average 
#'   count less than or equal to this value will be filtered out. Default is 5.
#'
#' @details 
#' The function performs the following steps:
#' \enumerate{
#'   \item \strong{Sample Selection}: Reorganizes matrix columns to follow the order of 
#'         \code{control} followed by \code{experiment}.
#'   \item \strong{Feature Classification}: Rows matching \code{te_pattern} are TE
#'         entries; all other rows are gene entries.
#'   \item \strong{Low-count Filtering}: Calculates the mean for each row and removes
#'         noise data falling below the \code{min_mean} threshold.
#' }
#'
#' @param infile Character. Path to the input count matrix file (TSV format).
#' @param control Character vector. Names of control samples (must match column names in the matrix).
#' @param experiment Character vector. Names of experiment samples (must match column names in the matrix).
#' @param min_mean Numeric. The threshold for row mean counts. Rows with an average
#'   count less than or equal to this value will be filtered out. Default is 5.
#' @param te_pattern Character. Regex identifying TE rows by colon-separated fields.
#'   TEcount subfamily rows have 3 fields ("str1:str2:str3"); TElocal locus rows
#'   have 4 fields. Default matches TEcount.
#'
#' @return A list containing three elements:
#'   \itemize{
#'     \item \code{Gene_TE}: The raw subsetted matrix containing all selected samples
#'           (before low-expression filtering).
#'     \item \code{TE}: A data frame containing only the filtered TE entries.
#'     \item \code{Gene}: A data frame containing only the filtered Gene entries.
#'   }
#' @export
TEcount_read <- function(infile, control, experiment, min_mean = 5,
                         te_pattern = "^[^:]+:[^:]+:[^:]+$"){
  if (!file.exists(infile)) {
    log_msg("ERROR",paste0(infile," don't exist"), quit = TRUE)
  }
  TEcount <- read.csv(infile, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  samples_sel <- c(control, experiment)
  if(!all(samples_sel %in% colnames(TEcount))){
    log_msg("ERROR", "Some samples in control/experiment not present in matrix columns.", quit = TRUE)
  }
  TEcount <- TEcount[, samples_sel, drop = FALSE]

  df_TE <- TEcount[grepl(te_pattern, rownames(TEcount)), , drop = FALSE]
  df_TE <- df_TE[rowMeans(df_TE, na.rm = TRUE) > min_mean, , drop = FALSE]

  df_Gene <- TEcount[!grepl(te_pattern, rownames(TEcount)), , drop = FALSE]
  df_Gene <- df_Gene[rowMeans(df_Gene, na.rm = TRUE) > min_mean, , drop = FALSE]

  list(Gene_TE = TEcount, TE = df_TE, Gene = df_Gene)
}

#' Read and Filter TElocal Locus-Specific Expression Matrix
#'
#' This function reads a TElocal output matrix, selects specific samples for
#' downstream analysis, and splits rows into locus-specific Transposable Element
#' (TE) entries (quadruple-colon naming convention) and gene entries, filtering
#' both subsets for low expression. Behaves like \code{\link{TEcount_read}} but
#' defaults to the TElocal 4-field row naming convention.
#'
#' @param infile Character. Path to the input TElocal matrix file (TSV format).
#' @param control Character vector. Names of control samples to include.
#' @param experiment Character vector. Names of experiment samples to include.
#' @param min_mean Numeric. The threshold for row mean counts. Features with an
#'   average count across selected samples less than or equal to this value
#'   will be removed. Default is 5.
#' @param te_pattern Character. Regex identifying TE loci rows. Default matches
#'   the TElocal naming convention with four segments
#'   (e.g. \code{TE_Name:Family:Class:Locus}).
#'
#' @return A list containing three elements:
#'   \itemize{
#'     \item \code{Gene_TE}: The raw subsetted matrix containing all selected samples
#'           (before low-expression filtering).
#'     \item \code{TE}: A data frame containing only the filtered TE loci entries.
#'     \item \code{Gene}: A data frame containing only the filtered Gene entries.
#'   }
#' @export
TElocal_read <- function(infile, control, experiment, min_mean = 5,
                         te_pattern = "^[^:]+:[^:]+:[^:]+:[^:]+$"){
  TEcount_read(infile = infile, control = control, experiment = experiment,
               min_mean = min_mean, te_pattern = te_pattern)
}

#' Read and Filter a Gene Expression Count Matrix
#'
#' This function imports a count matrix from a tab-separated file, subsets the 
#' columns based on specified control and experiment samples, and filters out 
#' low-expression features based on a mean count threshold. If non-integer 
#' (decimal) values are detected, they are automatically rounded to the nearest 
#' integer to ensure compatibility with DESeq2.
#'
#' @param infile Character. The file path to the input count matrix (TSV format).
#' @param control Character vector. Names of the control samples to be extracted.
#' @param experiment Character vector. Names of the experiment samples to be extracted.
#' @param min_mean Numeric. The threshold for the average expression across selected samples. 
#'   Rows with a mean value less than or equal to this threshold will be removed. Default is 5.
#'
#' @details 
#' The function ensures that all requested samples exist in the matrix headers 
#' before processing. It reorders the columns to match the order of \code{control} 
#' and \code{experiment}. Automatic rounding is performed if any cell in the 
#' subsetted matrix is not an integer.
#'
#' @return A data frame containing the subsetted, filtered, and rounded expression counts.
#' @export
count_read <- function(infile, control, experiment, min_mean = 5) {
  if (!file.exists(infile)) {
    log_msg("ERROR", paste0(infile, " does not exist"), quit = TRUE)
  }
  
  count <- read.csv(infile, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  
  samples_sel <- c(control, experiment)
  if (!all(samples_sel %in% colnames(count))) {
    log_msg("ERROR", "Some samples in control/experiment not present in matrix columns.", quit = TRUE)
  }
  count <- count[, samples_sel, drop = FALSE]

  if (!all(count == floor(count), na.rm = TRUE)) {
    log_msg("INFO", "Non-integer values detected in count matrix. Rounding to nearest integer for DESeq2 compatibility.")
    count <- round(count)
  }

  count <- count[rowMeans(count, na.rm = TRUE) > min_mean, , drop = FALSE]
  
  print(head(count))
  return(count)
}

#' DESeq2 Analysis CLI Parser
#'
#' Defines the command-line interface for performing differential expression
#' analysis on TEcount or TElocal output matrices.
#'
#' @section Arguments:
#' \itemize{
#'   \item \bold{-m, --mode}: \emph{Character}. The input data type. 
#'     \code{TEcount} for gene/subfamily level and \code{TElocal} for locus-specific level or \code{Count} for gene.
#'   \item \bold{-i, --matrix}: \emph{Character}. Path to the raw count matrix file (TSV).
#'   \item \bold{-g, --group}: \emph{Character}. Path to the sample metadata file. 
#'     Must be a tab-separated file with columns: \code{sample} and \code{group}.
#'   \item \bold{-p, --pattern}: \emph{Character (Two values)}. Defines the comparison logic. 
#'     Example: \code{--pattern control experiment}.
#'   \item \bold{-f, --figure}: \emph{Character (Multiple allowed)}. Types of plots to generate. 
#'     Options: \code{heatmap}, \code{volcano}, \code{pca}.
#'   \item \bold{-o, --outdir}: \emph{Character}. Directory path where results and plots will be saved.
#'   \item \bold{-a, --annotation}: \emph{Character}. Path to the gene annotation file 
#'     containing \code{gene_id}, \code{gene_type}, and \code{gene_name}.
#'   \item \bold{-Tcm, --TEcountMode}: \emph{Character (Multiple allowed)}. Only applicable 
#'     for \code{TEcount} mode. Specifies which data subsets to analyze. 
#'     Default is \code{all}.
#'   \item \bold{--format}: \emph{Character (Multiple allowed)}. Output image
#'     format(s), repeatable (\code{--format png --format pdf}); no short flag
#'     because \code{-f} is taken by \code{--figure}. Valid: png pdf svg ps eps
#'     tif tiff jpg jpeg (pgf/raw/rgba are matplotlib-only and fail fast).
#'     Default is \code{png}.
#' }
parser <- ArgumentParser(description = 'DESeq2 analysis for TEcount and TElocal')
parser$add_argument('-m', '--mode', type = 'character', required = TRUE, choices = c("TEcount", "TElocal", "Count"))
parser$add_argument('-i', '--matrix', type = 'character', required = TRUE)
parser$add_argument('-g', '--group', type = 'character', required = TRUE,
                    help = 'path to group file, tab-separated, two columns: sample,group')
parser$add_argument('-p', '--pattern', type = 'character', required = TRUE, nargs = 2,
                    help = '--pattern control experiment')
parser$add_argument('-f', '--figure', type = 'character', required = TRUE, nargs = "+",
                    choices = c("heatmap","volcano","pca"))
parser$add_argument('-o', '--outdir', type = 'character', required = TRUE)
parser$add_argument('-px', '--prefix', type = 'character', required = FALSE, default = "", help="file prefix for output files")
parser$add_argument('-a', '--annotation', type = 'character', required = FALSE, default = NULL,
                    help = 'annotation file: gene_id, gene_type, gene_name (tab)')
parser$add_argument('-Tcm', '--TEcountMode', type = 'character', nargs = "+",
                    choices = c("all","Gene_TE","TE","Gene"), default = "all")
parser$add_argument('--format', type = 'character', action = 'append',
                    help = 'output image format, repeatable (default: png); valid: png pdf svg ps eps tif tiff jpg jpeg')
args <- parser$parse_args()
fmts <- check_formats(if (is.null(args$format)) "png" else as.character(unlist(args$format)))

df_group <- read.csv(args$group, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
if(!all(c("sample","group") %in% colnames(df_group))){
  log_msg("ERROR","Group file must contain columns 'sample' and 'group'.",quit = TRUE)
}
controlStr <- args$pattern[1]
experimentStr <- args$pattern[2]
control <- as.vector(df_group$sample[df_group$group == controlStr])
experiment <- as.vector(df_group$sample[df_group$group == experimentStr])
if(length(control) == 0 || length(experiment) == 0){
  log_msg("ERROR","No samples found for control or experiment. Check group file and --pattern.", quit = TRUE)
}
sample_anno <- subset(df_group, group %in% c(controlStr, experimentStr))
sample_anno <- sample_anno[order(factor(sample_anno$group, levels = c(controlStr, experimentStr))), ]
rownames(sample_anno) <- sample_anno$sample
sample_anno$sample <- NULL
sample_anno$group <- factor(sample_anno$group)

colData <- data.frame(condition = factor(c(rep('control', length(control)),
                                          rep('experiment', length(experiment))),
                                        levels = c('control','experiment')))

colDataPca <- factor(c(rep(controlStr, length(control)), rep(experimentStr, length(experiment))))

outdir <- args$outdir
prefix <- args$prefix
ensure_dir(outdir)
log_msg("INFO","Mode:", args$mode, "Matrix:", args$matrix, "Outdir:", outdir, "Prefix:", prefix, "Figures:", paste(args$figure, collapse = ","))

if(args$mode %in% c("TEcount", "TElocal")){
  te_pattern <- if(args$mode == "TElocal") "^[^:]+:[^:]+:[^:]+:[^:]+$" else "^[^:]+:[^:]+:[^:]+$"
  mode_tag <- args$mode
  if(args$mode == "TElocal"){
    dfList <- TElocal_read(args$matrix, control, experiment, te_pattern = te_pattern)
  } else {
    dfList <- TEcount_read(args$matrix, control, experiment, te_pattern = te_pattern)
  }

  if("pca" %in% args$figure){
    df_for_pca <- dfList[["Gene"]]
    res_for_pca <- DESeq2Analysis(df_for_pca, colData, normMethods = "cpm")
    pca_out <- file.path(outdir, "PCA", paste0(prefix, ".cpmPCA"))
    ensure_dir(dirname(pca_out))
    normalized <- res_for_pca$normalized
    pca_mat <- prcomp(t(normalized))
    pca_df <- data.frame(PC1 = pca_mat$x[,1], PC2 = pca_mat$x[,2],
                         Sample = colnames(normalized),
                         Group = colDataPca)
    p <- ggplot(pca_df, aes(PC1, PC2, color = Group, label = Sample)) +
      geom_point(size = 3) + geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
      labs(x = paste0("PC1: ", round(summary(pca_mat)$importance[2,1]*100,2), "%"),
           y = paste0("PC2: ", round(summary(pca_mat)$importance[2,2]*100,2), "%"),
           title = "PCA Plot of CPM Values") +
      theme_minimal()
    for (fmt in fmts) {
      out_file <- sprintf("%s.%s", pca_out, fmt)
      ggsave(out_file, plot = p, width = 8, height = 6, dpi = 300)
      log_msg("INFO","PCA saved to:", out_file)
    }
  }

  modes_to_run <- if("all" %in% args$TEcountMode) names(dfList) else args$TEcountMode

  for(dataType in modes_to_run){
    if(!(dataType %in% names(dfList))){
      log_msg("INFO","Skipping unknown dataType:", dataType); next
    }
    df <- dfList[[dataType]]
    outfile <- file.path(outdir, paste0(prefix, ".", mode_tag, "_", dataType, ".tsv"))
    ensure_dir(dirname(outfile))
    resDds <- DESeq2Analysis(df, colData, resfile = outfile, normMethods = "DESeq2")
    res_comb <- resDds$results
    if(dataType == "Gene_TE"){
      res_comb_TE <- TEFilter(res_comb, pattern = te_pattern)
    } else {
      res_comb_TE <- res_comb
    }

    if("heatmap" %in% args$figure){
      heat_dir <- file.path(outdir, "heatmap")
      ensure_dir(heat_dir)
      sf <- ScreenFeature(res_comb_TE,
                          upfile = file.path(outdir,"upDown",paste0(prefix,".",mode_tag,"_",dataType,"_up.tsv")),
                          downfile = file.path(outdir,"upDown",paste0(prefix,".",mode_tag,"_",dataType,"_down.tsv")),
                          updownfile = file.path(outdir,"upDown",paste0(prefix,".",mode_tag,"_",dataType,"_updown.tsv")))
      if(nrow(sf$up) >= 2){
        up_file <- file.path(heat_dir, paste0(prefix,".", mode_tag, "_", dataType, "_up"))
        updown_heatmap(sf$up, up_file, coldata = sample_anno, formats = fmts)
      } else log_msg("INFO","Less than 2 up genes; skip heatmap for up.")
      if(nrow(sf$down) >= 2){
        down_file <- file.path(heat_dir, paste0(prefix,".", mode_tag, "_", dataType, "_down"))
        updown_heatmap(sf$down, down_file, coldata = sample_anno, formats = fmts)
      } else log_msg("INFO","Less than 2 down genes; skip heatmap for down.")
      if(nrow(rbind(sf$up, sf$down)) >= 2){
        all_file <- file.path(heat_dir, paste0(prefix,".", mode_tag, "_", dataType, "_updown"))
        updown_heatmap(rbind(sf$up, sf$down), all_file, coldata = sample_anno, formats = fmts)
      } else log_msg("INFO","Less than 2 dysregulated genes; skip combined heatmap.")
    }

    if("volcano" %in% args$figure){
      vol_dir <- file.path(outdir, "volcano")
      ensure_dir(vol_dir)
      is_gene_mode <- grepl("Gene", dataType)
      if(is_gene_mode){
        vol_out <- file.path(vol_dir, paste0(prefix,".",mode_tag,"_", dataType, "_volcano"))
        ma_out <- file.path(vol_dir, paste0(prefix,".",mode_tag,"_", dataType, "_MA"))
        plot_volcano(res_comb_TE, vol_out, mode = "gene", geneAnnotation = args$annotation, formats = fmts)
        plot_MA(res_comb_TE, ma_out, mode = "gene", geneAnnotation = args$annotation, formats = fmts)
      } else {
        vol_out <- file.path(vol_dir, paste0(prefix,".",mode_tag,"_", dataType, "_volcano"))
        ma_out <- file.path(vol_dir, paste0(prefix,".",mode_tag,"_", dataType, "_MA"))
        plot_volcano(res_comb_TE, vol_out, mode = "TE", formats = fmts)
        plot_MA(res_comb_TE, ma_out, mode = "TE", formats = fmts)
      }
    }
  }

} else if(args$mode == "Count"){
  df <- count_read(args$matrix, control, experiment)
  outfile <- file.path(outdir, paste0(prefix, ".Count.tsv"))
  resDds <- DESeq2Analysis(df, colData, resfile = outfile, normMethods = "DESeq2")
  res_comb <- resDds$results
  if("heatmap" %in% args$figure){
    plot_dir <- file.path(outdir, "heatmap")
    ensure_dir(plot_dir)
    sf <- ScreenFeature(res_comb,
                        upfile = file.path(outdir,"upDown",paste0(prefix,".Count_up.tsv")),
                        downfile = file.path(outdir,"upDown",paste0(prefix,".Count_down.tsv")),
                        updownfile = file.path(outdir,"upDown",paste0(prefix,".Count_updown.tsv")))
    if(nrow(sf$up) >= 2) updown_heatmap(sf$up, file.path(plot_dir, paste0(prefix, ".Count_up")), coldata = sample_anno, formats = fmts)
    if(nrow(sf$down) >= 2) updown_heatmap(sf$down, file.path(plot_dir, paste0(prefix, ".Count_down")), coldata = sample_anno, formats = fmts)
  }
  if("volcano" %in% args$figure){
    vol_dir <- file.path(outdir, "volcano")
    ensure_dir(vol_dir)
    plot_volcano(res_comb, file.path(vol_dir, paste0(prefix, ".Count_volcano")), mode = "gene", formats = fmts)
    plot_MA(res_comb, file.path(vol_dir, paste0(prefix, ".Count_MA")), mode = "gene", formats = fmts)
  }
} else {
  log_msg("ERROR","Unsupported mode:", args$mode, quit = TRUE)
}
outfile_done = file.path(outdir, "DESeq2.done")
writeLines("DESeq2 workflow completed successfully.", con = outfile_done)
log_msg("INFO","DESeq2 workflow finished.")
