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

# ---------------------------
# Utility helpers
# ---------------------------
ensure_dir <- function(dir){
  if(!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

msg <- function(...) message("[DESeq2Script] ", paste(..., collapse = " "))

# ---------------------------
# DESeq2Analysis: core function
# ---------------------------
DESeq2Analysis <- function(counts_df, colData, resfile = NULL,
                           normMethods = c("cpm","DESeq2"),
                           fitType = "mean", minReplicatesForReplace = 7,
                           parallel = FALSE){
  #' Perform DESeq2 differential expression and return results + dds
  #' counts_df: raw counts (rows = features, cols = samples)
  #' colData: data.frame with 'condition' factor matching columns order
  #' resfile: optional path to write combined results
  #' normMethods: "DESeq2" or "cpm"
  normMethods <- match.arg(normMethods)
  stopifnot(is.data.frame(counts_df) || is.matrix(counts_df))
  stopifnot(is.data.frame(colData))
  if(ncol(counts_df) != nrow(colData)){
    stop("Number of samples in count matrix (", ncol(counts_df),
         ") does not match number of rows in colData (", nrow(colData), ").")
  }

  msg("Building DESeqDataSet...")
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts_df),
                                colData = colData,
                                design = ~ condition)

  msg("Running DESeq() ...")
  dds <- tryCatch({
    DESeq(dds, fitType = fitType, minReplicatesForReplace = minReplicatesForReplace,
          parallel = parallel)
  }, error = function(e){
    stop("DESeq() failed: ", e$message)
  })

  # normalized counts
  if(normMethods == "DESeq2"){
    normalized_counts <- counts(dds, normalized = TRUE)
  } else if(normMethods == "cpm"){
    raw_counts <- counts(dds, normalized = FALSE)
    y <- edgeR::DGEList(counts = raw_counts)
    y <- edgeR::calcNormFactors(y, method = "TMM")
    normalized_counts <- edgeR::cpm(y, normalized.lib.sizes = TRUE, prior.count = 2)
  }

  normalized_counts <- as.data.frame(normalized_counts, check.names = FALSE)
  # DESeq2 results: experiment vs control (assumes factor levels: control, experiment)
  res <- results(dds, contrast = c("condition", "experiment", "control"))
  res_df <- as.data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

  # combine results and normalized counts; preserve rownames as feature id
  combined <- cbind(res_df, normalized_counts)
  rownames(combined) <- rownames(res_df)

  if(!is.null(resfile) && nzchar(resfile)){
    msg("Writing results to: ", resfile)
    ensure_dir(dirname(resfile))
    write.table(combined, file = resfile, sep = "\t", quote = FALSE, col.names = NA)
  }

  return(list(results = combined, dds = dds, normalized = normalized_counts))
}

# ---------------------------
# TEFilter: select TE-style rows by pattern
# ---------------------------
TEFilter <- function(df){
  #' Keep features whose rownames match "str1:str2:str3" (three colon-separated fields)
  stopifnot(!is.null(rownames(df)))
  df_TE <- df[grepl("^[^:]+:[^:]+:[^:]+$", rownames(df)), , drop = FALSE]
  return(df_TE)
}

# ---------------------------
# ScreenFeature: label up/down/none and export
# ---------------------------
ScreenFeature <- function(res_df, upfile = NULL, downfile = NULL, updownfile = NULL,
                          lfc_cut = 0.58, padj_cut = 0.05){
  #' Add a 'sig' column (up/down/none) based on thresholds and optionally write files.
  stopifnot(is.data.frame(res_df))
  # ensure padj and log2FoldChange exist
  if(!("padj" %in% colnames(res_df)) || !("log2FoldChange" %in% colnames(res_df))){
    stop("res_df must contain 'padj' and 'log2FoldChange' columns.")
  }

  # sort by padj asc then log2FC desc
  res_df <- res_df[order(res_df$padj, -res_df$log2FoldChange, na.last = TRUE), , drop = FALSE]

  # handle NA padj by assigning 'none'
  res_df$sig <- "none"
  res_df$sig[ which(!is.na(res_df$padj) & res_df$log2FoldChange >= lfc_cut & res_df$padj < padj_cut) ] <- "up"
  res_df$sig[ which(!is.na(res_df$padj) & res_df$log2FoldChange <= -lfc_cut & res_df$padj < padj_cut) ] <- "down"

  # move 'sig' after padj (if exists)
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
      msg("Wrote file:", path)
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

# ---------------------------
# Heatmap for up/down lists
# ---------------------------
updown_heatmap <- function(res_df, outfile, coldata = NULL, brewer_palette = "RdBu"){
  #' res_df: combined res + normalized counts; assumes counts start after columns from DESeq2 (i.e. >7)
  stopifnot(is.data.frame(res_df))
  # find numeric columns (counts) by detecting columns with numeric type
  heatmat <- as.matrix(res_df[, -c(1:7), drop = FALSE])
  # scale rows for heatmap
  if (nrow(heatmat) < 1) stop("Empty matrix for heatmap.")
  ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer_palette))(256)
  ramp <- rev(ramp) # original used reversed palette

  ensure_dir(dirname(outfile))
  # compute sensible row/column font sizes
  n_genes <- nrow(heatmat)
  row_font <- ifelse(n_genes <= 50, 8, max(3, 200 / n_genes))
  show_rows <- n_genes <=50
  # use pheatmap; write to file via png device to avoid device issues
  png(filename = outfile, width = 2000, height = 1500, res = 200)
  pheatmap::pheatmap(heatmat, color = ramp, scale = "row",
                     annotation_col = if(!is.null(coldata)) coldata else NA,
                     fontsize_col = 10,
                     cluster_cols = FALSE,
                     cluster_rows = FALSE,
                     show_rownames = show_rows,
                     fontsize_row = row_font,
                     border_color = NA)
  dev.off()
  msg("Heatmap saved to:", outfile)
}

# ---------------------------
# Volcano and MA plotting (shared helper to avoid duplication)
# ---------------------------
.add_gene_name <- function(res_df, mode = c("gene","TE"), geneAnnotation = NULL){
  mode <- match.arg(mode)
  res <- res_df %>% rownames_to_column(var = "index")
  if(mode == "TE"){
    res <- res %>% mutate(gene_name = stringr::word(index, 1, sep = ":"))
  } else {
    if(is.null(geneAnnotation) || !nzchar(geneAnnotation)){
      stop("Mode 'gene' requires geneAnnotation file.")
    }
    ann <- read.csv(geneAnnotation, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
    res <- res %>% mutate(
      gene_name = ifelse(str_count(index, ":") >= 2,
                         str_extract(index, "^[^:]+"),
                         ann$gene_name[match(index, rownames(ann))])
    )
  }
  return(res)
}

plot_volcano <- function(res_df, outfile, mode = c("gene","TE"), geneAnnotation = NULL,
                         lfc_cut = 0.58, padj_cut = 0.05, nlabel = 10, label.by = c("padj","log2FoldChange")){
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

  # prepare labels
  if(label.by == "padj"){
    up_genes <- res %>% filter(gene_status == "Upregulated genes") %>% arrange(padj) %>% head(nlabel)
    down_genes <- res %>% filter(gene_status == "Downregulated genes") %>% arrange(padj) %>% head(nlabel)
  } else {
    up_genes <- res %>% filter(gene_status == "Upregulated genes") %>% arrange(desc(log2FoldChange)) %>% head(nlabel)
    down_genes <- res %>% filter(gene_status == "Downregulated genes") %>% arrange(log2FoldChange) %>% head(nlabel)
  }

  # counts for legend
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

  ensure_dir(dirname(outfile))
  ggsave(outfile, plot = p, width = 8, height = 6, dpi = 300)
  msg("Volcano plot saved to:", outfile)
}

plot_MA <- function(res_df, outfile, mode = c("gene","TE"), geneAnnotation = NULL,
                    lfc_cut = 0.58, padj_cut = 0.05, nlabel = 10, label.by = c("padj","log2FoldChange")){
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

  ensure_dir(dirname(outfile))
  ggsave(outfile, plot = p, width = 8, height = 6, dpi = 300)
  msg("MA plot saved to:", outfile)
}

# ---------------------------
# Input readers for TEcount / TElocal
# ---------------------------
TEcount_read <- function(infile, control, experiment, min_mean = 5){
  #' Read TEcount (gene + TE) matrix, subset columns and split TE/Gene by rowname pattern
  stopifnot(file.exists(infile))
  TEcount <- read.csv(infile, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  # Subset columns by sample lists (preserve original column order of control then experiment)
  samples_sel <- c(control, experiment)
  if(!all(samples_sel %in% colnames(TEcount))){
    stop("Some samples in control/experiment not present in matrix columns.")
  }
  TEcount <- TEcount[, samples_sel, drop = FALSE]

  df_TE <- TEcount[grepl("^[^:]+:[^:]+:[^:]+$", rownames(TEcount)), , drop = FALSE]
  df_TE <- df_TE[rowMeans(df_TE, na.rm = TRUE) > min_mean, , drop = FALSE]

  df_Gene <- TEcount[!grepl("^[^:]+:[^:]+:[^:]+$", rownames(TEcount)), , drop = FALSE]
  df_Gene <- df_Gene[rowMeans(df_Gene, na.rm = TRUE) > min_mean, , drop = FALSE]

  return(list(Gene_TE = TEcount, TE = df_TE, Gene = df_Gene))
}

TElocal_read <- function(infile, control, experiment, min_mean = 5){
  stopifnot(file.exists(infile))
  TElocal <- read.csv(infile, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  samples_sel <- c(control, experiment)
  if(!all(samples_sel %in% colnames(TElocal))){
    stop("Some samples in control/experiment not present in matrix columns.")
  }
  TElocal <- TElocal[, samples_sel, drop = FALSE]
  df_TE <- TElocal[grepl("^[^:]+:[^:]+:[^:]+:[^:]+$", rownames(TElocal)), , drop = FALSE]
  df_TE <- df_TE[rowMeans(df_TE, na.rm = TRUE) > min_mean, , drop = FALSE]
  return(df_TE)
}

# ---------------------------
# Argument parsing and main workflow (keeps your original CLI but more robust)
# ---------------------------
parser <- ArgumentParser(description = 'DESeq2 analysis for TEcount and TElocal')
parser$add_argument('-m', '--mode', type = 'character', required = TRUE, choices = c("TEcount", "TElocal"))
parser$add_argument('-i', '--matrix', type = 'character', required = TRUE)
parser$add_argument('-g', '--group', type = 'character', required = TRUE,
                    help = 'path to group file, tab-separated, two columns: sample,group')
parser$add_argument('-p', '--pattern', type = 'character', required = TRUE, nargs = 2,
                    help = '--pattern control experiment')
parser$add_argument('-f', '--figure', type = 'character', required = TRUE, nargs = "+",
                    choices = c("heatmap","volcano","pca"))
parser$add_argument('-o', '--outdir', type = 'character', required = TRUE)
parser$add_argument('-a', '--annotation', type = 'character', required = TRUE,
                    help = 'annotation file: gene_id, gene_type, gene_name (tab)')
parser$add_argument('-Tcm', '--TEcountMode', type = 'character', nargs = "+",
                    choices = c("all","Gene_TE","TE","Gene"), default = "all")
args <- parser$parse_args()

# read group and extract sample lists
df_group <- read.csv(args$group, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
if(!all(c("sample","group") %in% colnames(df_group))){
  stop("Group file must contain columns 'sample' and 'group'.")
}
controlStr <- args$pattern[1]
experimentStr <- args$pattern[2]
control <- as.vector(df_group$sample[df_group$group == controlStr])
experiment <- as.vector(df_group$sample[df_group$group == experimentStr])
if(length(control) == 0 || length(experiment) == 0){
  stop("No samples found for control or experiment. Check group file and --pattern.")
}
# prepare annotation for heatmap align withe read function
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
ensure_dir(outdir)
msg("Mode:", args$mode, "Matrix:", args$matrix, "Outdir:", outdir)

# main branch
if(args$mode == "TEcount"){
  dfList <- TEcount_read(args$matrix, control, experiment)

  if("pca" %in% args$figure){
    df_for_pca <- dfList[["Gene"]]
    res_for_pca <- DESeq2Analysis(df_for_pca, colData, normMethods = "cpm")
    # res_for_pca$results contains combined table; pass normalized to PCA plotting
    pca_out <- file.path(outdir, "DESeq2", "plot", "cpmPCA.png")
    ensure_dir(dirname(pca_out))
    # create a simple PCA plot using normalized counts
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
    ggsave(pca_out, plot = p, width = 8, height = 6, dpi = 300)
    msg("PCA saved to:", pca_out)
  }

  # which TEcount parts to process
  modes_to_run <- if("all" %in% args$TEcountMode) names(dfList) else args$TEcountMode

  for(dataType in modes_to_run){
    if(!(dataType %in% names(dfList))){
      msg("Skipping unknown dataType:", dataType); next
    }
    df <- dfList[[dataType]]
    outfile <- file.path(outdir, "DESeq2", paste0("TEcount_", dataType, ".tsv"))
    ensure_dir(dirname(outfile))
    resDds <- DESeq2Analysis(df, colData, resfile = outfile, normMethods = "DESeq2")
    res_comb <- resDds$results
    # TEFilter only for Gene_TE type (to extract TE-like rows)
    if(dataType == "Gene_TE"){
      res_comb_TE <- TEFilter(res_comb)
    } else {
      res_comb_TE <- res_comb
    }

    # plotting
    if("heatmap" %in% args$figure){
      heat_dir <- file.path(outdir, "DESeq2", "heatmap")
      ensure_dir(heat_dir)
      # create up/down sets and heatmaps
      sf <- ScreenFeature(res_comb_TE,
                          upfile = file.path(outdir,"DESeq2","upDown",paste0("TEcount_",dataType,"_up.tsv")),
                          downfile = file.path(outdir,"DESeq2","upDown",paste0("TEcount_",dataType,"_down.tsv")),
                          updownfile = file.path(outdir,"DESeq2","upDown",paste0("TEcount_",dataType,"_updown.tsv")))
      if(nrow(sf$up) >= 2){
        up_file <- file.path(heat_dir, paste0("TEcount_", dataType, "_up.png"))
        updown_heatmap(sf$up, up_file, coldata = sample_anno)
      } else msg("Less than 2 up genes; skip heatmap for up.")
      if(nrow(sf$down) >= 2){
        down_file <- file.path(heat_dir, paste0("TEcount_", dataType, "_down.png"))
        updown_heatmap(sf$down, down_file, coldata = sample_anno)
      } else msg("Less than 2 down genes; skip heatmap for down.")
      if(nrow(rbind(sf$up, sf$down)) >= 2){
        all_file <- file.path(heat_dir, paste0("TEcount_", dataType, "_updown.png"))
        updown_heatmap(rbind(sf$up, sf$down), all_file, coldata = sample_anno)
      } else msg("Less than 2 dysregulated genes; skip combined heatmap.")
    }

    if("volcano" %in% args$figure){
      vol_dir <- file.path(outdir, "DESeq2", "volcano")
      ensure_dir(vol_dir)
      # decide gene vs TE based on name pattern or dataType
      is_gene_mode <- grepl("Gene", dataType)
      if(is_gene_mode){
        vol_out <- file.path(vol_dir, paste0("TEcount_", dataType, "_volcano.png"))
        ma_out <- file.path(vol_dir, paste0("TEcount_", dataType, "_MA.png"))
        plot_volcano(res_comb_TE, vol_out, mode = "gene", geneAnnotation = args$annotation)
        plot_MA(res_comb_TE, ma_out, mode = "gene", geneAnnotation = args$annotation)
      } else {
        vol_out <- file.path(vol_dir, paste0("TEcount_", dataType, "_volcano.png"))
        ma_out <- file.path(vol_dir, paste0("TEcount_", dataType, "_MA.png"))
        plot_volcano(res_comb_TE, vol_out, mode = "TE")
        plot_MA(res_comb_TE, ma_out, mode = "TE")
      }
    }
  } # end loop dataType

} else if(args$mode == "TElocal"){
  df <- TElocal_read(args$matrix, control, experiment)
  outfile <- file.path(outdir, "DESeq2", "TElocal_TE.tsv")
  resDds <- DESeq2Analysis(df, colData, resfile = outfile, normMethods = "DESeq2")
  res_comb <- resDds$results
  # plotting
  if("heatmap" %in% args$figure){
    plot_dir <- file.path(outdir, "DESeq2", "heatmap")
    ensure_dir(plot_dir)
    sf <- ScreenFeature(res_comb,
                        upfile = file.path(outdir,"DESeq2","upDown","TElocal_TE_up.tsv"),
                        downfile = file.path(outdir,"DESeq2","upDown","TElocal_TE_down.tsv"),
                        updownfile = file.path(outdir,"DESeq2","upDown","TElocal_TE_updown.tsv"))
    if(nrow(sf$up) >= 2) updown_heatmap(sf$up, file.path(plot_dir, "TElocal_up.png"), coldata = sample_anno)
    if(nrow(sf$down) >= 2) updown_heatmap(sf$down, file.path(plot_dir, "TElocal_down.png"), coldata = sample_anno)
  }
  if("volcano" %in% args$figure){
    vol_dir <- file.path(outdir, "DESeq2", "volcano")
    ensure_dir(vol_dir)
    plot_volcano(res_comb, file.path(vol_dir, "TElocal_volcano.png"), mode = "TE")
    plot_MA(res_comb, file.path(vol_dir, "TElocal_MA.png"), mode = "TE")
  }
} else {
  stop("Unsupported mode: ", args$mode)
}

msg("DESeq2 workflow finished.")
