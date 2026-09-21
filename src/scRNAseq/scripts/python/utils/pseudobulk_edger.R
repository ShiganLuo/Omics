#!/usr/bin/env Rscript
## ============================================================
## Pseudo-bulk DEG via edgeR (adapted from DEG.py)
## ============================================================
## Reads pseudo-bulk counts CSVs (rows=samples, cols=genes) produced by
## pseudobulk_prepare.py and runs edgeR QL-FTest pairwise comparisons.
##
## Usage:
##   Rscript pseudobulk_edger.R \
##     --counts-dir /path/to/pseudobulk_counts \
##     --out-dir /path/to/pseudobulk_DEG \
##     --samples-ovaries luanchao-21310-10XSC3,luanchao-11238-10XSC3 \
##     --comparisons-uterus spec.yaml
## ============================================================

suppressPackageStartupMessages({
  library(edgeR)
})

run_edger_pairwise <- function(
  counts_file,
  output_prefix,
  group_a,
  group_b
) {
  if (!file.exists(counts_file)) {
    cat("[skip] missing:", counts_file, "\n")
    return(NULL)
  }

  df <- read.csv(counts_file, row.names = 1, check.names = FALSE)
  # sample column was eaten as index; recover from rownames
  df$sample <- rownames(df)
  meta_cols <- c("sample", "cell_type", "n_cells", "group")
  gene_cols <- setdiff(colnames(df), meta_cols)

  if (length(gene_cols) < 100) {
    cat("[skip] too few genes:", counts_file, "\n")
    return(NULL)
  }

  # counts: rows=samples, cols=genes; edgeR wants genes x samples -> transpose
  counts <- t(as.matrix(df[, gene_cols, drop = FALSE]))
  storage.mode(counts) <- "numeric"

  keep_samples <- rownames(df)[df$sample %in% c(group_a, group_b)]
  if (length(keep_samples) < 2) {
    cat("[skip] too few samples for", group_b, "vs", group_a,
        "in", counts_file, "\n")
    return(NULL)
  }
  counts <- counts[, keep_samples, drop = FALSE]
  group <- factor(df[keep_samples, "sample"], levels = c(group_a, group_b))

  y <- DGEList(counts = counts, group = group)
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = FALSE]
  if (nrow(y) < 10) {
    cat("[skip] too few expressed genes:", counts_file, "\n")
    return(NULL)
  }

  y <- calcNormFactors(y)

  # With only 1 sample per group (no biological replicates), use
  # exactTest with a fixed common dispersion.
  # prior.count=2 shrinks log2FC for low-count genes (avoids ±13 artifacts).
  # dispersion=0.5 is realistic for pseudobulk from single-cell data.
  et <- exactTest(y, pair = c(group_a, group_b),
                  dispersion = 0.5, prior.count = 2)
  tt <- topTags(et, n = Inf, sort.by = "none")$table

  tt$gene <- rownames(tt)
  tt$cell_type <- df$cell_type[1]
  tt$comparison <- paste0(group_b, "_vs_", group_a)
  tt$log2FC <- tt$logFC
  # exactTest uses PValue/FDR; glmQLTest uses PValue/qval. Handle both.
  tt$pval <- tt$PValue
  if ("qval" %in% colnames(tt)) {
    tt$pval_adj <- tt$qval
  } else {
    tt$pval_adj <- tt$FDR
  }
  tt$significant <- (tt$pval_adj < 0.05) & (abs(tt$log2FC) > 1)

  write.csv(tt, paste0(output_prefix, "_DEG.csv"))
  cat(sprintf("  %-50s %4d genes, %4d significant\n",
              paste0(basename(counts_file), ": ",
                     group_b, " vs ", group_a),
              nrow(tt), sum(tt$significant)))
  invisible(tt)
}

# ------------------------------------------------------------
# CLI
# ------------------------------------------------------------
parser <- argparse::ArgumentParser()
parser$add_argument("--counts-dir", required = TRUE)
parser$add_argument("--out-dir", required = TRUE)
parser$add_argument("--label", required = TRUE,
                    help = "Dataset label, e.g. 'ovaries' or 'uterus'")
parser$add_argument("--comparisons", required = TRUE,
                    help = "Comma-separated group_a:group_b pairs, e.g. 'A:B,C:D'")
args <- parser$parse_args()

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

# Parse comparisons
comp_pairs <- list()
for (pair in strsplit(args$comparisons, ",")[[1]]) {
  ab <- strsplit(pair, ":")[[1]]
  if (length(ab) == 2) {
    comp_pairs[[length(comp_pairs) + 1]] <- list(a = ab[1], b = ab[2])
  }
}
if (length(comp_pairs) == 0) {
  stop("No valid comparisons parsed from: ", args$comparisons)
}

cat(sprintf("\n=== %s DEG ===\n", toupper(args$label)))
cat(sprintf("Counts dir: %s\n", args$counts_dir))
cat(sprintf("Out dir:    %s\n", args$out_dir))
cat(sprintf("Comparisons: %s\n\n", args$comparisons))

# Process all counts files
counts_files <- list.files(
  args$counts_dir,
  pattern = paste0("^", args$label, "_.*_counts\\.csv$"),
  full.names = TRUE
)
cat(sprintf("Found %d counts files\n", length(counts_files)))

for (f in counts_files) {
  ct_raw <- sub(paste0("^", args$label, "_"), "",
                 basename(f), perl = TRUE)
  ct_raw <- sub("_counts\\.csv$", "", ct_raw)
  ct_safe <- gsub("[ /]", "-", ct_raw)
  for (cp in comp_pairs) {
    out_prefix <- file.path(args$out_dir, paste0(args$label, "_",
                                                paste0(cp$a, "_vs_", cp$b),
                                                "_", ct_safe))
    run_edger_pairwise(f, out_prefix, cp$a, cp$b)
  }
}

# Summary
deg_files <- list.files(args$out_dir, pattern = "_DEG\\.csv$", full.names = TRUE)
cat(sprintf("\nTotal DEG outputs: %d\n", length(deg_files)))
sig_count <- 0
for (f in deg_files) {
  d <- read.csv(f)
  n <- sum(d$significant)
  if (n > 0) {
    cat(sprintf("  %s: %d significant\n", basename(f), n))
    sig_count <- sig_count + 1
  }
}
cat(sprintf("\nFiles with >=1 significant DEG: %d\n", sig_count))