#!/usr/bin/env Rscript

## =========================
## GO and KEGG enrichment for SV-affected gene lists
## Derived from modules/function/bin/go-kegg.r (same engine style), but
## consumes a flat gene list (SVs have no up/down direction) extracted from
## a VEP-annotated SV tab file's SYMBOL column.
## =========================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(org.Mmu.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(argparse)
})

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

## Read VEP tab: locate the #Uploaded_variation header, take unique SYMBOLs
read_vep_genes <- function(tab_file, impact_filter = c("HIGH","MODERATE","LOW")) {
  header_line <- NULL
  con <- file(tab_file, "r")
  on.exit(close(con))
  i <- 0L
  while (TRUE) {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    if (startsWith(line, "#Uploaded_variation")) { header_line <- i; break }
    i <- i + 1L
  }
  if (is.null(header_line))
    log_msg("ERROR", "No #Uploaded_variation header in ", tab_file, quit = TRUE)

  df <- read.delim(tab_file, sep = "\t", skip = header_line,
                   header = TRUE, quote = "", comment.char = "",
                   check.names = FALSE)
  if (!"IMPACT" %in% colnames(df) || !"SYMBOL" %in% colnames(df))
    log_msg("ERROR", "VEP tab lacks IMPACT/SYMBOL columns: ", tab_file, quit = TRUE)

  df_f <- df[df$IMPACT %in% impact_filter, ]
  genes <- unique(na.omit(df_f$SYMBOL))
  # VEP writes "-" for missing SYMBOL; drop it
  genes <- genes[genes != "-" & nzchar(trimws(genes))]
  list(genes = genes, total = nrow(df), filtered = nrow(df_f))
}

run_go_kegg <- function(
  genes,
  species = c("human", "mouse", "macaque"),
  type = c("go", "kegg")
) {
  species <- match.arg(species)
  type <- match.arg(type)
  empty <- data.frame(
    ID = character(), Description = character(), GeneRatio = character(),
    BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
    qvalue = numeric(), geneID = character(), Count = integer())

  if (length(genes) == 0) return(empty)

  if (species == "human") {
    OrgDb <- org.Hs.eg.db; kegg_org <- "hsa"
  } else if (species == "macaque") {
    OrgDb <- org.Mmu.eg.db; kegg_org <- "mcc"
  } else {
    OrgDb <- org.Mm.eg.db; kegg_org <- "mmu"
  }

  gene_df <- tryCatch(
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb, drop = TRUE),
    error = function(e) data.frame(SYMBOL = character(0), ENTREZID = character(0))
  )
  entrez <- unique(gene_df$ENTREZID)
  log_msg("INFO", length(genes), " symbols -> ", length(entrez),
          " ENTREZID via ", species, " OrgDb")
  if (length(entrez) == 0) return(empty)

  if (type == "go") {
    res <- enrichGO(
      gene = entrez, OrgDb = OrgDb, ont = "BP",
      pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE
    )
  } else {
    res <- enrichKEGG(
      gene = entrez, organism = kegg_org,
      pvalueCutoff = 0.05, qvalueCutoff = 0.05
    )
  }
  if (is.null(res)) return(empty)
  out <- as.data.frame(res@result)
  if (nrow(out) > 0 && "geneID" %in% colnames(out)) {
    # readable=FALSE for KEGG: map ENTREZ back to symbols for readability
    if (type == "kegg") {
      id2sym <- setNames(gene_df$SYMBOL, gene_df$ENTREZID)
      out$geneID <- vapply(strsplit(as.character(out$geneID), "/"),
                           function(ids) paste(unique(na.omit(id2sym[ids])), collapse = "/"),
                           character(1))
    }
  }
  out
}

clean_pathway_description <- function(description) {
  description <- as.character(description)
  description[is.na(description) | !nzchar(trimws(description))] <- "Unnamed pathway"
  sub(
    " - (Mus musculus \\\\(house mouse\\\\)|Homo sapiens \\\\(human\\\\)|Macaca mulatta \\\\(Rhesus macaque\\\\))$",
    "", trimws(description), perl = TRUE
  )
}

wrap_plot_labels <- function(labels, width = 42) {
  newline <- intToUtf8(10)
  vapply(
    labels,
    function(label) {
      label <- gsub("\\\\n", " ", label, fixed = FALSE)
      paste(strwrap(label, width = width), collapse = newline)
    },
    character(1)
  )
}

plot_bar <- function(res_df, top = 10, title = "", outfile) {
  if (is.null(res_df) || nrow(res_df) == 0) {
    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No enriched terms", size = 6) +
      labs(title = title, x = NULL, y = NULL) +
      theme_void() + theme(plot.title = element_text(hjust = 0.5))
    ggsave(outfile, p, width = 10, height = 4, dpi = 300, bg = "white")
    return(invisible(NULL))
  }
  top_df <- res_df %>% arrange(pvalue) %>% slice_head(n = top)
  top_df$plot_label <- wrap_plot_labels(clean_pathway_description(top_df$Description))
  p <- ggplot(top_df, aes(x = reorder(plot_label, -log10(pvalue)),
                          y = -log10(pvalue))) +
    geom_col(fill = "#4A72B7", width = 0.7) +
    coord_flip() +
    labs(title = title, x = NULL, y = expression(-log[10](pvalue))) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.text.y = element_text(size = 10, lineheight = 0.9),
      plot.title = element_text(hjust = 0.5)
    )
  n_labels <- length(unique(top_df$plot_label))
  ggsave(outfile, p, width = 12,
         height = max(6, min(18, 2.5 + 0.32 * n_labels)),
         dpi = 300, bg = "white")
}

## =========================
## CLI
## =========================
parser <- ArgumentParser(description = "GO/KEGG enrichment for SV-affected genes (VEP tab input)")

parser$add_argument("-i", "--input", required = TRUE, type = "character",
                    help = "VEP-annotated SV tab file (#Uploaded_variation header)")
parser$add_argument("-o", "--outdir", required = TRUE, type = "character",
                    help = "Output directory")
parser$add_argument("-s", "--species", default = "macaque", type = "character",
                    choices = c("human", "mouse", "macaque"),
                    help = "Species: human, mouse or macaque (default: macaque)")
parser$add_argument("--impact", type = "character", default = "HIGH,MODERATE,LOW",
                    help = "Comma-separated IMPACT levels to keep (default: HIGH,MODERATE,LOW)")
parser$add_argument("--top", type = "integer", default = 10,
                    help = "Top N pathways to plot (default: 10)")
parser$add_argument("--format", type = "character", default = "png",
                    choices = c("png", "pdf"), help = "Output image format (default: png)")
args <- parser$parse_args()

## =========================
## Main
## =========================
log_msg("INFO", "Input: ", args$input)
log_msg("INFO", "Outdir: ", args$outdir)
log_msg("INFO", "Species: ", args$species)
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

impact <- strsplit(args$impact, ",")[[1]]
info <- read_vep_genes(args$input, impact)
log_msg("INFO", "Total SVs: ", info$total,
        "  Filtered (IMPACT in ", paste(impact, collapse = ","), "): ",
        info$filtered, "  Unique genes: ", length(info$genes))

# Gene list is always written (empty list -> empty enrichment outputs, not a
# missing file: downstream rules declare these outputs)
writeLines(info$genes, file.path(args$outdir, "sv_genes.txt"))

for (type in c("go", "kegg")) {
  message("Running ", toupper(type), " enrichment ...")
  res <- run_go_kegg(info$genes, args$species, type)
  write.csv(res, file.path(args$outdir, paste0(type, "_enrichment.csv")),
            row.names = FALSE)
  plot_bar(res, top = args$top,
           title = paste0(toupper(type), " enrichment (SV-affected genes)"),
           outfile = file.path(args$outdir, paste0(type, "_enrichment.", args$format)))
}

log_msg("INFO", "SV GO/KEGG analysis completed.")
