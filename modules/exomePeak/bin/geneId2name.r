#!/usr/bin/env Rscript
options(width = 60)

get_file_signature <- function(file_path) {
  info <- file.info(file_path)
  list(
    path = normalizePath(file_path, winslash = "/", mustWork = TRUE),
    size = info$size,
    mtime = as.numeric(info$mtime)
  )
}

save_cache <- function(cache_path, gene_map, signature) {
  saveRDS(list(signature = signature, gene_map = gene_map), cache_path)
}

load_cache_if_valid <- function(cache_path, gtf_signature) {
  if (!file.exists(cache_path)) {
    return(NULL)
  }
  cache <- tryCatch(readRDS(cache_path), error = function(e) NULL)
  if (is.null(cache)) {
    return(NULL)
  }
  cached_sig <- cache$signature
  if (is.null(cached_sig)) {
    return(NULL)
  }
  if (identical(cached_sig, gtf_signature)) {
    message("cache matched, loading directly: ", cache_path)
    return(cache$gene_map)
  }
  message("cache mismatch detected -> ignoring cache and rebuilding")
  NULL
}

parse_gtf_gene_map <- function(gtf_path) {
  gene_env <- new.env(hash = TRUE, parent = emptyenv())

  extract_attr <- function(info, key) {
    pattern <- paste0(key, " \\\"([^\\\"]+)\\\"")
    match <- regexpr(pattern, info, perl = TRUE)
    if (match[1] == -1) {
      return(NULL)
    }
    hit <- regmatches(info, match)
    sub(pattern, "\\1", hit, perl = TRUE)
  }

  con <- file(gtf_path, open = "r")
  on.exit(close(con), add = TRUE)

  repeat {
    lines <- readLines(con, n = 20000, warn = FALSE)
    if (length(lines) == 0) {
      break
    }
    for (line in lines) {
      if (startsWith(line, "#")) {
        next
      }
      if (!grepl("\tgene\t", line, fixed = TRUE)) {
        next
      }
      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(fields) < 9) {
        next
      }
      info <- fields[9]
      gid <- extract_attr(info, "gene_id")
      gname <- extract_attr(info, "gene_name")
      if (!is.null(gid) && !is.null(gname)) {
        gene_env[[gid]] <- gname
      }
    }
  }

  gene_map <- as.list(gene_env)
  message("GTF parse done, genes: ", length(gene_map))
  gene_map
}

load_gtf_gene_map <- function(gtf_path, cache_path = "gene_map.rds") {
  gtf_signature <- get_file_signature(gtf_path)
  gene_map <- load_cache_if_valid(cache_path, gtf_signature)
  if (!is.null(gene_map)) {
    return(gene_map)
  }
  message("begin parsing GTF to build gene map...")
  gene_map <- parse_gtf_gene_map(gtf_path)
  save_cache(cache_path, gene_map, gtf_signature)
  message("cache written: ", cache_path)
  gene_map
}

translate_gene_ids <- function(value, gene_map) {
  if (is.na(value) || is.null(value)) {
    return(value)
  }
  value <- trimws(as.character(value))
  if (value == "") {
    return(value)
  }
  ids <- strsplit(value, ",", fixed = TRUE)[[1]]
  names <- vapply(ids, function(gid) {
    if (!is.null(gene_map[[gid]])) gene_map[[gid]] else gid
  }, character(1))
  paste(names, collapse = ",")
}

convert_gene_id_to_name <- function(infile, gtf_path, gene_id_col = "name") {
  if (!file.exists(infile)) {
    stop("input file not found: ", infile)
  }
  if (!file.exists(gtf_path)) {
    stop("GTF not found: ", gtf_path)
  }
  cache_path <- file.path(dirname(infile), "gene_map.rds")
  gene_map <- load_gtf_gene_map(gtf_path, cache_path)
  df <- read.table(infile,
                   sep = "\t",
                   header = TRUE,
                   quote = "",
                   comment.char = "",
                   stringsAsFactors = FALSE,
                   check.names = FALSE)
  if (!(gene_id_col %in% colnames(df))) {
    stop("input missing column: ", gene_id_col)
  }
  df$gene_name <- vapply(df[[gene_id_col]], translate_gene_ids, gene_map = gene_map, character(1))
  df
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  res <- list(infile = NULL, gtf = NULL, outfile = NULL, gene_id_col = "name")
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "-i" || key == "--infile") {
      i <- i + 1
      res$infile <- args[[i]]
    } else if (key == "--gtf") {
      i <- i + 1
      res$gtf <- args[[i]]
    } else if (key == "-o" || key == "--outfile") {
      i <- i + 1
      res$outfile <- args[[i]]
    } else if (key == "--gene_id_col") {
      i <- i + 1
      res$gene_id_col <- args[[i]]
    }
    i <- i + 1
  }
  res
}

main <- function() {
  args <- parse_args()
  if (is.null(args$infile) || is.null(args$gtf) || is.null(args$outfile)) {
    stop("usage: geneId2name.R -i <infile> --gtf <gtf> -o <outfile> [--gene_id_col <col>]")
  }
  df <- convert_gene_id_to_name(args$infile, args$gtf, args$gene_id_col)
  write.table(df, file = args$outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  message("done: ", args$outfile)
}

if (sys.nframe() == 0) {
  main()
}
