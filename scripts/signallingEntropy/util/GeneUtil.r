box::use(
  data.table[...],
  bi = biomaRt,
  stats[setNames]
)
#' Convert gene symbols in an expression matrix to Entrez Gene IDs
#'
#' @param expr_mat Gene expression matrix/data.frame/data.table.
#'   Rows (or first column for data.table) are gene symbols,
#'   columns are samples.
#' @param species Species name: "human", "mouse", or "rat".
#' @param keep_multi Logical; keep one-to-many symbol → Entrez mappings.
#' @param verbose Logical; print mapping statistics.
#'
#' @return
#' Expression matrix (same type as input) with Entrez Gene IDs
#' replacing gene symbols.
#'
#' @export
symbol_expr_to_entrez <- function(
  expr_mat,
  species = c("human", "mouse", "rat"),
  keep_multi = TRUE,
  verbose = TRUE
) {
  # ---- load required packages ----
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package 'AnnotationDbi' is required but not installed.")
  }
  # load organism DBs
  species <- match.arg(species)
  db <- switch(
    species,
    human = { require(org.Hs.eg.db); org.Hs.eg.db },
    mouse = { require(org.Mm.eg.db); org.Mm.eg.db },
    rat   = { require(org.Rn.eg.db); org.Rn.eg.db }
  )

  # ---- detect input type ----
  is_dt <- is.data.table(expr_mat)

  # ---- extract gene symbols ----
  if (is_dt) {
    gene_col <- colnames(expr_mat)[1]
    symbols <- expr_mat[[gene_col]]
  } else {
    if (is.null(rownames(expr_mat))) {
      stop("expr_mat must have gene symbols as row names.")
    }
    symbols <- rownames(expr_mat)
  }

  # ---- symbol -> Entrez mapping ----
  map_df <- AnnotationDbi::select(
    db,
    keys = unique(symbols),
    keytype = "SYMBOL",
    columns = "ENTREZID"
  )
  colnames(map_df) <- c("symbol", "entrez_id")
  if (!keep_multi) {
    map_df <- map_df[!duplicated(map_df$symbol), ]
  }

  # ---- print mapping stats ----
  if (verbose) {
    n_input   <- length(unique(symbols))
    n_mapped  <- length(unique(map_df$symbol[!is.na(map_df$entrez_id)]))
    n_unmap   <- n_input - n_mapped
    n_multi   <- sum(table(map_df$symbol) > 1)

    message(
      sprintf(
        "[Entrez mapping] input=%d, mapped=%d, unmapped=%d, multi-mapped=%d",
        n_input, n_mapped, n_unmap, n_multi
      )
    )
  }

  # ---- drop unmapped genes ----
  map_df <- map_df[!is.na(map_df$entrez_id), ]

  # ---- expand expression matrix ----
  if (is_dt) {
    dt <- expr_mat[
      map_df$symbol,
      on = setNames(colnames(expr_mat)[1], "symbol"),
      allow.cartesian = TRUE
    ]
    data.table::setnames(dt, colnames(expr_mat)[1], "entrez_id")
    return(dt)
  } else {
    expr_out <- expr_mat[map_df$symbol, , drop = FALSE]
    rownames(expr_out) <- map_df$entrez_id
    return(expr_out)
  }
}


#' Build adjacency matrix from STRING interaction file
#'
#' Convert a STRING three-column interaction file
#' (protein1, protein2, combined_score) into a weighted
#' gene-level adjacency matrix.
#'
#' Optionally supports conversion from ENSEMBL Protein IDs
#' (e.g. 9606.ENSP00000000233) to Entrez Gene IDs for
#' multi-species analyses.
#'
#' @param file_path Path to STRING interaction file.
#'   Must contain columns: protein1, protein2, combined_score.
#'
#' @param convert_to_entrez Logical; whether to convert ENSEMBL
#'   Protein IDs to Entrez Gene IDs.
#'
#' @param species Species identifier, either NCBI taxonomy ID
#'   (e.g. 9606, 10090) or species name ("human", "mouse").
#'
#' @param score_scale Numeric vector of length 2 specifying
#'   the target range for edge weights (e.g. c(0, 1)).
#'
#' @return
#' A symmetric weighted adjacency matrix with gene IDs
#' as row and column names.
#'
#' @export
string_to_adj <- function(file_path,
                          convert_to_entrez = TRUE,
                          species = 9606,
                          score_scale = c(0, 1)) {
  # 1. 读 STRING 文件
  dt <- fread(file_path, header = TRUE)
  dt[, c("protein1", "protein2") :=
      lapply(.SD, sub, pattern = "^[0-9]+\\.", replacement = ""),
   .SDcols = c("protein1", "protein2")]

  if(!all(c("protein1","protein2","combined_score") %in% names(dt))) {
    log_msg("ERROR", "STRING file must have columns: protein1, protein2, combined_score",quit=TRUE)
  }

  # 2. 可选 ID 转换
  if(convert_to_entrez) {
    log_msg("INFO","Converting ENSEMBL Protein IDs to Entrez Gene IDs via biomaRt...")

    # 2a. 判断 species
    mart_dataset <- switch(as.character(species),
                           "9606" = "hsapiens_gene_ensembl",
                           "human" = "hsapiens_gene_ensembl",
                           "10090" = "mmusculus_gene_ensembl",
                           "mouse" = "mmusculus_gene_ensembl",
                           log_msg("ERROR","Unsupported species. Please provide Taxonomy ID or human/mouse",quit=TRUE))

    tryCatch({
      ensembl <- bi$useMart("ensembl", dataset = mart_dataset)
    }, error = function(e) {
      log_msg("ERROR","Failed to connect to Ensembl BioMart. Please check your internet connection.",quit=TRUE)
    })

    res <- bi$getBM(
      attributes = c("ensembl_peptide_id","ensembl_gene_id"),
      mart = ensembl,
    )
    # 2b. 批量获取映射
    mapping <- bi$getBM(
      attributes = c("ensembl_peptide_id","ensembl_gene_id"),
      filters = "ensembl_peptide_id",
      values = unique(c(dt$protein1, dt$protein2)),
      mart = ensembl
    )
    
    map_tbl <- setNames(mapping$ensembl_gene_id, mapping$ensembl_peptide_id)

    dt[, protein1 := map_tbl[protein1]]
    dt[, protein2 := map_tbl[protein2]]

    # 删除没有映射的行
    dt <- dt[!is.na(protein1) & !is.na(protein2)]
  }

  # 3. 统一 protein 列名
  proteins <- unique(c(dt$protein1, dt$protein2))

  # 4. 创建空矩阵
  adj.m <- matrix(0, nrow = length(proteins), ncol = length(proteins),
                  dimnames = list(proteins, proteins))

  # 5. 填充邻接矩阵
  dt[, score := combined_score]

  # 线性缩放 combined_score 到 score_scale
  min_s <- min(dt$score, na.rm = TRUE)
  max_s <- max(dt$score, na.rm = TRUE)
  target_min <- score_scale[1]
  target_max <- score_scale[2]
  dt[, score := (score - min_s) / (max_s - min_s) * (target_max - target_min) + target_min]

  # 填充矩阵（无向图）
  for(i in seq_len(nrow(dt))) {
    p1 <- dt$protein1[i]
    p2 <- dt$protein2[i]
    s <- dt$score[i]
    adj.m[p1, p2] <- s
    adj.m[p2, p1] <- s
  }

  adj.m
}

#' Build a Binary Adjacency Matrix from a STRING Interaction File
#'
#' @description 
#' This function processes a STRING database interaction file to create a symmetric binary 
#' adjacency matrix. It includes data cleaning (removing taxonomy prefixes), confidence 
#' threshold filtering, and optional ID mapping from Ensembl protein IDs to Gene IDs.
#'
#' @param file_path Character. The file path to the STRING interaction data (typically a .txt or .csv).
#' @param cutoff Numeric. The confidence score threshold for filtering interactions. 
#' Defaults to 700 (STRING "High Confidence" threshold).
#' @param convert_to_entrez Logical. Whether to convert Ensembl protein IDs (ENSP) 
#' to Ensembl Gene IDs (ENSG) using biomaRt. Defaults to TRUE.
#' @param species Character or Numeric. Taxonomy ID or name (e.g., 9606 or "human") 
#' to specify the biomaRt dataset. Supported: Human (9606) and Mouse (10090).
#'
#' @return A symmetric binary adjacency matrix (0/1) where rows and columns are 
#' gene/protein identifiers.
#' 
#' @details 
#' The function ensures the output is an unweighted, undirected graph representation 
#' by symmetrizing the matrix. If ID conversion is enabled, proteins that cannot be 
#' mapped to a gene ID will be removed.
#'
#' @import data.table
#' @import biomaRt
#' @export
string_to_adj_binary <- function(file_path,
                                 cutoff = 700,
                                 convert_to_entrez = TRUE,
                                 species = 9606) {
  # 1. 读取并清洗数据
  dt <- fread(file_path, header = TRUE)
  # 移除 ID 前缀 (如 "9606.")
  dt[, c("protein1", "protein2") := 
       lapply(.SD, sub, pattern = "^[0-9]+\\.", replacement = ""),
     .SDcols = c("protein1", "protein2")]

  # 2. 阈值过滤 (核心修改：只保留高置信度的边)
  dt <- dt[combined_score >= cutoff]

  if(nrow(dt) == 0) stop("No interactions found above the specified cutoff.")

  # 3. ID 转换
  if(convert_to_entrez) {
    mart_dataset <- switch(as.character(species),
                           "9606" = "hsapiens_gene_ensembl",
                           "human" = "hsapiens_gene_ensembl",
                           "10090" = "mmusculus_gene_ensembl",
                           "mouse" = "mmusculus_gene_ensembl",
                           stop("Unsupported species"))
    tryCatch({
      ensembl <- bi$useMart("ensembl", dataset = mart_dataset)
    }, error = function(e) {
      log_msg("ERROR","Failed to connect to Ensembl BioMart. Please check your internet connection.",quit=TRUE)
    })

    mapping <- bi$getBM(
      attributes = c("ensembl_peptide_id", "ensembl_gene_id"),
      filters = "ensembl_peptide_id",
      values = unique(c(dt$protein1, dt$protein2)),
      mart = ensembl
    )
    
    map_tbl <- setNames(mapping$ensembl_gene_id, mapping$ensembl_peptide_id)
    dt[, protein1 := map_tbl[protein1]]
    dt[, protein2 := map_tbl[protein2]]
    dt <- dt[!is.na(protein1) & !is.na(protein2)]
  }

  # 4. 构建 0/1 矩阵
  proteins <- unique(c(dt$protein1, dt$protein2))
  adj.m <- matrix(0, nrow = length(proteins), ncol = length(proteins),
                  dimnames = list(proteins, proteins))

  # 5. 填充矩阵
  idx_matrix <- as.matrix(dt[, .(protein1, protein2)])
  adj.m[idx_matrix] <- 1
  
  # 对称化 (无向图)
  idx_matrix_rev <- as.matrix(dt[, .(protein2, protein1)])
  adj.m[idx_matrix_rev] <- 1
  
  adj.m
}
