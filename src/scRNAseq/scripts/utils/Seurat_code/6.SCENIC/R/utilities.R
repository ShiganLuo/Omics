#' Generate output file path based on input arguments
#'
#' @param adata An object of class 'AnnData'.
#' @param colnames A character vector of column names.
#' @param outdir A character string indicating the output directory.
#' @param sample_name A character string indicating the sample name.
#' @param property A character string indicating the property of the output
#'   file. Defaults to NULL.
#' @param method A character string indicating the method used to generate the
#'   output file.
#' @param extension A character string indicating the file extension.
#'
#' @return A character string indicating the output file path.
#'
#' @examples
#' # Generate an output file path with default arguments
#' generate_output(
#'   adata = my_adata, colnames = "gene_expression", outdir = "/path/to/output",
#'   sample_name = "my_sample", method = "pca", extension = "csv"
#' )
#' # Generate an output file path with property argument
#' generate_output(
#'   adata = my_adata, colnames = "gene_expression", outdir = "/path/to/output",
#'   sample_name = "my_sample", property = "normalized", method = "pca",
#'   extension = "csv"
#' )
#'
#' @export
generate_output <-
  function(adata,
           colnames,
           outdir,
           sample_name,
           property = NULL,
           method,
           extension) {
    if (colnames %in% colnames(adata)) {
      if (is.null(property)) {
        output <-
          file.path(outdir,
                    paste(sample_name, method, extension, sep = "."))
      } else {
        output <-
          file.path(outdir,
                    paste(sample_name, property, method, extension, sep = "."))
      }
      return(output)
    }
  }
#' Extract a column from an AnnData object
#'
#' @param adata An object of class 'AnnData'.
#' @param colnames A character vector of column names.
#'
#' @return A vector with the values of the column.
#'
#' @examples
#' # Extract the gene expression values from an AnnData object
#' gene_expression <- extract_column(adata = my_adata, colnames =
#' "gene_expression")
#'
#' @export
extract_column <- function(adata, colnames) {
  if (colnames %in% colnames(adata)) {
    return(adata[[colnames]])
  }
}

#' Insert a column into an AnnData object
#'
#' @param adata An object of class 'AnnData'.
#' @param exist A character string indicating the name of an existing column to
#'   insert the new column after.
#' @param colnames A character vector of column names.
#' @param value A vector of values to be inserted as a column.
#'
#' @return An updated AnnData object with the new column inserted after the
#'   specified existing column.
#'
#' @examples
#' # Insert a new column of gene expression values after the 'batch' column
#' new_adata <- insert_column(adata = my_adata, exist = "batch", colnames =
#' "gene_expression", value = gene_exp_values)
#'
#' @export
insert_column <- function(adata, exist, colnames, value) {
  if (exist %in% colnames(adata)) {
    adata[[colnames]] <- value
    return(adata)
  } else {
    return(adata)
  }
}

#' Create parameter-value string
#'
#' This function creates a parameter-value string for a command-line interface
#' based on the input parameter and value. If the value is TRUE, only the
#' parameter is returned. If the value is not NULL and not a logical value, the
#' parameter and the value are returned as a character vector.
#'
#' @param parameter A character string representing the parameter.
#' @param value An optional value for the parameter. If the value is TRUE, only
#'   the parameter is returned. If the value is not NULL and not a logical
#'   value, the parameter and the value are returned as a character vector.
#'
#' @return If the value is TRUE, returns the parameter as a character string. If
#'   the value is not NULL and not a logical value, returns the parameter and
#'   value as a character vector.
#'
#' @examples
#' pv("--input", "data.txt")
#' pv("--debug", TRUE)
#' pv("--output")
#'
#' @export
pv <- function(parameter, value = NULL) {
  if (!is.null(value) && value == TRUE) {
    return(parameter)
  }
  if (!is.null(value) && !is.logical(value)) {
    return(c(parameter,
             value))
  }
}

#' Extract filename without extension from a file path
#'
#' @param x A character string indicating a file path.
#'
#' @return A character string indicating the filename without extension.
#'
#' @examples
#' # Extract filename without extension from a file path
#' cat_filename("/path/to/file.txt")
#'
#' @export
cat_filename <- function(x) {
  x <-
    stringr::str_split_fixed(basename(x), "\\.", 2)[, 1]
  return(x)
}

#' Set a directory path and create it if it does not exist
#'
#' @param path A character string indicating a directory path.
#'
#' @return A character string indicating the input directory path.
#'
#' @examples
#' # Set a directory path and create it if it does not exist
#' hp_set_path("/path/to/directory")
#'
#' @export
hp_set_path <- function(path) {
  if (!file.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  return(path)
}

#' Check if files exist and optionally stop if not
#'
#' @param x A character vector indicating file paths.
#' @param stop A logical value indicating whether to stop execution if a file
#'   does not exist. Defaults to TRUE.
#'
#' @return If stop is TRUE, the function stops execution if any file listed in x
#'   does not exist. If stop is FALSE, the function returns a character vector
#'   indicating the files that do not exist.
#'
#' @examples
#' # Check if files exist and stop if not hp_check_file(x =
#' c("/path/to/file1.txt", "/path/to/file2.txt"))
#'
#' # Check if files exist and return a character vector indicating the files
#' that do not exist hp_check_file(x = c("/path/to/file1.txt",
#' "/path/to/file2.txt"), stop = FALSE)
#'
#' @export
hp_check_file <- function(x, stop = TRUE) {
  not_exist <- x[!file.exists(x)]
  if (length(not_exist) > 0) {
    if (stop) {
      stop(
        paste0(
          "Please check the file path, the file listed below does not exist:\n"
        ),
        paste0(not_exist, collapse = "\n")
      )
    } else {
      return(not_exist)
    }
  }
}

#' Write data to an Excel file with customization options
#'
#' This function writes data in the form of a data frame to an Excel file with
#' optional customization options.
#'
#' @param x A data frame to be written to an Excel file.
#' @param filename A string specifying the name of the output Excel file.
#' @param sheet_names A vector of strings specifying the names of the sheets in
#'   the Excel file. Default is NULL.
#' @return None
#'
#' @import WriteXLS
#'
#' @examples
#' hp_write_xlsx(iris, "iris.xlsx", sheet_names = c("Sheet1", "Sheet2"))
#'
#' @export
hp_write_xlsx <- function(x, filename, sheet_names = NULL) {
  WriteXLS(
    x = x,
    ExcelFileName = filename,
    SheetNames = sheet_names,
    AdjWidth = TRUE,
    row.names = FALSE,
    FreezeRow = 1,
    BoldHeaderRow = TRUE,
    AutoFilter = TRUE
  )
}

#' Perform differential expression analysis and annotate genes
#'
#' This function performs differential expression analysis on a Seurat object
#' based on specified cell types, groups, treatments, and controls. The function
#' outputs a data frame of differentially expressed genes with optional gene
#' annotation and Excel file writing.
#'
#' @param x A Seurat object.
#' @param cell_type A string specifying the name of the column in the metadata
#'   slot of \code{x} that contains the cell type information.
#' @param group A string specifying the name of the column in the metadata slot
#'   of \code{x} that contains the group information.
#' @param treatment A string specifying the treatment group for differential
#'   expression analysis.
#' @param control A string specifying the control group for differential
#'   expression analysis.
#' @param features A vector of gene names to be used for differential expression
#'   analysis. Default is NULL, which uses all genes in \code{x}.
#' @param test_use A string specifying the statistical test to use for
#'   differential expression analysis. Default is "wilcox".
#' @param logfc_threshold A numeric value specifying the threshold for log2
#'   fold-change to consider a gene differentially expressed. Default is 0.25.
#' @param p_value_threshold A numeric value specifying the threshold for p-value
#'   to consider a gene differentially expressed. Default is 0.05.
#' @param min_pct A numeric value specifying the minimum percentage of cells
#'   expressing a gene to consider it for differential expression analysis.
#'   Default is 0.1.
#' @param verbose A logical value specifying whether to print verbose output
#'   during differential expression analysis. Default is FALSE.
#' @param annotating A string specifying the gene annotation source to use for
#'   the output data frame. Default is NULL, which does not annotate the genes.
#' @param outdir A string specifying the directory to write the output Excel
#'   file to. Default is NULL, which does not write an output Excel file.
#' @return A data frame of differentially expressed genes with gene annotation
#'   and an optional Excel file output.
#'
#' @importFrom Seurat Idents FindMarkers
#' @importFrom dplyr group_by summarise desc inner_join mutate inner_join
#'   arrange if_else
#' @importFrom tibble rownames_to_column
#' @examples
#' cat_deg(
#'   x = pbmc_small,
#'   cell_type = "cell_type",
#'   group = "stim",
#'   treatment = "stimulated",
#'   control = "unstimulated",
#'   features = c("CD3D", "CD4", "CD8A", "CD19", "CD14", "FCGR3A"),
#'   test_use = "t",
#'   logfc_threshold = 0.5,
#'   p_value_threshold = 0.01,
#'   min_pct = 0.1,
#'   verbose = TRUE,
#'   annotating = "human",
#'   outdir = "deg_output"
#' )
#'
#' @export
hp_calculate_difference <-
  function(x,
           treatment,
           control,
           group = NULL,
           cluster = NULL,
           features = NULL,
           test_use = "wilcox",
           logfc_threshold = 0.25,
           p_value_threshold = 0.05,
           min_pct = 0.1,
           min_unit_group = 3,
           verbose = FALSE,
           annotating = NULL,
           outdir = NULL,
           ...) {
    if (!is.null(annotating) &&
        !requireNamespace("annotables", quietly = TRUE)) {
      stop(
        "The 'lisi' package is not installed.
           Please install it using
           devtools::install_github('stephenturner/annotables')."
      )
    }

    SeuratObject::Idents(x) <- cluster

    freq <-
      rowSums(table(x@meta.data[[cluster]],
                    x@meta.data[[group]]) > min_unit_group)
    clusters <- names(freq[freq == 2])

    logger::log_info("Differentially expressed gene analysis...")
    adata <- lapply(clusters, function(cluster) {
      tic(cluster)
      deg <- FindMarkers(
        x,
        features = features,
        subset.ident = cluster,
        ident.1 = treatment,
        ident.2 = control,
        group.by = group,
        test.use = test_use,
        logfc.threshold = logfc_threshold,
        verbose = verbose,
        min.pct = min_pct,
        min.cells.group = min_unit_group,
        ...
      ) |>
        rownames_to_column("gene") |>
        mutate(
          cluster = cluster,
          group = group,
          case = paste(treatment, "vs",
                       control,
                       sep = "_"),
          test_use = test_use
        )
      toc()
      return(deg)
    })
    adata <- adata |> Reduce(f = rbind)
    adata <- adata |>
      mutate(change = if_else(
        abs(avg_log2FC) >= logfc_threshold & p_val <= p_value_threshold,
        if_else(
          p_val <= p_value_threshold &
            avg_log2FC >= logfc_threshold,
          "Upregulated",
          "Downregulated"
        ),
        "Stable"
      )) |>
      arrange(cluster, desc(change), desc(avg_log2FC))
    if (!is.null(annotating)) {
      if (annotating == "human") {
        annotation_tables <- annotables::grch38
      } else if (annotating == "mouse") {
        annotation_tables <- annotables::grcm38
      }
      adata <- adata |>
        inner_join(annotation_tables, by = c("gene" = "symbol"))
    }
    if (!is.null(outdir)) {
      outdir <- hp_set_path(outdir)
      hp_write_xlsx(
        x = adata,
        filename = file.path(outdir,
                             paste("deg",
                                   test_use,
                                   group,
                                   "xlsx",
                                   sep = ".")),
        sheet_names = "deg"
      )
    }
    print(table(adata[["cluster"]], adata[["change"]]))
    return(adata)
  }

#' Differential Gene Expression Analysis using DESeq2
#'
#' This function performs differential gene expression analysis using the DESeq2
#' package.
#'
#' @param x A Seurat object containing RNA-seq data
#' @param control A character string specifying the name of the control group
#' @param treatment A character string specifying the name of the treatment
#'   group
#' @param log2fc_threshold A numeric value specifying the threshold for log2
#'   fold change
#' @param p_value_threshold A numeric value specifying the threshold for p-value
#' @param outdir A character string specifying the output directory
#'
#' @return A data frame containing differentially expressed genes and their fold
#'   change values and statistical significance
#'
#' @importFrom dplyr mutate pull
#' @importFrom tibble rownames_to_column
#' @importFrom readr write_csv
#'
#' @examples
#' cat_deseq2(seurat_object, "Control", "Treatment")
#'
#' @export
cat_deseq2 <- function(x,
                       control,
                       treatment,
                       log2fc_threshold = 0.5,
                       p_value_threshold = 0.05,
                       outdir = "./deseq2") {
  outdir <- hp_set_path(file.path(x@board, "deseq2"))
  x <- x |>
    subset(group %in% c(control, treatment))
  counts <- as.data.frame(x@assays$RNA@counts)
  metadata <- x@meta.data |>
    mutate(group = factor(group))

  if (all(pull(metadata, sample) %in% colnames(counts))) {
    counts <- counts[rowSums(counts) > 0,]
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                          colData = metadata,
                                          design = ~ group)
    dds <- DESeq2::DESeq(dds, parallel = TRUE)
    saveRDS(dds, file.path(outdir, paste0(
      treatment,
      "_vs_",
      control, "_dds.rds"
    )))
    # 差异基因
    res <-
      DESeq2::results(dds, contrast = c("group", treatment, control))

    res_ordered <- res[order(res$pvalue),]
    deg <- na.omit(as.data.frame(res_ordered))

    deg <- deg |>
      mutate(change = ifelse(
        pvalue <= p_value_threshold &
          abs(log2FoldChange) >= log2fc_threshold,
        ifelse(
          log2FoldChange > log2fc_threshold,
          "Upregulated",
          "Downregulated"
        ),
        "Stable"
      )) |>
      rownames_to_column("gene")

    write_csv(deg, file.path(outdir, paste0(
      treatment,
      "_vs_",
      control, "_deg.csv"
    )))
  }
  return(deg)
}

#' Gene Enrichment Analysis using clusterProfiler
#'
#' This function performs gene enrichment analysis using the clusterProfiler
#' package.
#'
#' @param adata A data frame containing gene expression data
#' @param gene_clusters A character vector specifying the gene clusters to be
#'   analyzed
#' @param databases A character string specifying the database to be used,
#'   either "GOBP" or "KEGG"
#' @param species A character string specifying the species, either "human" or
#'   "mouse"
#' @param prefix A character string specifying the prefix to be added to the
#'   output file names
#' @param outdir A character string specifying the output directory
#'
#' @return A clusterProfiler object containing enriched terms and their
#'   statistical significance
#'
#' @importFrom dplyr full_join rename select
#'
#' @examples
#' hp_enrich(adata, "cluster1", "GOBP", "human", "prefix", "outdir")
#'
#' @export
hp_enrich <- function(adata,
                      gene_clusters,
                      databases = NULL,
                      species = "human",
                      prefix = NULL,
                      outdir = NULL) {
  adata <- as.data.frame(adata)
  # TODO: require install org.Hs.eg.db or org.Mm.eg.db
  # species
  org_db <- switch(species,
                   human = "org.Hs.eg.db",
                   mouse = "org.Mm.eg.db")
  organism <- switch(species,
                     human = "hsa",
                     mouse = "mmu")

  log_info("{databases} database gene enrichment analysis...")
  # GOBP
  if (databases == "GOBP") {
    x <- clusterProfiler::compareCluster(
      geneClusters = gene_clusters,
      data = adata,
      fun = "enrichGO",
      OrgDb = org_db,
      pvalueCutoff = 0.05,
      keyType = "SYMBOL",
      ont = "BP"
    )
  }

  # KEGG
  if (databases == "KEGG") {
    gene <- all.vars(gene_clusters)[1]
    gid <-
      clusterProfiler::bitr(unique(adata[, gene]), "SYMBOL", "ENTREZID",
                            OrgDb = org_db)
    adata <-
      dplyr::full_join(adata, gid, by = c(gene = "SYMBOL")) |>
      dplyr::select(-gene) |>
      dplyr::rename(gene = ENTREZID)
    x <-
      clusterProfiler::compareCluster(
        gene_clusters,
        data = adata,
        fun = "enrichKEGG",
        pvalueCutoff = 1,
        organism = organism
      )
    x <- clusterProfiler::setReadable(x,
                                      OrgDb = org_db,
                                      keyType = "ENTREZID")
  }
  if (!is.null(outdir)) {
    outdir <- hp_set_path(outdir)
    if (!is.null(prefix)) {
      prefix <- paste0(".", prefix)
    }
    saveRDS(x, file.path(outdir,
                         paste0(
                           "enrichment.", databases, prefix, ".rds"
                         )))
    hp_write_xlsx(
      x = x@compareClusterResult,
      filename = file.path(outdir,
                           paste0(
                             "enrichment.", databases, prefix, ".xlsx"
                           )),
      sheet_names = databases
    )
  }
  return(x)
}


#' Check if one or more packages are installed
#'
#' This function checks whether one or more packages are installed in the
#' current R environment, and raises an error if any of the packages is not
#' found.
#'
#' @param ... One or more character strings specifying the names of the packages
#'   to check.
#' @param error Logical indicating whether to raise an error if any of the
#'   packages is not found (default is TRUE).
#'
#' @return A logical vector indicating whether each package is installed or not.
#'
#' @examples
#' # Check if the "ggplot2" package is installed
#' hp_check_package("ggplot2")
#'
#' # Check if multiple packages are installed
#' hp_check_package("ggplot2", "dplyr", "tidyr")
#'
#' @export
hp_check_package <- function(..., error = TRUE) {
  pkgs <- unlist(list(...))
  package_installed <-
    suppressWarnings(sapply(pkgs, requireNamespace))
  if (error && any(!package_installed)) {
    missing_pkgs <- pkgs[!package_installed]
    stop(paste(
      "Cannot find the following packages:",
      paste(missing_pkgs, collapse = ", ")
    ),
    ". Please install.")
  }
  invisible(package_installed)
}


#' Create KEGG gene sets data frame
#'
#' This function creates a data frame of KEGG gene sets from a JSON file.
#'
#' @param json_path The path to the JSON file containing the KEGG gene sets.
#'
#' @return A data frame with the KEGG gene sets, including the levels and genes.
#'
#' @importFrom rjson fromJSON
#' @importFrom stringr str_sub str_replace_all
#'
#' @examples
#' # Download KEGG json file from:
#' https://www.genome.jp/kegg-bin/show_organism?org=hsa
#' # Create KEGG gene sets data frame
#' hp_create_kegg("hsa00001.json")
#'
#' @export
hp_create_kegg <- function(json_path) {
  result <- fromJSON(file = json_path)

  # Metabolism
  gene_sets <- data.frame()
  for (i in 1:8) {
    # 1111
    # i <- 1
    # print("===================")
    # print(result[["children"]][[i]]$name)
    for (a1 in 1:length(result[["children"]][[i]]$children)) {
      # print(result[["children"]][[i]]$children[[a1]]$name)
      for (a2 in 1:length(result[["children"]][[i]]$children[[a1]]$children)) {
        if (length(
          result[["children"]][[i]]$children[[a1]]$children[[a2]]) == 2
          ) {
          # level4 <-
          #   result[["children"]][[i]]$children[[a1]]$children[[a2]]$name
          genes <-
            lapply(
              result[["children"]][[i]]$children[[a1]]$children[[a2]]$children,
              function(x) {
              strsplit(strsplit(x$name, ";")[[1]][1], " ")[[1]][2]
            })
          genes <- unlist(genes)
          level1 <- result[["children"]][[i]]$name
          level2 <- result[["children"]][[i]]$children[[a1]]$name
          level3 <-
            result[["children"]][[i]]$children[[a1]]$children[[a2]]$name
          # print(level3)
          temp_gene_sets <- data.frame(
            level1 = level1,
            level2 = level2,
            level3 = level3,
            # level4 = level4,
            genes = genes
          )
          gene_sets <- rbind(gene_sets, temp_gene_sets)
        }
      }
    }
  }
  gene_sets$level1 <- str_sub(gene_sets$level1, 7, -1)
  gene_sets$level2 <- str_sub(gene_sets$level2, 7, -1)
  gene_sets$level3 <- str_sub(gene_sets$level3, 7, -1)
  gene_sets$level3 <-
    str_replace_all(gene_sets$level3, " \\[.+?\\]", "")
  return(gene_sets)
}

#' Calculate correlation between two features in a data set
#'
#' This function calculates the correlation between two features (columns) in a
#' data set using a specified correlation method.
#'
#' @param x A data frame or matrix containing the data set.
#' @param feature_x The name or index of the first feature (column) to calculate
#'   the correlation for.
#' @param feature_y The name or index of the second feature (column) to
#'   calculate the correlation for.
#' @param method The correlation method to use. Defaults to "pearson".
#' @param orientation The orientation of the data set. If "column", the
#'   correlation is calculated between two columns. If "row", the correlation is
#'   calculated between two rows.
#'
#' @return A tibble with the following columns:
#' \describe{
#'   \item{feature_x}{The name or index of the first feature.}
#'   \item{feature_y}{The name or index of the second feature.}
#'   \item{p_value}{The p-value of the correlation test.}
#'   \item{estimate}{The estimate of the correlation coefficient.}
#'   \item{n}{The number of observations used in the correlation test.}
#'   \item{method}{The correlation method used.}
#' }
#'
#' @examples
#' data(mtcars)
#' .calculate_correlation(mtcars, "mpg", "wt")
.calculate_correlation <- function(x,
                                   feature_x,
                                   feature_y,
                                   method = "pearson",
                                   orientation = "column") {
  if (orientation == "column") {
    cor_test_res <-
      cor.test(x[[feature_x]], x[[feature_y]], method = method)
    p_value <- cor_test_res$p.value
    estimate <- cor_test_res$estimate
    n <- nrow(x)
  } else if (orientation == "row") {
    cor_test_res <-
      cor.test(x[feature_x,], x[feature_y,], method = method)
    p_value <- cor_test_res$p.value
    estimate <- cor_test_res$estimate
    n <- ncol(x)
  } else {
    stop("Invalid orientation parameter. Use 'column' or 'row'.")
  }

  return(
    tibble(
      feature_x = feature_x,
      feature_y = feature_y,
      p_value = p_value,
      estimate = estimate,
      n = n,
      method = method
    )
  )
}

#' Calculate pairwise correlations between features in a Seurat object
#'
#' This function calculates pairwise correlations between features in a Seurat
#' object for a given assay.
#'
#' @param x A Seurat object or a data.frame containing the expression data.
#' @param feature_x A character vector of feature names to use as the x-axis in
#'   the correlation calculation.
#' @param feature_y A character vector of feature names to use as the y-axis in
#'   the correlation calculation. If NULL, the same feature names as feature_x
#'   will be used.
#' @param method The correlation coefficient to use. Must be one of "pearson",
#'   "spearman", or "kendall".
#' @param assay The assay to use for the correlation calculation. Default is
#'   "RNA".
#' @param slot The slot in the Seurat object where the assay data is stored.
#'   Default is "data".
#' @param ncores The number of CPU cores to use for parallel computation.
#'   Default is 4.
#'
#' @return A data.frame containing a pairwise comparison of each feature in
#'   feature_x and feature_y, with the correlation coefficient and p-value for
#'   each comparison.
#'
#' @examples
#' data("pbmc")
#' hp_calculate_correlation(x = pbmc,
#' feature_x = c("CD3D", "CD4", "CD8A", "CD8B"),
#' method = "pearson")
#'
#' @importFrom Seurat GetAssayData
#' @importFrom Matrix as.matrix
#' @importFrom furrr future_map_dfr
#' @importFrom progressr progressor with_progress
#' @importFrom future plan multisession sequential
#' @export
hp_calculate_correlation <- function(x,
                                     feature_x,
                                     feature_y = NULL,
                                     method = "pearson",
                                     asssay = "RNA",
                                     slot = "data",
                                     ncores = 1,
                                     parallel_method = "furrr",
                                     ...) {
  # Simplify input parameters
  if (is.null(feature_y)) {
    feature_y <- feature_x
  }

  # Extract data matrix
  if (inherits(x, what = "Seurat")) {
    data_matrix <-
      GetAssayData(x, slot = slot, assay = asssay)
  } else {
    data_matrix <- x
  }
  data_matrix <- as.matrix(data_matrix)

  # Generate pairwise combinations of features
  feature_combinations <-
    expand.grid(feature_x = feature_x, feature_y = feature_y)
  if (parallel_method == "plyr") {
    print("Use plyr...")
    x <- feature_combinations[["feature_x"]]
    y <- feature_combinations[["feature_y"]]
    correlation_results <-
      plyr::ldply(
        .data = seq_along(x),
        .fun = function(i) {
          .calculate_correlation(
            x = data_matrix,
            feature_x = x[i],
            feature_y = x[i],
            method = method,
            orientation = "row"
          )
        },
        ...
      )
  } else if (parallel_method == "furrr") {
    print("Use furrr...")
    if (ncores == 1) {
      plan(sequential)
    } else {
      plan(multisession, workers = ncores)
    }
    with_progress({
      p <- progressor(steps = nrow(feature_combinations))
      correlation_results <-
        future_pmap_dfr(
          .l = feature_combinations,
          .f = function(feature_x, feature_y) {
            p()
            .calculate_correlation(
              x = data_matrix,
              feature_x = feature_x,
              feature_y = feature_y,
              method = method,
              orientation = "row"
            )
          },
          ...
        )
    })
  }

  return(correlation_results)
}



ff <- function(feature_combinations, data_matrix, method, ...) {
  p <- progressor(steps = nrow(feature_combinations))
  correlation_results <-
    future_pmap_dfr(
      .l = feature_combinations,
      .f = function(feature_x, feature_y) {
        p()
        .calculate_correlation(
          x = data_matrix,
          feature_x = feature_x,
          feature_y = feature_y,
          method = method,
          orientation = "row"
        )
      },
      ...
    )
}
