#' @include utilities.R
NULL

#' Run CellphoneDB analysis on normalized gene expression data
#'
#' This function runs CellphoneDB analysis on normalized gene expression data for either human or mouse species. It converts mouse gene names to human gene names using biomaRt package. It creates the necessary input files for CellphoneDB analysis and runs the analysis using the specified number of cores.
#'
#' @param x A Seurat object containing normalized gene expression data
#' @param outdir Directory to save CellphoneDB output files
#' @param prefix Prefix for CellphoneDB output files
#' @param cellphonedb Path to CellphoneDB executable. If NULL, default path "/opt/cellphonedb/3.1.0/bin/cellphonedb" is used.
#' @param species Species for analysis. Either "human" or "mouse". Default is "human".
#' @param counts_data Type of gene identifier used in the counts data. Default is "hgnc_symbol".
#' @param ncores Number of cores to use for the analysis. Default is 10.
#' @param host Host for biomaRt package. Default is "https://www.ensembl.org".
#'
#' @return NULL
#'
#' @importFrom data.table fwrite
#' @importFrom dplyr select
#' @importFrom logger log_info
#' @importFrom Seurat Idents NormalizeData
#' @importFrom tictoc tic toc
#'
#' @examples
#' hp_run_cellphonedb(x = seurat_obj, outdir = "./cellphonedb_output")
#'
#' @export
hp_run_cellphonedb <- function(x,
                               outdir = "./",
                               prefix = NULL,
                               cellphonedb = NULL,
                               species = "human",
                               counts_data = "hgnc_symbol",
                               ncores = 10,
                               host = "https://www.ensembl.org") {
  if (is.null(cellphonedb)) {
    cellphonedb <- "/opt/cellphonedb/3.1.0/bin/cellphonedb"
    hp_check_file(cellphonedb)
  }
  outdir <- hp_set_path(outdir)
  x <- NormalizeData(x, assay = "RNA")
  normalize_data <- x@assays$RNA@data
  if (species == "mouse") {
    human <-
      biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host)
    mouse <-
      biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host)

    # Basic function to convert mouse to human gene names
    mouse_genes <- rownames(normalize_data)
    genesV2 <- biomaRt::getLDS(
      attributes = c("mgi_symbol"),
      filters = "mgi_symbol",
      values = mouse_genes,
      mart = mouse,
      attributesL = c("hgnc_symbol"),
      martL = human,
      uniqueRows = T
    )
    print(head(genesV2))
    sp1_counts <-
      data.frame(gene = rownames(normalize_data),
                 normalize_data,
                 check.names = F)
    # 转小鼠基因名为人类基因名
    sp1_counts$Gene <-
      genesV2[match(sp1_counts$gene, genesV2[, 1]), 2]
    sp1_counts <- subset(sp1_counts, Gene != 'NA')
    sp1_counts <- dplyr::select(sp1_counts, Gene, everything())
    sp1_counts <- sp1_counts[, !(colnames(sp1_counts) %in% 'gene')]
    dim(sp1_counts)
    sp1_counts[1:2, 1:2]
    normalize_data <- sp1_counts
  }
  tictoc::tic("Run CellphoneDB")
  logger::log_info("Run CellphoneDB...")
  normalize_data_path <- file.path(outdir, "data.txt")
  metadata_path <- file.path(outdir, "metadata.tsv")
  if (!file.exists(metadata_path)) {
    data.table::fwrite(normalize_data, normalize_data_path, sep = "\t")
    # normalize_data <- x@assays$RNA@data
    # writeMM(normalize_data, file = file.path(normalize_data_path,
    #                                          "matrix.mtx"))
    # write(
    #   x = rownames(normalize_data),
    #   file = file.path(normalize_data_path,
    #                    "features.tsv")
    # )
    # write(
    #   x = colnames(normalize_data),
    #   file = file.path(normalize_data_path,
    #                    "barcodes.tsv")
    # )
    if (!("cell_type" %in% colnames(x@meta.data))) {
      x$cell_type <- Idents(x)
    }
    x@meta.data$cell_id <- rownames(x@meta.data)
    metadata <- x@meta.data[, c("cell_id", "cell_type")]
    write.table(
      metadata,
      file = metadata_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  significant_means <- file.path(outdir, "significant_means.txt")
  if (!file.exists(significant_means)) {
    system2(
      command = cellphonedb,
      args = c(
        "method",
        "statistical_analysis",
        metadata_path,
        normalize_data_path,
        "--output-path",
        outdir,
        "--counts-data",
        counts_data,
        "--threads",
        ncores
      )
    )
  }
  tictoc::toc()
}

#' Automated cell type annotation for scRNA-seq datasets
#'
#' @param x The seurat object contains counts
#' @param celltypist The absolute path of the celltypist
#' @param model Model used for predictions. If not provided,
#' default to using the `Immune_All_High.pkl` model.
#' Detailed model information can be found at
#' `https://www.celltypist.org/models`.
#' @param outdir Directory to store the output files.
#'
#' @return Seurat object
#' @importFrom tibble column_to_rownames
#' @importFrom Seurat AddMetaData
#' @importFrom data.table fwrite
#' @importFrom tictoc tic toc
#' @importFrom logger log_info
#' @importFrom tibble column_to_rownames
#' @export
#'
#' @examples
#' \dontrun{
#' library(hephaestus)
#' data("pbmc")
#' pbmc <- hp_run_celltypist(
#'   x = pbmc,
#'   celltypist = "/opt/celltypist/1.3.0/bin/celltypist"
#' )
#' DimPlot(
#'   pbmc,
#'   group.by = "Immune_All_High_predicted_labels",
#'   label = TRUE,
#'   repel = TRUE
#' )
#' }
hp_run_celltypist <-
  function(x,
           celltypist = NULL,
           model = "Immune_All_High.pkl",
           outdir = NULL) {
    celltypist <- celltypist %||% "/opt/celltypist/1.3.0/bin/celltypist"
    hp_check_file(celltypist)

    tic("Run CellTypist")
    log_info("Run CellTypist...")
    # celltypist --update-models
    if (is.null(outdir)) {
      outdir <- tempdir()
    } else {
      outdir <- hp_set_path(outdir)
    }
    if (inherits(x, what = "Seurat")) {
      input_data <- file.path(outdir, "counts.csv")
      if (!file.exists(input_data)) {
        counts <- x@assays$RNA@counts |>
          as.data.frame()
        # write.csv(counts, input_data, quote = FALSE)
        fwrite(
          x = counts,
          file = input_data,
          quote = FALSE,
          row.names = TRUE
        )
      }
    } else {
      input_data <- x
    }
    predicted_labels_path <-
      file.path(outdir, paste0(model, ".", "predicted_labels.csv"))
    if (!file.exists(predicted_labels_path)) {
      system2(
        command = celltypist,
        args = c(
          "--majority-voting",
          "--indata",
          input_data,
          "--model",
          model,
          "--outdir",
          outdir,
          "--transpose-input",
          "--prefix",
          paste0(model, ".")
        )
      )
      predicted_labels <- read.csv(predicted_labels_path)
    } else {
      predicted_labels <- read.csv(predicted_labels_path)
    }
    toc()
    metadata <- column_to_rownames(predicted_labels, "X")
    colnames(metadata) <-
      c(
        paste0(gsub(".pkl", "", model), "_predicted_labels"),
        paste0(gsub(".pkl", "", model), "_over_clustering"),
        paste0(gsub(".pkl", "", model), "_majority_voting")
      )
    x <-
      AddMetaData(x, metadata = metadata)
    return(x)
  }

#' A lightning-fast python implementation of the SCENIC pipeline (Single-Cell
#' rEgulatory Network Inference and Clustering) which enables biologists to
#' infer transcription factors, gene regulatory networks and cell types from
#' single-cell RNA-seq data
#'
#' @param x
#' @param outdir
#' @param pyscenic
#' @param species
#' @param tfs_fname
#' @param database_fname
#' @param annotations_fname
#' @param mode
#' @param ncores
#' @param seed
#'
#' @return
#' @importFrom data.table fwrite
#' @importFrom logger log_info
#' @importFrom tictoc tic toc
#' @importFrom stringr str_split_fixed
#' @importFrom dplyr select rename count mutate left_join
#' @export
#'
#' @examples
hp_run_pyscenic <-
  function(x,
           outdir = "./pyscenic",
           pyscenic = NULL,
           species = "human",
           features = NULL,
           tfs_fname = NULL,
           database_fname = NULL,
           annotations_fname = NULL,
           mode = "custom_multiprocessing",
           ncores = 10,
           seed = 717) {
    if (is.null(pyscenic)) {
      pyscenic <- "/opt/pyscenic/0.12.1/bin/pyscenic"
      hp_check_file(pyscenic)
    }
    if (species == "human") {
      tfs_fname <- "/DATA/public/cistarget/tf_lists/allTFs_hg38.txt"
      database_fname <-
        c(
          "/DATA/public/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
          "/DATA/public/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
        )
      annotations_fname <-
        "/DATA/public/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    } else if (species == "mouse") {
      tfs_fname <- "/DATA/public/cistarget/tf_lists/allTFs_mm.txt"
      database_fname <-
        c(
          "/DATA/public/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
          "/DATA/public/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
        )
      annotations_fname <-
        "/DATA/public/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl"
    } else {
      stop("Please check if the species is human or mouse!")
    }
    hp_check_file(x = c(tfs_fname, database_fname, annotations_fname))
    outdir <- hp_set_path(outdir)
    expression_mtx_fname <- file.path(outdir, "expression_mtx.csv")
    module_fname <-
      file.path(outdir, "expression_mtx.adjacencies.tsv")
    signatures_fname <- file.path(outdir, "regulons.gmt")
    auc_mtx_fname <- file.path(outdir, "auc_mtx.csv")
    tfs_target_fname <- file.path(outdir, "tfs_targer.tsv")
    auc_g_mtx_fname <- file.path(outdir, "auc_g_mtx.csv")
    # pre-processing
    if (!is.null(features)) {
      # x <- x[features,]
      tfs <- read.table(tfs_fname)[, 1]
      filtered_tfs <- intersect(tfs, features)
      tfs_fname <- file.path(outdir, "tfs.txt")
      write.table(
        x = filtered_tfs,
        file = tfs_fname,
        row.names = FALSE,
        quote = FALSE,
        col.names = FALSE
      )
      counts <- x@assays$RNA@counts[features,] |>
        as.data.frame()
    } else {
      counts <- x@assays$RNA@counts |>
        as.data.frame()
    }
    if (!file.exists(expression_mtx_fname)) {
      log_info("Write counts to a csv file...")
      tic("Write csv")
      fwrite(
        x = counts,
        file = expression_mtx_fname,
        quote = FALSE,
        row.names = TRUE
      )
      toc()
    }
    # grn
    if (!file.exists(module_fname)) {
      log_info("Derive co-expression modules from expression matrix...")
      tic("Run grn")
      system2(
        command = pyscenic,
        args = c(
          "grn",
          "--transpose",
          "--num_workers",
          ncores,
          "--seed",
          seed,
          "--output",
          module_fname,
          expression_mtx_fname,
          tfs_fname
        )
      )
      toc()
    }
    # ctx
    if (!file.exists(signatures_fname)) {
      log_info(
        "Find enriched motifs for a gene signature and optionally prune targets from this signature based on cis-regulatory cues..."
      )
      tic("Run ctx")
      system2(
        command = pyscenic,
        args = c(
          "ctx",
          "--transpose",
          "--output",
          signatures_fname,
          "--mode",
          mode,
          "--annotations_fname",
          annotations_fname,
          "--num_workers",
          ncores,
          "--expression_mtx_fname",
          expression_mtx_fname,
          module_fname,
          database_fname
        )
      )
      toc()
    }
    # aucell
    if (!file.exists(auc_mtx_fname)) {
      log_info("Quantify activity of gene signatures across single cells...")
      tic("Run aucell")
      system2(
        command = pyscenic,
        args = c(
          "aucell",
          "--transpose",
          "--output",
          auc_mtx_fname,
          "--num_workers",
          ncores,
          expression_mtx_fname,
          signatures_fname
        )
      )
      toc()
    }

    # tfs_target
    if (!file.exists(tfs_target_fname)) {
      signatures <-
        clusterProfiler::read.gmt(signatures_fname)
      signatures_count <- signatures |>
        count(term)
      signatures <-
        left_join(x = signatures, y = signatures_count,
                  by = "term")
      signatures <- signatures |>
        mutate(term = str_split_fixed(
          string = term,
          pattern = "\\(",
          n = 2
        )[, 1]) |>
        mutate(symbol = paste0(term, " (", n, "g)")) |>
        select(symbol, term, gene) |>
        rename(target_gene = gene, tf = term)
      write.csv(signatures,
                tfs_target_fname,
                row.names = FALSE)
    } else {
      signatures <- read.csv(tfs_target_fname)
    }
    # auc_g_mtx
    if (!file.exists(auc_g_mtx_fname)) {
      auc_mtx <-
        read.csv(auc_mtx_fname,
                 row.names = 1,
                 check.names = FALSE)
      rownames(auc_mtx) <- unique(signatures$tf)
      fwrite(
        x = auc_mtx,
        file = auc_g_mtx_fname,
        quote = FALSE,
        row.names = TRUE
      )
    } else {
      auc_mtx <- read.csv(auc_g_mtx_fname,
                          row.names = 1,
                          check.names = FALSE)
    }
    colnames(auc_mtx) <-
      gsub(pattern = "\\.",
           replacement = "-",
           x = colnames(auc_mtx))
    x[["scenic"]] <- CreateAssayObject(counts = auc_mtx)

    return(x)
  }


#' Run Cell Ranger count pipeline
#'
#' This function runs the Cell Ranger count pipeline for single-cell RNA sequencing data analysis.
#'
#' @param fastqs A character vector specifying the path to the fastq files
#' @param sample A character string specifying the sample name
#' @param species A character string specifying the species, either "human" or "mouse"
#' @param cellranger A character string specifying the path to the Cell Ranger executable, if not specified, "/opt/cellranger/7.1.0/bin/cellranger" will be used
#' @param index A character string specifying the path to the reference transcriptome index, if not specified, the default index for the corresponding species will be used
#' @param memory An integer specifying the amount of memory to be used (in GB)
#' @param ncores An integer specifying the number of CPU cores to be used
#' @param outdir A character string specifying the output directory
#'
#' @examples
#' hp_run_cellranger(
#'   fastqs = "/path/to/fastqs",
#'   sample = "sample_name",
#'   species = "human",
#'   cellranger = "/path/to/cellranger",
#'   index = "/path/to/reference/transcriptome/index",
#'   memory = 60,
#'   ncores = 20,
#'   outdir = "./cellranger_count"
#' )
#'
#' @export
hp_run_cellranger <-
  function (fastqs,
            sample = NULL,
            species = "human",
            cellranger = NULL,
            index = NULL,
            memory = 60,
            ncores = 20,
            outdir = "./cellranger_count",
            log = NULL)
  {
    if (is.null(cellranger)) {
      cellranger <- "/opt/cellranger/7.1.0/bin/cellranger"
    }
    if (is.null(index) && species == "human") {
      index <- "/DATA/public/references/refdata-gex-GRCh38-2020-A/"
    }
    else if (is.null(index) && species == "mouse") {
      index <- "/DATA/public/references/refdata-gex-mm10-2020-A/"
    }
    if (is.null(sample)) {
      sample <- basename(fastqs)
    }
    outdir <- hp_set_path(outdir)
    wd <- normalizePath(outdir)
    processx::run(
      command = cellranger,
      args = c(
        "count",
        "--id",
        sample,
        "--fastqs",
        normalizePath(fastqs),
        "--sample",
        sample,
        "--transcriptome",
        index,
        "--localmem",
        memory,
        "--localcores",
        ncores
      ),
      wd = wd,
      echo_cmd = TRUE,
      spinner = TRUE,
      stderr_to_stdout = TRUE,
      stdout = log
    )
  }


#' Calculate cell proportions for each cluster and sample
#'
#' This function calculates the proportion of cells in each cluster for each sample in a metadata table.
#'
#' @param x A metadata table or a Seurat object containing a metadata table.
#' @param cluster The name of the column in the metadata table that contains the cluster information.
#' @param sample The name of the column in the metadata table that contains the sample information.
#' @param additional_cols A character vector of additional column names to include in the output.
#' @return A data frame containing the cell proportions for each cluster, sample, and additional column.
#' @export
#'
#' @examples
#' # Load example metadata table
#' metadata <- read.csv("metadata.csv")
#'
#' # Calculate cell proportions
#' cell_props <- hp_calculate_cell_proportion(metadata, cluster = "cell_type", sample = "sample", additional_cols = c("treatment", "replicate"))
#'
#' # View the resulting data frame
#' head(cell_props)
#'
#' @import Seurat
#' @import dplyr
#' @import purrr
#' @importFrom rlang sym
hp_calculate_cell_proportion <-
  function(x,
           cluster = "cell_type",
           sample = "sample",
           additional_cols = NULL) {
    # Check if input is a Seurat object or a metadata table
    if (inherits(x, what = "Seurat")) {
      metadata <- x@meta.data
    } else {
      metadata <- x
    }

    # Check that the required columns are present in the metadata table
    stopifnot(all(c(cluster, sample) %in% colnames(metadata)))

    # Subset the metadata table to the required columns
    metadata <- metadata[, c(sample, cluster, additional_cols)]

    # Get the unique sample names
    samples <- unique(metadata[[sample]])

    # Calculate cell proportions for each sample
    cell_props <- purrr::map_dfr(samples, function(x) {
      # Subset the metadata table to the current sample
      cell_types <- metadata |>
        dplyr::filter(!!dplyr::sym(sample) == x) |>
        dplyr::pull({
          {
            cluster
          }
        })

      # Calculate the cell proportions for the current sample and round to 2 digits
      prop <-
        round(prop.table(table(cell_types)) * 100, digits = 2) |>
        as.data.frame()
      prop$sample <- x

      # Add additional columns to the output if specified
      if (!is.null(additional_cols)) {
        prop_additional <- metadata |>
          filter(!!sym(sample) == x) |>
          summarise(across(all_of(additional_cols), ~ first(.))) |>
          as.list() |>
          unlist()
        prop_additional <-
          prop_additional[names(prop_additional) %in% additional_cols]
        prop <- c(prop_additional, prop)
      }

      return(prop)
    })

    # Rename columns of the output data frame
    colnames(cell_props) <-
      c(additional_cols, "cluster", "proportion", "sample")

    # Return the output data frame
    return(cell_props)
  }

#' Calculate odds ratio and association statistics
#'
#' This function calculates the odds ratio and association statistics for two categorical variables in a given dataset. It also performs tissue-level normalization using Startrac.
#'
#' @param x An object of class Seurat or a data.frame containing the metadata with the variables of interest
#' @param variable A character string specifying the variable of interest
#' @param cluster A character string specifying the cluster variable, default is "cell_type"
#' @param sample A character string specifying the sample variable, default is "sample"
#' @param cuts A numeric vector specifying the cutoff values for binning the tissue-level expression levels, default is c(0, 0.8, 1.2, Inf)
#' @param levels A character vector specifying the labels for the tissue-level expression level bins, default is c("-", "+/-", "+")
#'
#' @importFrom data.table as.data.table dcast :=
#'
#' @examples
#' hp_calculate_odds_ratio(
#'   x = pbmc,
#'   variable = "CD8A",
#'   cluster = "cell_type",
#'   sample = "sample",
#'   cuts = c(0, 0.8, 1.2, Inf),
#'   levels = c("-", "+/-", "+")
#' )
#'
#' @export
hp_calculate_odds_ratio <-
  function(x,
           variable,
           cluster = NULL,
           sample = NULL,
           cuts = c(0, 0.8, 1.2, Inf),
           levels = c("-", "+/-", "+")) {
    if (inherits(x, what = "Seurat")) {
      metadata <- x@meta.data
    } else {
      metadata <- x
    }
    cluster <- cluster %||% "cell_type"
    sample <- sample %||% "sample"
    if (!all(c(variable, cluster, sample) %in% colnames(metadata))) {
      stop("Some columns are missing in metadata.")
    }

    colnames(metadata)[which(colnames(metadata) == variable)] <-
      "variable"
    colnames(metadata)[which(colnames(metadata) == cluster)] <-
      "cluster"
    metadata <- as.data.table(metadata)
    metadata[, cell_type := as.character(cell_type)]
    loc.avai.vec <- na.omit(unique(metadata[["variable"]]))
    count.dist <-
      unclass(metadata[, table(cell_type, variable)])[, loc.avai.vec]
    freq.dist <- sweep(count.dist, 1, rowSums(count.dist), "/")
    freq.dist.bin <- floor(freq.dist * 100 / 10)
    print(freq.dist.bin)

    count.dist.melt.ext.tb <-
      hp_calculate_association(count_matrix = count.dist)
    p.dist.tb <-
      dcast(count.dist.melt.ext.tb, rid ~ cid, value.var = "p_value")
    OR.dist.tb <-
      dcast(count.dist.melt.ext.tb, rid ~ cid, value.var = "odds_ratio")
    OR.dist.mtx <- as.matrix(OR.dist.tb[, -1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]

    startrac.dist <-
      unclass(
        Startrac::calTissueDist(
          metadata,
          colname.cluster = "cluster",
          colname.patient = sample,
          colname.tissue = "variable"
        )
      )
    startrac.dist <- startrac.dist[, loc.avai.vec]

    startrac.dist.bin.values <-
      factor(x = levels, levels = levels)
    startrac.dist.bin <-
      matrix(startrac.dist.bin.values[findInterval(startrac.dist, cuts)],
             ncol = ncol(startrac.dist))
    colnames(startrac.dist.bin) <- colnames(startrac.dist)
    rownames(startrac.dist.bin) <- rownames(startrac.dist)

    return(
      list(
        "count.dist.melt.ext.tb" = count.dist.melt.ext.tb,
        "p.dist.tb" = p.dist.tb,
        "OR.dist.tb" = OR.dist.tb,
        "OR.dist.mtx" = OR.dist.mtx,
        "startrac.dist.bin" = startrac.dist.bin,
        "startrac.dist" = startrac.dist
      )
    )
  }


#' Calculate association statistics for a contingency table
#'
#' This function calculates association statistics (odds ratio and p-value) for a contingency table using Fisher's exact test. The function takes a contingency table in matrix format and returns a data.table with one row for each row/column combination, along with the odds ratio and adjusted p-value.
#'
#' @param count_matrix A contingency table in matrix format, where rows represent one variable and columns represent another variable. Each cell contains the count of observations for that combination of variables.
#' @param min_row_sum The minimum row sum required for a row to be included in the analysis. Default is 0, meaning all rows are included.
#'
#' @return A data.table with the following columns:
#' \itemize{
#'  \item \code{rid}: Row name (variable 1)
#'  \item \code{cid}: Column name (variable 2)
#'  \item \code{count}: Count of observations for that row/column combination
#'  \item \code{odds_ratio}: Odds ratio for that row/column combination
#'  \item \code{p_value}: P-value for Fisher's exact test for that row/column combination
#'  \item \code{adj_p_value}: Adjusted p-value using the Benjamini-Hochberg method
#' }
#'
#' @importFrom data.table setDT melt as.data.table :=
#' @importFrom plyr ldply
#'
#' @examples
#' # Create a contingency table
#' count_matrix <- matrix(c(10, 20, 30, 40), nrow = 2, dimnames = list(c("Group A", "Group B"), c("Outcome 1", "Outcome 2")))
#'
#' # Calculate association statistics
#' hp_calculate_association(count_matrix)
#'
#' @export
hp_calculate_association <-
  function(count_matrix, min_row_sum = 0) {
    # Subset rows with row sums greater than or equal to min_row_sum
    count_matrix <-
      count_matrix[rowSums(count_matrix) >= min_row_sum, , drop = FALSE]

    # Calculate column and row sums
    col_sums <- colSums(count_matrix)
    row_sums <- rowSums(count_matrix)

    # Convert count_matrix to data.frame and set rownames as a separate column
    count_df <- as.data.frame(count_matrix)
    setDT(count_df, keep.rownames = TRUE)

    # Melt count_df to long format
    count_melt <- melt(count_df, id.vars = "rn")
    colnames(count_melt) <- c("rid", "cid", "count")

    # Calculate Fisher's exact test for each row/col combination
    count_melt_ext <-
      as.data.table(ldply(seq_len(nrow(count_melt)), function(i) {
        this_row <- count_melt$rid[i]
        this_col <- count_melt$cid[i]
        this_count <- count_melt$count[i]
        other_col_count <- col_sums[this_col] - this_count
        this_matrix <- matrix(
          c(
            this_count,
            row_sums[this_row] - this_count,
            other_col_count,
            sum(col_sums) - row_sums[this_row] - other_col_count
          ),
          ncol = 2
        )
        res_test <- fisher.test(this_matrix)
        data.frame(
          rid = this_row,
          cid = this_col,
          p_value = res_test$p.value,
          odds_ratio = res_test$estimate
        )
      }))

    # Merge count_melt_ext with original count_melt
    count_melt_ext <-
      merge(count_melt, count_melt_ext, by = c("rid", "cid"))

    # Adjust p-values using Benjamini-Hochberg method
    count_melt_ext[, adj_p_value := p.adjust(p_value, "BH")]

    # Return final data.table
    return(count_melt_ext)
  }

#' Calculate LISI score
#'
#' This function calculates the Local Intrinsic Dimensionality-based Outlier score (Lisi score) for a Seurat object.
#'
#' @param object A Seurat object.
#' @param reduction The dimensionality reduction method used to generate the embeddings (default is "pca").
#' @param label The column name in object@meta.data that specifies the sample labels (default is "sample").
#'
#' @return A data frame with the Lisi score for each cell, along with the cell ID.
#'
#' @importFrom Seurat Embeddings
#' @importFrom tibble rownames_to_column
#'
#' @examples
#' # Calculate Lisi score for a Seurat object
#' hp_calculate_lisi(x = seurat_obj)
#'
#' @export
hp_calculate_lisi <-
  function(x,
           reduction = "pca",
           label = "sample") {
    if (!requireNamespace("lisi", quietly = TRUE)) {
      stop(
        "The 'lisi' package is not installed.
           Please install it using
           devtools::install_github('immunogenomics/lisi')."
      )
    }
    lisi_score <- lisi::compute_lisi(
      X = Embeddings(x, reduction = reduction),
      meta_data = x@meta.data,
      label_colnames = label
    ) |> rownames_to_column("cell_id")
    return(lisi_score)
  }


#' Calculate ROGUE score
#'
#' This function calculates the ROGUE score for a Seurat object.
#'
#' @param x A Seurat object.
#' @param cluster The name of the column in x@meta.data that specifies the cell clusters (default is NULL).
#' @param sample The name of the column in x@meta.data that specifies the sample labels (default is NULL).
#' @param span The span parameter used in the ROGUE calculation (default is 0.9).
#'
#' @return A data frame with the ROGUE score for each cell, along with the cell ID.
#'
#' @importFrom Seurat GetAssayData
#' @importFrom Matrix as.matrix
#' @examples
#' # Calculate ROGUE score for a Seurat object
#' hp_calculate_rogue(x, cluster = "cluster", sample = "sample")
#'
#' @export
hp_calculate_rogue <-
  function(x,
           cluster = NULL,
           sample = NULL,
           span = 0.9) {
    if (!requireNamespace("ROGUE", quietly = TRUE)) {
      stop(
        "The 'ROGUE' package is not installed.
           Please install it using
           devtools::install_github('PaulingLiu/ROGUE')."
      )
    }
    counts <- GetAssayData(x, slot = "counts")
    counts <- as.matrix(counts)
    metadata <- x@meta.data
    counts <-
      ROGUE::matr.filter(counts, min.cells = 10, min.genes = 10)
    entropy <- ROGUE::SE_fun(counts)
    rogue <- ROGUE::CalculateRogue(entropy, platform = "UMI")
    if (!is.null(cluster) && !is.null(sample)) {
      rogue <-
        ROGUE::rogue(
          counts,
          labels = as.character(metadata[[cluster]]),
          samples = as.character(metadata[[sample]]),
          platform = "UMI",
          span = span
        )
    }
    return(rogue)
  }

##### From raw count matrices to high-quality cellular data --------------------
# 1. Filtering low-quality cells and noise correction
# SoupX
# Cellbender
# scDblFinder
# 2. Normalization and variance stabilization
# 3. Removing confounding sources of variation
# 4. Selecting informative features and reducing dimensionality
#### From clusters to cell identities ------------------------------------------
# 1. From single cells to clusters
# 2. Mapping cell clusters to cell identities
# 3. From discrete states to continuous processes
#### Revealing mechanisms ------------------------------------------------------
# 1. Differential gene expression analysis
# 2. Gene set enrichment analysis
# 3. Deciphering changes in cell composition
# 4. Inferring perturbation effects
# 5. Communication events across cells
