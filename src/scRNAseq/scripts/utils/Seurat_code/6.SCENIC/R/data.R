#' The PBMC dataset
#'
#' 10X Genomics' 3k PBMC dataset
#'
#' @format A Seurat object with the following slots filled
#' \describe{
#'   \item{assays}{
#'   \itemize{Currently only contains one assay
#'    ("RNA" - scRNA-seq expression data)
#'   \item{counts - Raw expression data}
#'   \item{data - Normalized expression data}
#'   \item{scale.data - Scaled expression data}
#'   \item{var.features - names of the current features selected as variable}
#'   \item{meta.features - Assay level metadata such as mean and variance}
#'    }}
#'   \item{meta.data}{Cell level metadata}
#'   \item{active.assay}{Current default assay}
#'   \item{active.ident}{Current default idents}
#'   \item{graphs}{Neighbor graphs computed, currently stores the SNN}
#'   \item{reductions}{Dimensional reductions: currently PCA and tSNE}
#'   \item{version}{Seurat version used to create the object}
#'   \item{commands}{Command history}
#' }
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k}
#'
"pbmc"


#' Get metadata for RNA-seq or MeRIP-seq samples
#'
#' This function returns a data frame containing metadata for RNA-seq or
#' MeRIP-seq samples, including sample names, grouping information, and file
#' paths for input and IP fastq files.
#'
#' @param type A character string indicating the type of sequencing data.
#'   Default is "rna-seq".
#' @return A data frame containing metadata for RNA-seq or MeRIP-seq samples.
#'
#' @examples
#' metadata_rna <- get_metadata(type = "rna-seq")
#' metadata_merip <- get_metadata(type = "merip-seq")
#'
#' @export
get_metadata <- function(type = "rna-seq") {
  fastq_files <- system.file("extdata/fastq", package = "hephaestus")
  list.files(fastq_files)
  # sample
  sample <- c("WT1", "WT2", "KO1", "KO2")
  group <- c("WT", "WT", "KO", "KO")
  # input
  input_r1 <- c(
    file.path(fastq_files, "WT1.input.R1.fastq.gz"),
    file.path(fastq_files, "WT2.input.R1.fastq.gz"),
    file.path(fastq_files, "KO1.input.R1.fastq.gz"),
    file.path(fastq_files, "KO2.input.R1.fastq.gz")
  )
  input_r2 <- c(
    file.path(fastq_files, "WT1.input.R2.fastq.gz"),
    file.path(fastq_files, "WT2.input.R2.fastq.gz"),
    file.path(fastq_files, "KO1.input.R2.fastq.gz"),
    file.path(fastq_files, "KO2.input.R2.fastq.gz")
  )
  if (type == "rna-seq") {
    metadata <- data.frame(
      sample = sample,
      group = group,
      input_r1 = input_r1,
      input_r2 = input_r2
    )
    return(metadata)
  } else if (type == "merip-seq") {
    ip_r1 <- c(
      file.path(fastq_files, "WT1.ip.R1.fastq.gz"),
      file.path(fastq_files, "WT2.ip.R1.fastq.gz"),
      file.path(fastq_files, "KO1.ip.R1.fastq.gz"),
      file.path(fastq_files, "KO2.ip.R1.fastq.gz")
    )
    ip_r2 <- c(
      file.path(fastq_files, "WT1.ip.R2.fastq.gz"),
      file.path(fastq_files, "WT2.ip.R2.fastq.gz"),
      file.path(fastq_files, "KO1.ip.R2.fastq.gz"),
      file.path(fastq_files, "KO2.ip.R2.fastq.gz")
    )
    metadata <- data.frame(
      sample = sample,
      group = group,
      input_r1 = input_r1,
      input_r2 = input_r2,
      ip_r1 = ip_r1,
      ip_r2 = ip_r2
    )
    return(metadata)
  }
}
