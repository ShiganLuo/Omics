#' @include utilities.R
NULL
#' FastQC - A high throughput sequence QC analysis tool
#'
#' @param x fastq files
#' @param fastqc The absolute path of the FastQC
#' @param outdir Directory to store the output files
#' @param ncores Max of threads to process the data
#'
#' @return
#' @importFrom tictoc tic toc
#' @importFrom logger log_info
#' @export
#'
#' @examples
#' library(hephaestus)
#'
#' x <- create_hephaestus()
#' cat_fastqc(
#'   x = x,
#'   fastqc = "/opt/fastqc/0.11.9/bin/fastqc"
#' )
cat_fastqc <- function(x,
                       fastqc = NULL,
                       outdir = "./fastqc",
                       ncores = 1) {
  if (is.null(fastqc)) {
    fastqc <- "/opt/fastqc/0.11.9/bin/fastqc"
    hp_check_file(fastqc)
  }
  # x可以是Hephaestus对象，也可以是一个fastq文件路径
  if (class(x) == "Hephaestus") {
    board <- x@board
    fastq <- c(
      "input_r1",
      "input_r2",
      "ip_r1",
      "ip_r2",
      "input_clean_r1",
      "input_clean_r2",
      "ip_clean_r1",
      "ip_clean_r2"
    )
    fastq <- fastq[fastq %in% colnames(board)]
    input <- lapply(fastq, function(x) {
      board[[x]]
    }) |>
      unlist()
    outdir <- file.path(x@outdir, "fastqc")
  } else {
    input <- x
  }
  outdir <- hp_set_path(outdir)
  # prefix <- cat_filename(input)
  prefix <- basename(gsub(".fastq.gz|.fq.gz", "", input))
  output <- file.path(outdir, paste0(prefix, "_fastqc.html"))
  adata <- data.frame(prefix = prefix,
                      input = input,
                      output = output)
  hp_check_file(adata$input)
  output_not_exists <- hp_check_file(adata$output, stop = FALSE)
  input <- adata[adata$output %in% output_not_exists,]$input

  n <- length(input)
  if (n > 0) {
    if (ncores == 1 && n <= 20) {
      ncores <- n
    } else {
      ncores <- 20
    }
    log_info("Data quality assessment on {n} files using {ncores} threads")
    # tic(msg = "Data quality assessment")
    system2(
      command = fastqc,
      args = c("--outdir", outdir,
               "--threads", ncores,
               input)
    )
    # toc()
  }
  if (class(x) == "Hephaestus") {
    for (i in fastq) {
      x@board[[paste0(i, "_fastqc")]] <- gsub(
        pattern = ".fastq.gz",
        replacement = "_fastqc.html",
        x = file.path(outdir, basename(board[[i]]))
      )
      x@board[[paste0(i, "_fastqc")]] <- gsub(
        pattern = ".fq.gz",
        replacement = "_fastqc.html",
        x = file.path(outdir, basename(board[[i]]))
      )
    }
    return(x)
  } else {
    return(output)
  }
}
#' A wrapper tool around Cutadapt and FastQC to consistently apply quality and
#' adapter trimming to FastQ files, with some extra functionality for
#' MspI-digested RRBS-type (Reduced Representation Bisufite-Seq) libraries.
#'
#' @param x fastq files
#' @param trim_galore The absolute path of the Trim Galore!
#' @param cutadapt The absolute path of the cutadapt
#' @param outdir Directory to store the output files
#' @param ncores Max of threads to process the data
#'
#' @return clean data
#' @importFrom logger log_info
#' @importFrom tictoc tic toc
#' @export
#'
#' @examples
#' library(hephaestus)
#'
#' x <- create_hephaestus()
#' x <- cat_trim_galore(
#'   x,
#'   trim_galore = "/opt/trim-galore/0.6.10/bin/trim_galore",
#'   cutadapt = "/opt/cutadapt/4.2/bin/cutadapt"
#' )
cat_trim_galore <- function(x = NULL,
                            trim_galore = NULL,
                            cutadapt = NULL,
                            outdir = "./clean",
                            ncores = 8) {
  if (is.null(trim_galore)) {
    trim_galore <- "/opt/trim-galore/0.6.10/bin/trim_galore"
    hp_check_file(trim_galore)
  }
  if (is.null(cutadapt)) {
    cutadapt <- "/opt/cutadapt/4.2/bin/cutadapt"
    hp_check_file(cutadapt)
  }

  # outdir
  outdir <- hp_set_path(file.path(x@outdir, "clean"))
  # board
  board <- x@board
  # sample name
  sample_name <- x$sample
  # input
  input_r1 <- board$input_r1
  input_r2 <- board$input_r2
  input_basename <- paste(sample_name, "input", sep = "_")
  input_clean_r1 <-
    file.path(outdir, paste(input_basename, "val_1.fq.gz", sep = "_"))
  input_clean_r2 <-
    file.path(outdir, paste(input_basename, "val_2.fq.gz", sep = "_"))
  # ip
  if (all(c("ip_r1", "ip_r2") %in% colnames(board))) {
    ip_r1 <- board$ip_r1
    ip_r2 <- board$ip_r2
    ip_basename <- paste(sample_name, "ip", sep = "_")
    ip_clean_r1 <-
      file.path(outdir, paste(ip_basename, "val_1.fq.gz", sep = "_"))
    ip_clean_r2 <-
      file.path(outdir, paste(ip_basename, "val_2.fq.gz", sep = "_"))
  } else {
    ip_r1 <- NULL
    ip_r2 <- NULL
    ip_basename <- NULL
    ip_clean_r1 <- NULL
    ip_clean_r2 <- NULL
  }

  basename <- c(input_basename, ip_basename)
  r1 <- c(input_r1, ip_r1)
  r2 <- c(input_r2, ip_r2)
  clean_r1 <- c(input_clean_r1, ip_clean_r1)
  clean_r2 <- c(input_clean_r2, ip_clean_r2)
  # PREFERRED_NAME_val_1.fq(.gz) and PREFERRED_NAME_val_2.fq(.gz)
  # sample_name_input_val_1.fq.gz
  # sample_name_ip_val_1.fq.gz
  for (i in seq_along(basename)) {
    temp_r1 <- r1[i]
    temp_r2 <- r2[i]
    if (!all(file.exists(clean_r1[i], clean_r2[i]))) {
      log_info("Start quality and adapter trimming...")
      system2(
        command = trim_galore,
        args = c(
          # "--fastqc",
          # "--fastqc_args",
          # fastqc_args,
          "--basename",
          basename[i],
          "--cores",
          ncores,
          "--paired",
          temp_r1,
          temp_r2,
          "--path_to_cutadapt",
          cutadapt,
          "--gzip",
          "--output_dir",
          outdir
        )
      )
    } else {
      log_info("The file already exists, please check if it is complete?")
    }
  }

  x@board$input_clean_r1 <- input_clean_r1
  x@board$input_clean_r2 <- input_clean_r2
  if (all(c("ip_r1", "ip_r2") %in% colnames(x@board))) {
    x@board$ip_clean_r1 <- ip_clean_r1
    x@board$ip_clean_r2 <- ip_clean_r2
  }
  return(x)
}

#' Reading/writing/editing/indexing/viewing SAM/BAM/CRAM format
#'
#' @param samtools
#' @param command
#' @param quality
#' @param ncores
#' @param bam
#' @param sam
#'
#' @return
#' @export
#'
#' @examples
cat_samtools <-
  function(samtools = NULL,
           command,
           quality = NULL,
           ncores = 10,
           bam,
           sam) {
    if (is.null(samtools)) {
      samtools <- "/opt/samtools/1.16.1/bin/samtools"
      hp_check_file(samtools)
    }
    match.arg(arg = command, choices = c("view", "sort", "index"))
    # veiw
    if (command == "view") {
      system2(
        command = samtools,
        args = c(
          command,
          "--bam",
          "--with-header",
          "-S",
          if (!is.null(quality)) {
            c("--min-MQ",
              quality)
          },
          "--threads",
          ncores,
          sam,
          "--output",
          bam
        )
      )
    }
    # sort
    if (command == "sort") {
      system2(
        command = samtools,
        args = c(command,
                 "--threads",
                 ncores,
                 "-o",
                 bam,
                 sam)
      )
    }
    # index
    if (command == "index") {
      system2(command = samtools,
              args = c(command,
                       bam))
    }
  }

#' Graph-based alignment of next generation sequencing reads to
#' a population of genomes
#'
#' @param x
#' @param hisat2
#' @param species
#' @param index
#' @param ncores
#' @param workers
#' @param outdir
#'
#' @return
#' @importFrom furrr future_walk
#' @importFrom future plan multisession
#' @importFrom logger log_info
#' @export
#'
#' @examples
cat_hisat2 <- function(x,
                       hisat2 = NULL,
                       species = "human",
                       index = NULL,
                       ncores = 5,
                       workers = "auto",
                       outdir = "./align") {
  # software
  if (is.null(hisat2)) {
    hisat2 <- "/opt/hisat2/0.2.2/bin/hisat2"
    hp_check_file(hisat2)
  }

  # index
  if (is.null(index)) {
    if (species == "human") {
      index <- "/DATA/public/index/hisat2/homo_sapiens/GRCh38.p13/genome"
    } else if (species == "mouse") {
      index <- "/DATA/public/index/hisat2/mus_musculus/GRCm38.p6/genome"
    }
  }
  hp_check_file(paste0(index, ".1.ht2"))

  if (class(x) == "Hephaestus") {
    # outdir
    outdir <- hp_set_path(file.path(x@outdir, "align"))
    # sample name
    sample <- x$sample
    # input
    input_clean_r1 <- x@board$input_clean_r1
    input_clean_r2 <- x@board$input_clean_r2
    input_sample <- paste(sample, "input", sep = "_")
    input_sam <- file.path(outdir, paste0(input_sample, ".sam"))
    # ip
    if (all(c("ip_clean_r1", "ip_clean_r1") %in% colnames(x@board))) {
      ip_clean_r1 <- x@board$ip_clean_r1
      ip_clean_r2 <- x@board$ip_clean_r2
      ip_sample <- paste(sample, "ip", sep = "_")
      ip_sam <- file.path(outdir, paste0(ip_sample, ".sam"))
    } else {
      ip_clean_r1 <- NULL
      ip_clean_r2 <- NULL
      ip_sample <- NULL
      ip_sam <- NULL
    }
    sample <- c(input_sample, ip_sample)
    clean_r1 <- c(input_clean_r1, ip_clean_r1)
    clean_r2 <- c(input_clean_r2, ip_clean_r2)
    sam <- c(input_sam, ip_sam)
  } else {
    stop("I'm not hephaestus!")
  }
  # print(clean_r1)
  # print(clean_r2)
  # Set a "plan" for how the code should run.
  n <- length(sam[!file.exists(sam)])
  if (n > 0) {
    if (workers == "auto" && n <= 4) {
      workers <- n
    } else if (workers == "auto" && n > 4) {
      workers <- 4
    }
    log_info("{n} samples start aligning to genomes...")
    log_info("Use {workers} workers, each assigned {ncores} cores...")
    plan(multisession, workers = workers)
  }
  # This does run in parallel!
  future_walk(
    .x = seq_along(sample),
    .f = function(i) {
      if (!file.exists(sam[i])) {
        # log_info("{sample[i]} start aligning to genomes...")
        system2(
          command = hisat2,
          args = c(
            "--threads",
            ncores,
            "--dta",
            "-x",
            index,
            "--summary-file",
            file.path(outdir, paste0(sample[i], ".summary.txt")),
            "-1",
            clean_r1[i],
            "-2",
            clean_r2[i],
            "-S",
            sam[i]
          ),
          stdout = NULL,
          stderr = NULL
        )
      }
    }
  )
  # log_tictoc("{n} samples aligned completed!")
  # mapping rate
  input_mapping_rate <- sapply(input_sample, function(i) {
    mapping_rate <-
      read.table(file = file.path(outdir, paste0(i, ".summary.txt")),
                 skip = 14)[, 1]
    mapping_rate <- as.numeric(gsub("%", "", mapping_rate))
    return(mapping_rate)
  })
  x@board$input_sam <- input_sam
  x$input_mapping_rate <- input_mapping_rate
  if (!is.null(ip_sample)) {
    ip_mapping_rate <- sapply(ip_sample, function(i) {
      mapping_rate <-
        read.table(file = file.path(outdir, paste0(i, ".summary.txt")),
                   skip = 14)[, 1]
      mapping_rate <- as.numeric(gsub("%", "", mapping_rate))
      return(mapping_rate)
    })
    x@board$ip_sam <- ip_sam
    x$ip_mapping_rate <- ip_mapping_rate
  }

  return(x)
}

#' Sambamba: process your BAM data faster!
#'
#' @param sambamba
#' @param command
#' @param input
#' @param sam_input
#' @param output
#' @param ncores
#'
#' @return
#' @importFrom furrr future_walk
#' @importFrom future plan multisession
#' @importFrom logger log_info
#' @export
#'
#' @examples
cat_sambamba <- function(x,
                         sambamba = NULL,
                         outdir = "align",
                         quality = 20,
                         format = "bam",
                         remove_duplicates = FALSE,
                         workers = "auto") {
  # software
  if (is.null(sambamba)) {
    sambamba <- "/opt/sambamba/1.0/bin/sambamba"
    hp_check_file(sambamba)
  }

  # config
  if (class(x) == "Hephaestus") {
    # outdir
    outdir <- hp_set_path(file.path(x@outdir, "align"))
    # sample name
    sample <- x$sample
    # quality
    if (!is.null(quality)) {
      quality_label <- paste0(".quality", quality)
    } else {
      quality_label <- NULL
    }
    # input
    input_sample <- paste(sample, "input", sep = "_")
    input_sam <-
      file.path(outdir, paste0(input_sample, ".sam"))
    input_bam <-
      file.path(outdir, paste0(input_sample, quality_label, ".bam"))
    input_sorted_bam <-
      file.path(outdir,
                paste0(input_sample,
                       quality_label,
                       if (remove_duplicates) {
                         ".rmdup"
                       },
                       ".sorted.bam"))
    # ip
    if ("ip_sam" %in% colnames(x@board)) {
      ip_sample <- paste(sample, "ip", sep = "_")
      ip_sam <- file.path(outdir, paste0(ip_sample, ".sam"))
      ip_bam <-
        file.path(outdir, paste0(ip_sample, quality_label, ".bam"))
      ip_sorted_bam <-
        file.path(outdir,
                  paste0(ip_sample,
                         quality_label,
                         if (remove_duplicates) {
                           ".rmdup"
                         },
                         ".sorted.bam"))
    } else {
      ip_sample <- NULL
      ip_sam <- NULL
      ip_sorted_bam <- NULL
    }
    sample <- c(input_sample, ip_sample)
    sam <- c(input_sam, ip_sam)
    bam <- c(input_bam, ip_bam)
    sorted_bam <- c(input_sorted_bam,
                    ip_sorted_bam)
  } else {
    stop("I'm not hephaestus!")
  }

  # Set a "plan" for how the code should run.
  n <- length(sorted_bam[!file.exists(sorted_bam)])
  if (n > 0) {
    if (workers == "auto" && n <= 4) {
      workers <- n
    } else if (workers == "auto" && n > 4) {
      workers <- 4
    }
    log_info("Process SAM/BAM data...")
    log_info("Use {workers} workers...")
    plan(multisession, workers = workers)
  }
  # run
  furrr::future_walk(
    .x = seq_along(sample),
    .f = function(i) {
      # sam to bam
      if (!file.exists(bam[i])) {
        # log_info("Convert {sample[i]} SAM file to BAM file...")
        system2(
          command = sambamba,
          args = c(
            "view",
            # "--nthreads",
            # ncores,
            "--format",
            format,
            if (!is.null(quality)) {
              c("--filter",
                paste0("'", "mapping_quality >= ", quality, "'"))
            },
            "--sam-input",
            sam[i],
            "--output-filename",
            bam[i]
          ),
          stdout = NULL,
          stderr = NULL
        )
      }
      # remove duplicate
      if (remove_duplicates) {
        x <- bam[i]
        rmdup_bam <- gsub(pattern = ".bam",
                          replacement = ".rmdup.bam",
                          x = x)
        bam[i] <- rmdup_bam
        if (!file.exists(rmdup_bam)) {
          # log_info("Remove {sample[i]} duplicates...")
          system2(
            command = sambamba,
            args = c("markdup",
                     "--remove-duplicates",
                     # "--nthreads",
                     # ncores,
                     x,
                     rmdup_bam),
            stdout = NULL,
            stderr = NULL
          )
        }
      }
      # bam sorted
      if (!file.exists(sorted_bam[i])) {
        # log_info("Start to sort {sample[i]} BAM genomic coordinates...")
        system2(
          command = sambamba,
          args = c("sort",
                   # "--nthreads",
                   # ncores,
                   bam[i]),
          stdout = NULL,
          stderr = NULL
        )
      }
    }
  )
  # log_tictoc("{n} samples process completed!")

  x@board$input_bam <- input_sorted_bam
  if ("ip_sam" %in% colnames(x@board)) {
    x@board$ip_bam <- ip_sorted_bam
  }
  return(x)
}

#' Graph-based alignment of next generation sequencing reads to
#' a population of genomes
#'
#' @param x
#' @param hisat2
#' @param samtools
#' @param species
#' @param index
#' @param ncores
#' @param quality
#' @param outdir
#'
#' @return
#' @importFrom stringr str_split_fixed
#' @importFrom logger log_info
#' @importFrom tictoc tic toc
#' @export
#'
#' @examples
cat_align <- function(x = NULL,
                      hisat2 = NULL,
                      sambamba = NULL,
                      species = "human",
                      index = NULL,
                      quality = 20,
                      format = "bam",
                      remove_duplicates = FALSE,
                      ncores = 8,
                      workers = "auto",
                      outdir = "./align") {
  # hisat2
  x <-
    cat_hisat2(
      x = x,
      hisat2 = hisat2,
      species = species,
      index = index,
      outdir = outdir,
      ncores = ncores,
      workers = workers
    )
  # sambamba
  x <- cat_sambamba(
    x = x,
    sambamba = sambamba,
    outdir = outdir,
    quality = quality,
    format = format,
    remove_duplicates = remove_duplicates,
    workers = workers
  )

  return(x)
}

#' A ultrafast and accurate read summarization program
#'
#' @param x
#' @param featurecounts
#' @param species
#' @param gtf
#' @param ncores
#' @param outdir
#'
#' @return
#' @importFrom SeuratObject CreateAssayObject
#' @importFrom dplyr ends_with select
#' @importFrom tibble column_to_rownames
#' @export
#'
#' @examples
cat_featurecounts <- function(x,
                              featurecounts = NULL,
                              species = "human",
                              gtf = NULL,
                              ncores = 10,
                              outdir = "./counts") {
  if (is.null(featurecounts)) {
    featurecounts <- "/opt/subread/2.0.3/bin/featureCounts"
    hp_check_file(featurecounts)
  }

  if (is.null(gtf)) {
    if (species == "human") {
      gtf <-
        "/DATA/public/genome/homo_sapiens/GRCh38.p13/gencode.v43.annotation.gtf"
    } else if (species == "mouse") {
      gtf <-
        "/DATA/public/genome/mus_musculus/GRCm38.p6/gencode.vM25.annotation.gtf"
    }
  }
  hp_check_file(gtf)

  if (class(x) == "Hephaestus") {
    bam <- x@board$input_bam
    outdir <- file.path(x@outdir, "counts")
  } else {
    bam <- x
  }
  outdir <- hp_set_path(outdir)
  read_counts <- file.path(outdir, "counts.txt")

  if (!file.exists(read_counts)) {
    system2(
      command = featurecounts,
      args = c(
        "-a",
        gtf,
        "-o",
        read_counts,
        bam,
        "-t",
        "exon",
        "-g",
        "gene_name",
        "-T",
        ncores,
        "-p"
      )
    )
  }

  if (class(x) == "Hephaestus") {
    adata <- read.table(read_counts, header = TRUE)
    counts <- adata |>
      select(Geneid, ends_with("bam")) |>
      column_to_rownames("Geneid")
    metadata <- adata |>
      select(-ends_with("bam")) |>
      column_to_rownames("Geneid")

    colnames(counts) <- colnames(x)
    counts <- counts[rowSums(counts) > 0,]
    x[["RNA"]] <-
      CreateAssayObject(counts = counts)
    x@feature[["RNA"]] <- metadata[rownames(counts),]
    return(x)
  } else {
    return(read_counts)
  }
}




#' Detect differential alternative splicing events from RNA-Seq data
#'
#' @param x
#' @param control
#' @param treatment
#' @param python
#' @param rmats
#' @param species
#' @param gtf
#' @param outdir
#' @param t
#' @param read_length
#' @param ncores
#'
#' @return
#' @importFrom purrr map_dfr
#' @importFrom dplyr mutate
#' @export
#'
#' @examples
cat_rmats <- function(x = NULL,
                      control = NULL,
                      treatment = NULL,
                      python = NULL,
                      rmats = NULL,
                      species = "human",
                      gtf = NULL,
                      outdir = "./rmats",
                      t = "paired",
                      read_length = 150,
                      psi_threshold = 0.1,
                      p_value_threshold = NULL,
                      fdr_threshold = 0.05,
                      ncores = 10) {
  if (is.null(rmats)) {
    rmats <- "/opt/rmats/4.1.2/bin/rmats.py"
    hp_check_file(rmats)
  }
  # Get python path
  if (is.null(python)) {
    python <- file.path(dirname(rmats), "python")
    hp_check_file(python)
  }
  if (!is.null(x) && class(x) == "Hephaestus") {
    bam <- x$bam
    case <- paste0(treatment, "_vs_", control)
    outdir <-
      file.path(x$outdir, "rmats", case)
    if (!is.null(control)) {
      control_index <- which(x$group %in% control)
      control <- x$bam[control_index]
    }
    if (!is.null(treatment)) {
      treatment_index <- which(x$group %in% treatment)
      treatment <- x$bam[treatment_index]
    }
    outdir <- file.path(x@board, "counts")
  }
  outdir <- hp_set_path(outdir)
  if (is.null(gtf)) {
    if (species == "human") {
      gtf <-
        "/DATA/public/genome/homo_sapiens/GRCh38.p13/gencode.v42.annotation.gtf"
    } else if (species == "mouse") {
      gtf <-
        "/DATA/public/genome/mus_musculus/GRCm39/gencode.vM31.annotation.gtf"
    }
  }
  # Preprocessing
  # b1
  if (!is.null(control)) {
    if (all(file.exists(control))) {
      b1 <- file.path(outdir, "b1.txt")
      write.table(
        paste0(control, collapse = ","),
        file = b1,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
      )
    } else {
      stop(
        paste0(
          "Please check the file path, the file listed below does not exist:\n"
        ),
        paste0(control[!file.exists(control)], collapse = "\n")
      )
    }
  }
  # b2
  if (!is.null(treatment)) {
    if (all(file.exists(treatment))) {
      b2 <- file.path(outdir, "b2.txt")
      write.table(
        paste0(treatment, collapse = ","),
        file = b2,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
      )
    } else {
      stop(
        paste0(
          "Please check the file path, the file listed below does not exist:\n"
        ),
        paste0(treatment[!file.exists(treatment)], collapse = "\n")
      )
    }
  }
  summary_file <- file.path(outdir, "rmats_summary.tsv")
  if (!file.exists(summary_file)) {
    # rmats
    system2(
      command = python,
      args = c(
        rmats,
        "--gtf",
        gtf,
        "--b1",
        b1,
        if (!is.null(treatment)) {
          c("--b2", b2)
        },
        "--od",
        outdir,
        "--tmp",
        file.path(outdir, "tmp_output"),
        "-t",
        t,
        "--readLength",
        read_length,
        "--variable-read-length",
        "--nthread",
        ncores,
        if (is.null(treatment)) {
          "--statoff"
        }
      )
    )

    # Organize results
    as_events <- c("SE", "A5SS", "A3SS", "MXE", "RI")
    adata <- map_dfr(as_events, function(x) {
      read.delim(file.path(outdir,
                           paste0(x,
                                  ".MATS.JC.txt"))) |>
        mutate(as_event = x)
    })
    if (!is.null(p_value_threshold)) {
      adata <- adata |>
        mutate(change = ifelse(
          PValue <= p_value_threshold &
            abs(IncLevelDifference) >= psi_threshold,
          ifelse(
            IncLevelDifference >= psi_threshold,
            "Upregulated",
            "Downregulated"
          ),
          "Stable"
        ))
    } else {
      adata <- adata |>
        mutate(change = ifelse(
          FDR <= fdr_threshold &
            abs(IncLevelDifference) >= psi_threshold,
          ifelse(
            IncLevelDifference >= psi_threshold,
            "Upregulated",
            "Downregulated"
          ),
          "Stable"
        ))
    }
    write.table(adata,
                summary_file,
                row.names = FALSE,
                quote = FALSE)
  } else {
    adata <- read.table(summary_file, header = TRUE)
  }

  if (!is.null(x) && class(x) == "Hephaestus") {
    x$rmart[[case]]$summary <- adata
    return(x)
  } else {
    return(adata)
  }
}

#' Dynamics analysis of alternative polyadenylation from multiple RNA-seq data
#'
#' @param x
#' @param outdir
#' @param python
#' @param bedtools
#' @param samtools
#' @param species
#' @param annotated_3utr
#' @param chr_list
#' @param ncores
#'
#' @return
#' @importFrom dplyr mutate
#' @importFrom tictoc tic toc
#' @importFrom stringr str_split_fixed
#' @export
#'
#' @examples
.dapars2 <-
  function(x,
           python = NULL,
           bedtools = NULL,
           samtools = NULL,
           species = "human",
           annotated_3utr = NULL,
           chr_list = NULL,
           sample_name = NULL,
           outdir = "./dapars2",
           ncores = 10) {
    if (is.null(python)) {
      python <- "/opt/anaconda3/bin/python"
      hp_check_file(python)
    }
    if (is.null(bedtools)) {
      bedtools <- "/opt/bedtools/2.30.0/bin/bedtools"
      hp_check_file(bedtools)
    }
    if (is.null(samtools)) {
      samtools <- "/opt/samtools/1.16.1/bin/samtools"
      hp_check_file(samtools)
    }
    if (!is.null(x) && class(x) == "Hephaestus") {
      bam <- x$align_data$bam
      outdir <- file.path(x$outdir, "dapars2")
    } else {
      bam <- x
    }
    wig_path <- hp_set_path(file.path(outdir, "wig"))
    script_dir <-
      system.file("script/dapars2", package = "hephaestus")
    extdata <-
      system.file("extdata/dapars2", package = "hephaestus")
    if (species == "human") {
      annotated_3utr <-
        file.path(extdata, "hg38.vM43.3UTR.annotation.bed")
      chr_list <- file.path(extdata, "hg38.chr.list.txt")
    } else if (species == "mouse") {
      annotated_3utr <-
        file.path(extdata, "mm10.vM25.3UTR.annotation.bed")
      chr_list <- file.path(extdata, "mm10.chr.list.txt")
    }
    mapping_bam_location_with_depth_path <-
      file.path(outdir, "mapping_bam_location_with_depth.txt")
    summary_file <- file.path(outdir, "dapars2_summary.csv")
    if (!file.exists(summary_file)) {
      # pre-processing
      # run
      log_info("Convert Bam files into bedgraph format using bedtools...")
      depth_list <- lapply(bam, function(i) {
        if (is.null(sample_name)) {
          s_sample_name <- cat_filename(i)
        } else {
          s_sample_name <- sample_name[which(i == bam)]
        }
        wig_filename <- paste0(s_sample_name, ".wig")
        if (!file.exists(file.path(wig_path, wig_filename))) {
          tic(paste0("bedtools genomecov ", s_sample_name))
          system2(
            command = bedtools,
            args = c(
              "genomecov",
              "-ibam",
              i,
              "-bga",
              "-split",
              "-trackline",
              ">",
              file.path(wig_path, wig_filename)
            )
          )
          toc()
        }
        depth <- system2(
          command = samtools,
          args = c("flagstat",
                   i, "|head -n 1|cut -d+ -f1"),
          stdout = TRUE
        )
        depth <- as.numeric(depth)
        mapping_bam_location_with_depth <-
          data.frame(x = file.path(wig_path, wig_filename), y = depth)
        return(mapping_bam_location_with_depth)
      })
      mapping_bam_location_with_depth <-
        Reduce(f = rbind, x = depth_list)
      write.table(
        mapping_bam_location_with_depth,
        file = mapping_bam_location_with_depth_path,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )

      dapars2_configure_file_path <-
        file.path(outdir, "dapars2_configure_file")
      cat(
        c(
          "# Specify the reference of 3'UTR region",
          paste0("\nAnnotated_3UTR=", annotated_3utr),
          "\n# A comma separated list of wig files of all samples",
          paste0(
            "\nAligned_Wig_files=",
            paste(mapping_bam_location_with_depth[, 1], collapse = ",")
          ),
          paste0("\nOutput_directory=", file.path(outdir, "results/")),
          "\nOutput_result_file=dapars2",
          "\n# Specify Coverage threshold",
          "Coverage_threshold=10",
          "\n# Specify the number of threads to process the analysis",
          paste0("Num_Threads=", ncores),
          "\n# Provide sequencing depth file for normalization",
          paste0(
            "\nsequencing_depth_file=",
            mapping_bam_location_with_depth_path
          )
        ),
        file = dapars2_configure_file_path
      )
      system2(
        command = python,
        args = c(
          file.path(script_dir, "DaPars2_Multi_Sample_Multi_Chr.py"),
          dapars2_configure_file_path,
          chr_list
        )
      )
      # summary results
      chr <- read.table(chr_list)[, 1]
      adata <- lapply(chr, function(x) {
        read.delim(file.path(
          outdir,
          paste0("results_", x),
          paste0("dapars2_result_temp.", x, ".txt")
        ))
      })
      adata <- Reduce(f = rbind, x = adata)
      sample_name_index <- 5:ncol(adata)
      if (is.null(sample_name)) {
        colnames(adata)[sample_name_index] <- x$metadata$sample
      } else {
        colnames(adata)[sample_name_index] <- sample_name
      }

      gene_info <- str_split_fixed(adata$Gene, "\\|", 4)
      adata <- adata |>
        mutate(
          gene_id = gene_info[, 1],
          gene_symbol = gene_info[, 2],
          chr = gene_info[, 3],
          strand = gene_info[, 4]
        )
      write.csv(adata,
                file = summary_file,
                quote = FALSE,
                row.names = FALSE)
    } else {
      adata <- read.csv(file = summary_file)
    }
    if (!is.null(x) && class(x) == "Hephaestus") {
      x$dapars2$summary <- adata
      return(x)
    } else {
      return(adata)
    }
  }


#' Dynamics analysis of alternative polyadenylation from multiple RNA-seq data
#'
#' @param x
#' @param outdir
#' @param python
#' @param bedtools
#' @param samtools
#' @param species
#' @param annotated_3utr
#' @param chr_list
#' @param ncores
#'
#' @return
#' @importFrom dplyr mutate
#' @importFrom tictoc tic toc
#' @importFrom stringr str_split_fixed
#' @export
#'
#' @examples
.dapars <-
  function(control_bams,
           treatment_bams,
           control_sample_names = NULL,
           treatment_sample_names = NULL,
           python = NULL,
           bedtools = NULL,
           samtools = NULL,
           species = "mouse",
           annotated_3utr = NULL,
           chr_list = NULL,
           fdr_cutoff = 0.05,
           pdui_cutoff = 0.2,
           fold_change_cutoff = 0.59,
           outdir = "./dapars") {
    if (is.null(python)) {
      python <- "/opt/anaconda3/bin/python"
      hp_check_file(python)
    }
    if (is.null(bedtools)) {
      bedtools <- "/opt/bedtools/2.30.0/bin/bedtools"
      hp_check_file(bedtools)
    }
    if (is.null(samtools)) {
      samtools <- "/opt/samtools/1.16.1/bin/samtools"
      hp_check_file(samtools)
    }
    # if (!is.null(x) && class(x) == "Hephaestus") {
    #   bam <- x$align_data$bam
    #   outdir <- file.path(x$outdir, "dapars")
    # } else {
    #   bam <- x
    # }
    wig_path <- hp_set_path(file.path(outdir, "wig"))
    # script_dir <-
    #   system.file("script/dapars", package = "hephaestus")
    # extdata <-
    #   system.file("extdata/dapars2", package = "hephaestus")
    if (species == "human") {
      annotated_3utr <-
        file.path("/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/inst/extdata/dapars2/hg38.vM43.3UTR.annotation.bed")
    } else if (species == "mouse") {
      annotated_3utr <-
        file.path("/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/inst/extdata/dapars2/mm10.vM25.3UTR.annotation.bed")
    }
      # pre-processing
      # run
      log_info("Convert Bam files into bedgraph format using bedtools...")
      control_wigs <- sapply(control_bams, function(i) {
        s_sample_name <- control_sample_names[which(i == control_bams)]
        wig_filename <-
          file.path(wig_path, paste0(s_sample_name, ".wig"))

        if (!file.exists(wig_filename)) {
          tic(paste0("bedtools genomecov ", s_sample_name))
          system2(
            command = bedtools,
            args = c(
              "genomecov",
              "-ibam",
              i,
              "-bga",
              "-split",
              "-trackline",
              ">",
              wig_filename
            )
          )
          toc()
        }
        return(wig_filename)
      })
      # treatment
      treatment_wigs <- sapply(treatment_bams, function(i) {
        s_sample_name <- treatment_sample_names[which(i == treatment_bams)]
        # wig_filename <- paste0(s_sample_name, ".wig")
        wig_filename <-
          file.path(wig_path, paste0(s_sample_name, ".wig"))
        if (!file.exists(wig_filename)) {
          tic(paste0("bedtools genomecov ", s_sample_name))
          system2(
            command = bedtools,
            args = c(
              "genomecov",
              "-ibam",
              i,
              "-bga",
              "-split",
              "-trackline",
              ">",
              wig_filename
            )
          )
          toc()
        }
        return(wig_filename)
      })

      dapars2_configure_file_path <-
        file.path(outdir, "dapars_configure_file")
      cat(
        c(
          paste0("Annotated_3UTR=", annotated_3utr),
          paste0(
            "\nGroup1_Tophat_aligned_Wig=",
            paste(treatment_wigs, collapse = ",")
          ),
          paste0(
            "\nGroup2_Tophat_aligned_Wig=",
            paste(control_wigs, collapse = ",")
          ),
          paste0("\nOutput_directory=", file.path(outdir, "results/")),
          "\nOutput_result_file=dapars",
          "\n#Parameters",
          "\nNum_least_in_group1=1",
          "\nNum_least_in_group2=1",
          "\nCoverage_cutoff=30",
          paste0("\nFDR_cutoff=", fdr_cutoff),
          paste0("\nPDUI_cutoff=", pdui_cutoff),
          paste0("\nFold_change_cutoff=", fold_change_cutoff)
        ),
        file = dapars2_configure_file_path
      )
      dapars_results_path <- file.path(outdir,
                                       "results",
                                       "dapars_All_Prediction_Results.txt")
      if (!file.exists(dapars_results_path)) {
        processx::run(
          command = python,
          args = c(
            "/home/lixy/R/x86_64-pc-linux-gnu-library/4.3/hephaestus/inst/script/dapars/DaPars_main.py",
            "dapars_configure_file"
          ),
          echo_cmd = T,
          wd = outdir
        )
      }
      adata <- read.delim(dapars_results_path)

      # system2(command = python,
      #         args = c(
      #           "./inst/script/dapars/DaPars_main.py",
      #           # file.path(script_dir, "DaPars_main.py"),
      #           dapars2_configure_file_path
      #         ))
      #   # summary results
      #   chr <- read.table(chr_list)[, 1]
      #   adata <- lapply(chr, function(x) {
      #     read.delim(file.path(
      #       outdir,
      #       paste0("results_", x),
      #       paste0("dapars2_result_temp.", x, ".txt")
      #     ))
      #   })
      #   adata <- Reduce(f = rbind, x = adata)
      #   sample_name_index <- 5:ncol(adata)
      #   if (is.null(sample_name)) {
      #     colnames(adata)[sample_name_index] <- x$metadata$sample
      #   } else {
      #     colnames(adata)[sample_name_index] <- sample_name
      #   }
      #
      gene_info <- str_split_fixed(adata$Gene, "\\|", 4)
      adata <- adata |>
        mutate(
          gene_id = gene_info[, 1],
          gene_symbol = gene_info[, 2],
          chr = gene_info[, 3],
          strand = gene_info[, 4]
        )
      return(adata)
  }


#' fastp: an ultra-fast all-in-one FASTQ preprocessor
#'
#' @param software
#' @param input_r1
#' @param output_r1
#' @param input_r2
#' @param output_r2
#' @param html
#' @param json
#' @param ncores
#' @param log
#' @param ...
#'
#' @return
#' @importFrom processx run
#' @importFrom logger log_info
#' @export
#'
#' @examples
fastp <- function(software = NULL,
                  input_r1,
                  output_r1,
                  input_r2 = NULL,
                  output_r2 = NULL,
                  html = NULL,
                  json = NULL,
                  ncores = 3,
                  log = NULL,
                  ...) {
  # Use Sun's Lab server default software location
  if (is.null(software)) {
    software <- "/opt/fastp/0.23.2/bin/fastp"
  }
  # check command path
  if (!file.exists(software)) {
    stop(paste0("Command path not found: ", software))
  }
  # Run external command, and wait until finishes
  log_info("Start running external command:")
  results <- run(
    command = software,
    args = c(
      "--in1",
      input_r1,
      "--out1",
      output_r1,
      pv("--in2", input_r2),
      pv("--out2", output_r2),
      pv("--html", html),
      pv("--json", json),
      "--thread",
      ncores
    ),
    echo_cmd = TRUE,
    spinner = TRUE,
    stderr_to_stdout = TRUE,
    stdout = log,
    ...
  )
  return(results)
}

#' Title
#'
#' @param software
#' @param input_r1
#' @param input_r2
#' @param index
#' @param sam
#' @param ncores
#' @param summary_file
#' @param log
#' @param ...
#'
#' @return
#' @importFrom processx run
#' @importFrom logger log_info
#' @export
#'
#' @examples
hisat2 <- function(software = NULL,
                   input_r1,
                   input_r2 = NULL,
                   index,
                   sam,
                   ncores = 5,
                   summary_file = NULL,
                   log = NULL,
                   ...) {
  # Use Sun's Lab server default software location
  if (is.null(software)) {
    software <- "/opt/hisat2/0.2.2/bin/hisat2"
  }
  # Check command path
  if (!file.exists(software)) {
    stop(paste0("Command path not found: ", software))
  }
  # Pre-processing
  if (is.null(input_r2)) {
    input <- pv("-U", input_r1)
  } else {
    input <- c("-1", input_r1, "-2", input_r2)
  }
  # Run external command, and wait until finishes
  log_info("Start running external command:")
  results <- run(
    command = software,
    args = c(
      "-x",
      index,
      input,
      "-S",
      sam,
      pv("--summary-file", summary_file),
      "--downstream-transcriptome-assembly",
      "--threads",
      ncores
    ),
    echo_cmd = TRUE,
    spinner = TRUE,
    stderr_to_stdout = TRUE,
    stdout = log,
    ...
  )
  return(results)
}

#' Title
#'
#' @param software
#' @param command
#' @param input
#' @param output
#' @param quality
#' @param log
#' @param ...
#'
#' @return
#' @importFrom processx run
#' @importFrom logger log_info
#' @export
#'
#' @examples
samtools <- function(software = NULL,
                     command,
                     input,
                     output = NULL,
                     quality = NULL,
                     log = NULL,
                     ...) {
  # Use Sun's Lab server default software location
  if (is.null(software)) {
    software <- "/opt/samtools/1.16.1/bin/samtools"
  }
  # check command path
  if (!file.exists(software)) {
    stop(paste0("Command path not found: ", software))
  }
  # Run external command, and wait until finishes
  log_info("Start running external command:")
  # view
  if (command == "view") {
    results <- run(
      command = software,
      args = c(
        command,
        input,
        pv("--min-MQ", quality),
        "--with-header",
        "-S",
        "--bam",
        "--output",
        output
      ),
      echo_cmd = TRUE,
      spinner = TRUE,
      stderr_to_stdout = TRUE,
      stdout = log,
      ...
    )
  }
  # sort
  if (command == "sort") {
    results <- run(
      command = software,
      args = c(command,
               input,
               "-o",
               output),
      echo_cmd = TRUE,
      spinner = TRUE,
      stderr_to_stdout = TRUE,
      stdout = log,
      ...
    )
  }
  # index
  if (command == "index") {
    results <- run(
      command = software,
      args = c(command,
               input),
      echo_cmd = TRUE,
      spinner = TRUE,
      stderr_to_stdout = TRUE,
      stdout = log,
      ...
    )
  }
  return(results)
}

#' Run Picard command on a BAM file
#'
#' This function provides a simple interface for running Picard commands on a BAM file.
#'
#' @param software The path to the Picard software. If not provided, it defaults to "/opt/picard/2.27.5/bin/picard".
#' @param command The Picard command to run, default is "MarkDuplicates".
#' @param input The input BAM file.
#' @param output The output BAM file.
#' @param metric_file The path to the file where Picard will write metrics about the command.
#' If not provided, no metrics will be written.
#' @param remove_duplicates A logical value indicating whether to remove duplicates from the BAM file.
#' The default is `TRUE`.
#' @param log The path to the file where Picard will write log output. If not provided, log output
#' will be printed to the console.
#' @param ... Additional arguments to pass to the `run()` function that actually executes the Picard command.
#'
#' @return The results of the Picard command.
#'
#' @examples
#' picard(input = "input.bam", output = "output.bam")
#'
#' @importFrom processx run
#' @importFrom logger log_info
#' @export
picard <- function(software = NULL,
                   command = "MarkDuplicates",
                   input,
                   output,
                   metric_file = NULL,
                   remove_duplicates = TRUE,
                   log = NULL,
                   ...) {
  # Use Sun's Lab server default software location
  if (is.null(software)) {
    software <- "/opt/picard/2.27.5/bin/picard"
  }
  # check command path
  if (!file.exists(software)) {
    stop(paste0("Command path not found: ", software))
  }
  remove_duplicates <- ifelse(remove_duplicates, "true", "false")
  # Run external command, and wait until finishes
  log_info("Start running external command:")
  results <- run(
    command = software,
    args = c(
      "--INPUT",
      input,
      "--OUTPUT",
      output,
      pv("--METRICS_FILE", metric_file),
      "--REMOVE_DUPLICATES",
      true
    ),
    echo_cmd = TRUE,
    spinner = TRUE,
    stderr_to_stdout = TRUE,
    stdout = log,
    ...
  )
  return(results)
}
#' Title
#'
#' @param software
#' @param gtf
#' @param input
#' @param output
#' @param log
#' @param paired
#' @param ncores
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
featurecounts <- function(software = NULL,
                          gtf,
                          input,
                          output,
                          log = NULL,
                          paired = TRUE,
                          ncores = ncores,
                          ...) {
  # Use Sun's Lab server default software location
  if (is.null(software)) {
    software <- "/opt/subread/2.0.3/bin/featureCounts"
  }
  # Check command path
  if (!file.exists(software)) {
    stop(paste0("Command path not found: ", software))
  }
  # Run external command, and wait until finishes
  log_info("Start running external command:")
  results <- run(
    command = software,
    args = c(
      "-a",
      gtf,
      "-o",
      output,
      input,
      "-g",
      "gene_name",
      "-T",
      ncores,
      pv("-p", paired)
    ),
    echo_cmd = TRUE,
    spinner = TRUE,
    stderr_to_stdout = TRUE,
    stdout = log,
    ...
  )
  return(results)
}

stringtie <- function(software = NULL,
                      gtf,
                      input,
                      output,
                      ncores = 5,
                      log = NULL,
                      ...) {
  # Use Sun's Lab server default software location
  if (is.null(software)) {
    software <- "/opt/stringtie/2.2.1/bin/stringtie"
  }
  # check command path
  if (!file.exists(software)) {
    stop(paste0("Command path not found: ", software))
  }
  # Run external command, and wait until finishes
  log_info("Start running external command:")
  results <- run(
    command = software,
    args = c("-e",
             "-B",
             "-p", ncores,
             "-G", gtf,
             "-o", output,
             input),
    echo_cmd = TRUE,
    spinner = TRUE,
    stderr_to_stdout = TRUE,
    stdout = log,
    ...
  )
  return(results)
}

#' Title
#'
#' @param x
#' @param software
#' @param ncores
#' @param workers
#'
#' @return
#' @importFrom furrr future_walk furrr_options
#' @importFrom future plan multisession
#' @importFrom logger log_info
#' @export
#'
#' @examples
hp_control_quality <- function(x,
                               software = NULL,
                               ncores = 3,
                               workers = "auto") {
  # Extract data
  method <- "fastp"
  # outdir
  outdir <- hp_set_path(file.path(x@outdir, "clean"))
  # board
  board <- x@board
  # sample name
  sample_name <- x$sample
  # input
  input_r1 <- extract_column(board, "input_r1")
  input_r2 <- extract_column(board, "input_r2")
  ip_r1 <- extract_column(board, "ip_r1")
  ip_r2 <- extract_column(board, "ip_r2")
  # output
  input_clean_r1 <-
    generate_output(board,
                    "input_r1",
                    outdir,
                    sample_name,
                    "input.R1",
                    method,
                    "fastq.gz")
  input_clean_r2 <-
    generate_output(board,
                    "input_r2",
                    outdir,
                    sample_name,
                    "input.R2",
                    method,
                    "fastq.gz")
  input_html <-
    generate_output(board,
                    "input_r1",
                    outdir,
                    sample_name,
                    "input",
                    method,
                    "html")
  input_json <-
    generate_output(board,
                    "input_r1",
                    outdir,
                    sample_name,
                    "input",
                    method,
                    "json")
  ip_clean_r1 <-
    generate_output(board,
                    "ip_r1",
                    outdir,
                    sample_name,
                    "ip.R1",
                    method,
                    "fastq.gz")
  ip_clean_r2 <-
    generate_output(board,
                    "ip_r2",
                    outdir,
                    sample_name,
                    "ip.R2",
                    method,
                    "fastq.gz")
  ip_html <-
    generate_output(board,
                    "ip_r1", outdir, sample_name, "ip", method, "html")
  ip_json <-
    generate_output(board,
                    "ip_r1", outdir, sample_name, "ip", method, "json")
  # Combine
  combine_r1 <- c(input_r1, ip_r1)
  combine_r2 <- c(input_r2, ip_r2)
  combine_clean_r1 <-
    c(input_clean_r1, ip_clean_r1)
  combine_clean_r2 <-
    c(input_clean_r2, ip_clean_r2)
  combine_html <- c(input_html, ip_html)
  combine_json <- c(input_json, ip_json)

  # Execute program
  need_to_run <- combine_r1[!file.exists(combine_clean_r1)]
  n <- length(need_to_run)
  if (n > 0) {
    if (workers == "auto" && n <= 4) {
      workers <- n
    } else if (workers == "auto" && n > 4) {
      workers <- 4
    }
    log_info("{n} samples start quality control...")
    log_info("Use {workers} workers, each assigned {ncores} cores...")
    plan(multisession, workers = workers)
  }
  # This does run in parallel!
  future_walk(
    .x = which(combine_r1 %in% need_to_run),
    .f = function(i) {
      if (!file.exists(combine_clean_r1[i])) {
        results <- fastp(
          software = software,
          input_r1 = combine_r1[i],
          output_r1 = combine_clean_r1[i],
          input_r2 = combine_r2[i],
          output_r2 = combine_clean_r2[i],
          html = combine_html[i],
          json = combine_json[i],
          ncores = ncores
        )
      }
    },
    .options = furrr_options(seed = TRUE)
  )
  # Return object
  board <-
    insert_column(board, "input_r1", "input_clean_r1", input_clean_r1)
  board <-
    insert_column(board, "input_r2", "input_clean_r2", input_clean_r2)
  board <-
    insert_column(board, "ip_r1", "ip_clean_r1", ip_clean_r1)
  board <-
    insert_column(board, "ip_r2", "ip_clean_r2", ip_clean_r2)
  x@board <- board
  return(x)
}



#' Title
#'
#' @param x
#' @param hisat2
#' @param samtools
#' @param species
#' @param index
#' @param quality
#' @param ncores
#' @param workers
#'
#' @return
#' @importFrom furrr future_walk furrr_options
#' @importFrom future plan multisession
#' @importFrom logger log_info
#' @export
#'
#' @examples
hp_align_read <- function(x,
                          hisat2 = NULL,
                          samtools = NULL,
                          picard = NULL,
                          species = "human",
                          index = NULL,
                          quality = 20,
                          remove_duplicates = FALSE,
                          ncores = 5,
                          workers = "auto") {
  # index
  if (is.null(index)) {
    if (species == "human") {
      index <- "/DATA/public/index/hisat2/homo_sapiens/GRCh38.p13/genome"
    } else if (species == "mouse") {
      index <- "/DATA/public/index/hisat2/mus_musculus/GRCm38.p6/genome"
    }
  }
  hp_check_file(paste0(index, ".1.ht2"))

  # Extract data
  method <- "hisat2"
  # outdir
  outdir <- hp_set_path(file.path(x@outdir, "align"))
  # board
  board <- x@board
  # sample name
  sample_name <- x$sample
  # input
  input_clean_r1 <- extract_column(board, "input_clean_r1")
  input_clean_r2 <- extract_column(board, "input_clean_r2")
  ip_clean_r1 <- extract_column(board, "ip_clean_r1")
  ip_clean_r2 <- extract_column(board, "ip_clean_r2")
  # output
  input_sam <-
    generate_output(board,
                    "input_clean_r1",
                    outdir,
                    sample_name,
                    "input",
                    method,
                    "sam")
  input_bam <-
    generate_output(board,
                    "input_clean_r1",
                    outdir,
                    sample_name,
                    "input",
                    method,
                    "bam")
  input_sorted_bam <-
    generate_output(board,
                    "input_clean_r1",
                    outdir,
                    sample_name,
                    "input",
                    method,
                    "sorted.bam")
  input_rmdup_sorted_bam <-
    generate_output(board,
                    "input_clean_r1",
                    outdir,
                    sample_name,
                    "input",
                    method,
                    "rmdup.sorted.bam")
  input_metric <-
    generate_output(board,
                    "input_clean_r1",
                    outdir,
                    sample_name,
                    "input",
                    method,
                    "metric.txt")
  input_summary <-
    generate_output(board,
                    "input_clean_r1",
                    outdir,
                    sample_name,
                    "input",
                    method,
                    "summary.txt")
  ip_sam <-
    generate_output(board,
                    "ip_clean_r1",
                    outdir,
                    sample_name,
                    "ip",
                    method,
                    "sam")
  ip_bam <-
    generate_output(board,
                    "ip_clean_r1",
                    outdir,
                    sample_name,
                    "ip",
                    method,
                    "bam")
  ip_sorted_bam <-
    generate_output(board,
                    "ip_clean_r1",
                    outdir,
                    sample_name,
                    "ip",
                    method,
                    "sorted.bam")
  ip_rmdup_sorted_bam <-
    generate_output(board,
                    "ip_clean_r1",
                    outdir,
                    sample_name,
                    "ip",
                    method,
                    "rmdup.sorted.bam")
  ip_metric <-
    generate_output(board,
                    "ip_clean_r1",
                    outdir,
                    sample_name,
                    "ip",
                    method,
                    "metric.txt")

  ip_summary <-
    generate_output(board,
                    "ip_clean_r1",
                    outdir,
                    sample_name,
                    "ip",
                    method,
                    "summary.txt")
  # Combine
  combine_clean_r1 <-
    c(input_clean_r1, ip_clean_r1)
  combine_clean_r2 <-
    c(input_clean_r2, ip_clean_r2)
  combine_sam <- c(input_sam, ip_sam)
  combine_summary <- c(input_summary, ip_summary)
  combine_bam <- c(input_bam, ip_bam)
  combine_sorted_bam <- c(input_sorted_bam, ip_sorted_bam)
  combine_rmdup_sorted_bam <-
    c(input_rmdup_sorted_bam, ip_rmdup_sorted_bam)
  combine_metric <- c(input_metric, ip_metric)
  # print(combine_sorted_bam)
  # Execute program
  if (remove_duplicates) {
    need_to_run <-
      combine_clean_r1[!file.exists(combine_rmdup_sorted_bam)]
  } else {
    need_to_run <- combine_clean_r1[!file.exists(combine_sorted_bam)]
  }
  n <- length(need_to_run)
  if (n > 0) {
    if (workers == "auto" && n <= 4) {
      workers <- n
    } else if (workers == "auto" && n > 4) {
      workers <- 4
    }
    log_info("{n} samples start aligning to genomes...")
    log_info("Use {workers} workers, each assigned {ncores} cores...")
    plan(multisession, workers = workers)
  }
  # This does run in parallel!
  future_walk(
    .x = which(combine_clean_r1 %in% need_to_run),
    .f = function(i) {
      # align
      if (!file.exists(combine_sam[i])) {
        results <- hisat2(
          software = hisat2,
          input_r1 = combine_clean_r1[i],
          input_r2 = combine_clean_r2[i],
          index = index,
          sam = combine_sam[i],
          ncores = ncores,
          summary_file = combine_summary[i]
        )
      }
      # sam to bam
      # if (!file.exists(combine_bam[i])) {
      #   results <- samtools(
      #     software = samtools,
      #     command = "view",
      #     input = combine_sam[i],
      #     output = combine_bam[i],
      #     quality = quality
      #   )
      # }
      # bam sorted
      if (!file.exists(combine_sorted_bam[i])) {
        results <- samtools(
          software = samtools,
          command = "sort",
          input = combine_sam[i],
          output = combine_sorted_bam[i]
        )
        results <- samtools(software = samtools,
                            command = "index",
                            input = combine_sorted_bam[i],)
      }
      # rmdup
      print(combine_rmdup_sorted_bam[i])
      if (!file.exists(combine_rmdup_sorted_bam[i]) &&
          remove_duplicates) {
        results <- picard(
          software = picard,
          input = combine_sorted_bam[i],
          output = combine_rmdup_sorted_bam[i],
          metric_file = combine_metric[i]
        )
        results <- samtools(software = samtools,
                            command = "index",
                            input = combine_rmdup_sorted_bam[i])
      }
    },
    .options = furrr_options(seed = TRUE)
  )
  # Return object
  # mapping rate
  if ("input_clean_r1" %in% colnames(x@board)) {
    input_mapping_rate <- sapply(input_summary, function(i) {
      mapping_rate <-
        read.table(file = i,
                   skip = 14)[, 1]
      mapping_rate <- as.numeric(gsub("%", "", mapping_rate))
      return(mapping_rate)
    })
    x$input_mapping_rate <- input_mapping_rate
  }
  if ("ip_clean_r1" %in% colnames(x@board)) {
    ip_mapping_rate <- sapply(ip_summary, function(i) {
      mapping_rate <-
        read.table(file = i,
                   skip = 14)[, 1]
      mapping_rate <- as.numeric(gsub("%", "", mapping_rate))
      return(mapping_rate)
    })
    x$ip_mapping_rate <- ip_mapping_rate
  }

  board <-
    insert_column(board, "input_clean_r1", "input_sam", input_sam)
  board <-
    insert_column(board, "input_clean_r1", "input_bam", input_bam)
  board <-
    insert_column(board,
                  "input_clean_r1",
                  "input_sorted_bam",
                  input_sorted_bam)
  board <-
    insert_column(board,
                  "input_clean_r1",
                  "input_rmdup_sorted_bam",
                  input_rmdup_sorted_bam)
  board <-
    insert_column(board, "ip_clean_r1", "ip_sam", ip_sam)
  board <-
    insert_column(board, "ip_clean_r1", "ip_bam", ip_bam)
  board <-
    insert_column(board, "ip_clean_r1", "ip_sorted_bam", ip_sorted_bam)
  board <-
    insert_column(board,
                  "ip_clean_r1",
                  "ip_rmdup_sorted_bam",
                  ip_rmdup_sorted_bam)
  x@board <- board
  return(x)
}


#' Title
#'
#' @param x
#' @param featurecounts
#' @param ncores
#' @param species
#' @param gtf
#'
#' @return
#' @export
#'
#' @examples
hp_quantify_gene <- function(x,
                             featurecounts = NULL,
                             ncores = 8,
                             species = "human",
                             gtf = NULL) {
  if (is.null(gtf)) {
    if (species == "human") {
      gtf <-
        "/DATA/public/genome/homo_sapiens/GRCh38.p13/gencode.v43.annotation.gtf"
    } else if (species == "mouse") {
      gtf <-
        "/DATA/public/genome/mus_musculus/GRCm38.p6/gencode.vM25.annotation.gtf"
    }
  }
  hp_check_file(gtf)
  # Extract data
  method <- "featurecounts"
  # outdir
  outdir <- hp_set_path(file.path(x@outdir, "counts"))
  # board
  board <- x@board
  # sample name
  sample_name <- x$sample
  # input
  input_bam <- extract_column(board, "input_sorted_bam")
  if (all(c("input_r1", "input_r2") %in% colnames(board))) {
    paired <- TRUE
  } else {
    paired <- FALSE
  }
  # output
  read_counts <-
    file.path(outdir, paste0("counts.", method, ".txt"))
  log <- file.path(outdir, paste0("counts.", method, ".log"))
  # Execute program
  n <- length(input_bam)
  log_info("{n} samples start quantify gene...")
  log_info("Use {ncores} cores...")

  # This does run in parallel!
  if (!file.exists(read_counts)) {
    results <- featurecounts(
      software = featurecounts,
      gtf = gtf,
      input = input_bam,
      output = read_counts,
      paired = paired,
      ncores = ncores,
      log = log
    )
  }
  # Return object
  log_adata <- read.delim(log)
  start <-
    str_locate_all(log_adata, "\\(")[[1]][-c(1, 2, 3),][, 2] + 1
  end <- str_locate_all(log_adata, "%\\)")[[1]][, 1] - 1
  assigned_rate <- sapply(seq_along(start), function(i) {
    as.numeric(str_sub(log_adata, start = start[i], end = end[i]))
  })
  x$assigned_rate <- assigned_rate
  adata <- read.table(read_counts, header = TRUE)
  counts <- adata |>
    select(Geneid, ends_with("bam")) |>
    column_to_rownames("Geneid")
  metadata <- adata |>
    select(-ends_with("bam")) |>
    column_to_rownames("Geneid")

  colnames(counts) <- colnames(x)
  counts <- counts[rowSums(counts) > 0,]
  metadata <- metadata[rownames(counts),]
  identical(rownames(counts), rownames(metadata))
  x[["RNA"]] <-
    CreateAssayObject(counts = counts)
  x@feature[["RNA"]] <- metadata[rownames(counts),]
  return(x)
}


#' Title
#'
#' @param x
#' @param control
#' @param treatment
#' @param group
#' @param log2fc_threshold
#' @param p_value_threshold
#' @param fdr_threshold
#' @param method
#'
#' @return
#' @importFrom tibble column_to_rownames
#' @export
#'
#' @examples
hp_compare_difference <- function(x,
                                  control,
                                  treatment,
                                  group = "group",
                                  log2fc_threshold = 0.5,
                                  p_value_threshold = 0.05,
                                  fdr_threshold = NULL,
                                  method = "DESeq2") {
  outdir <- hp_set_path(file.path(x@outdir, "deg"))
  raw_x <- x
  x$case <- x[[group]]
  x <- x |>
    subset(case %in% c(control, treatment))
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
      control, ".dds.rds"
    )))
    # 差异基因
    res <-
      DESeq2::results(dds, contrast = c("group", treatment, control))

    res_ordered <- res[order(res$pvalue),]
    deg <- na.omit(as.data.frame(res_ordered))
    if (!is.null(fdr_threshold)) {
      deg <- deg |>
        mutate(change = ifelse(
          padj <= fdr_threshold &
            abs(log2FoldChange) >= log2fc_threshold,
          ifelse(
            log2FoldChange > log2fc_threshold,
            "Upregulated",
            "Downregulated"
          ),
          "Stable"
        ))
    } else {
      deg <- deg |>
        mutate(
          change = ifelse(
            pvalue <= p_value_threshold &
              abs(log2FoldChange) >= log2fc_threshold,
            ifelse(
              log2FoldChange > log2fc_threshold,
              "Upregulated",
              "Downregulated"
            ),
            "Stable"
          )
        )
    }
    write.csv(deg,
              file.path(
                outdir,
                paste0(treatment,
                       "_vs_",
                       control, ".deg.", method, ".csv")
              ),
              row.names = TRUE,
              quote = FALSE)
    features <- x@feature$RNA
    features <- merge(features, deg, by = 0)
    features <- features |>
      column_to_rownames("Row.names")
    features <- features[rownames(raw_x),]
    identical(rownames(features), rownames(raw_x))
    raw_x@feature$RNA <- features
  }
  return(raw_x)
}



# detect_significant_APA <- function(data,
#                                    control_cols,
#                                    treatment_cols,
#                                    fdr_threshold = 0.05,
#                                    mean_diff_threshold = 0.2,
#                                    mean_fold_change_threshold = 1.5) {
#   # Calculate PDUIs for control and treatment samples
#   control_pd_uis <- rowMeans(data[, control_cols])
#   treatment_pd_uis <- rowMeans(data[, treatment_cols])
#
#   # Calculate PDUI differences and fold-changes
#   pd_diff <- treatment_pd_uis - control_pd_uis
#   pd_fold_change <- treatment_pd_uis / control_pd_uis
#
#   # Perform Fisher's exact test
#   fisher_p_values <- mapply(function(control_row, treatment_row) {
#     matrix_data <- matrix(c(control_row["wiL"],
#                             treatment_row["wiL"],
#                             control_row["wiS"],
#                             treatment_row["wiS"]), nrow = 2)
#     fisher.test(matrix_data)$p.value
#   }, data[, control_cols], data[, treatment_cols], SIMPLIFY = TRUE)
#
#   # Adjust p-values using Benjamini-Hochberg procedure
#   adjusted_p_values <- p.adjust(fisher_p_values, method = "BH")
#
#   # Apply criteria
#   significant_APA <- adjusted_p_values <= fdr_threshold &
#     pd_diff >= mean_diff_threshold &
#     pd_fold_change >= mean_fold_change_threshold
#
#   # Create output data frame
#   output_data <- data.frame(
#     gene = rownames(data),
#     control_pd_ui = control_pd_uis,
#     treatment_pd_ui = treatment_pd_uis,
#     fold_change = pd_fold_change,
#     p_value = fisher_p_values,
#     adjusted_p_value = adjusted_p_values
#   )
#
#   # Return significant APA events
#   return(output_data[significant_APA, ])
# }
