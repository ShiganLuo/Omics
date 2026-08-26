#' Call peaks from alignment results
#'
#' @param x
#' @param macs2
#' @param species
#' @param format
#' @param genome_size
#' @param keep_dup
#' @param output
#' @param prefix
#' @param qvalue
#' @param pvalue
#' @param fe_cutoff
#' @param seed
#'
#' @return
#' @importFrom furrr future_walk
#' @importFrom future plan multisession
#' @export
#'
#' @examples
cat_macs2 <-
  function(x,
           macs2 = NULL,
           bedtools = NULL,
           species = "human",
           format = "BEDPE",
           genome_size = NULL,
           keep_dup = "all",
           output = "./macs2",
           prefix = NULL,
           qvalue = 0.05,
           slocal = 1000,
           pvalue = NULL,
           fe_cutoff = 1,
           workers = "auto") {
    if (is.null(macs2)) {
      macs2 <- "/opt/macs2/2.2.7.1/bin/macs2"
      hp_check_file(macs2)
    }
    if (is.null(bedtools)) {
      bedtools <- "/opt/bedtools/2.30.0/bin/bedtools"
      hp_check_file(bedtools)
    }

    if (species == "human") {
      genome_size <- "hs"
    } else if (species == "mouse") {
      genome_size <- "mm"
    }
    # outdir
    outdir <- hp_set_path(file.path(x@board, "macs2"))
    # sample name
    sample_name <- x$sample
    # bam
    input_bam <- x$input_bam
    ip_bam <- x$ip_bam
    # bed
    input_bed <-
      gsub(
        pattern = ".bam",
        replacement = ".bed",
        x = input_bam
      )
    ip_bed <- gsub(
      pattern = ".bam",
      replacement = ".bed",
      x = ip_bam
    )
    # narrow peak
    narrow_peak <-
      file.path(outdir, paste0(sample_name, "_peaks.narrowPeak"))
    # bed
    bed <-
      file.path(outdir, paste0(sample_name, ".narrow_peak.bed"))
    # Set a "plan" for how the code should run.
    n <- length(bed[!file.exists(bed)])
    if (n > 0) {
      if (workers == "auto" && n <= 4) {
        workers <- n
      } else if (workers == "auto" && n > 4) {
        workers <- 4
      }
      log_info("Use {workers} workers...")
      plan(multisession, workers = workers)
    }
    # This does run in parallel!
    future_walk(
      .x = seq_along(sample_name),
      .f = function(i) {
        if (!file.exists(bed[i])) {
          # bam to bed
          system2(
            command = bedtools,
            args = c(
              "bamtobed",
              "-i",
              input_bam[i],
              "|",
              "grep",
              "'^chr'",
              ">",
              input_bed[i]
            )
          )
          system2(
            command = bedtools,
            args = c(
              "bamtobed",
              "-i",
              ip_bam[i],
              "|",
              "grep",
              "'^chr'",
              ">",
              ip_bed[i]
            )
          )
          # macs2
          system2(
            command = macs2,
            args = c(
              "callpeak",
              "--treatment",
              ip_bed[i],
              "--control",
              input_bed[i],
              "--format",
              # {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,
              # ELANDEXPORT,BOWTIE,BAMPE,BEDPE}
              format,
              "--gsize",
              genome_size,
              "--nomodel",
              "--slocal",
              slocal,
              # Effective genome size. It can be 1.0e+9 or 1000000000, or
              # shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9),
              # 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8)
              "--keep-dup",
              keep_dup,
              # --keep-dup: 保留重复。默认MACS(auto)会使用二项分布估计每个位置上
              # 是否存在重复（默认是1，也就是每个位置上出现一个read的概率最大）。
              # 如果你前期已经去重，那就使用all省了这一步.
              "--outdir",
              outdir,
              "--name",
              sample_name[i],
              # "--bdg",
              # "--trackline",
              # "--SPMR",
              # "--nomodel",
              # "--extsize",
              # extsize,
              if (is.null(pvalue)) {
                c(
                  "--qvalue",
                  qvalue
                )
              },
              if (!is.null(pvalue)) {
                c(
                  "--pvalue",
                  pvalue
                )
              }
              # "--fe-cutoff",
              # fe_cutoff
            )
            # stdout = file.path(outdir, paste0(sample_name[i], ".out")),
            # stderr = file.path(outdir, paste0(sample_name[i], ".err"))
          )
          # narrow_peak_bed <-
          #   read.table(narrow_peak[i])[, c(4, 1, 2, 3, 6)]
          narrow_peak_bed <-
            read.table(narrow_peak[i])[, c(1, 2, 3)]
          narrow_peak_bed$V4 <- "+"
          # narrow_peak_bed <-
          #   narrow_peak_bed[grep("chr", narrow_peak_bed[, 1]),]
          write.table(
            x = narrow_peak_bed,
            file = bed[i],
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t"
          )
        }
      }
    )
    x$bed <- bed
    return(x)
  }


#' Program will find de novo and known motifs in regions in the genome
#'
#' @param x
#' @param find_motifs_genome
#' @param shufflebed
#' @param bg
#' @param bed
#' @param genome
#' @param outdir
#' @param len
#' @param ncores
#'
#' @return
#' @importFrom furrr future_walk
#' @importFrom future plan multisession
#' @export
#'
#' @examples
cat_find_motifs_genome <- function(x,
                                   find_motifs_genome = NULL,
                                   shufflebed = NULL,
                                   bg = FALSE,
                                   bed = NULL,
                                   species = "human",
                                   genome = NULL,
                                   outdir = "./homer",
                                   rounds = 10,
                                   len = c(6, 7),
                                   ncores = 10,
                                   workers = "auto") {
  if (is.null(find_motifs_genome)) {
    find_motifs_genome <- "/opt/homer/4.11/bin/findMotifsGenome.pl"
    hp_check_file(find_motifs_genome)
  }

  if (species == "human") {
    genome <- "hg38"
  } else if (species == "mouse") {
    genome <- "mm10"
  }

  sample_name <- x$sample
  bed <- x$bed
  outdir <- file.path(x@board, "homer", sample_name)
  homer_results <- file.path(outdir, "homerResults.html")

  # Set a "plan" for how the code should run.
  # n <- length(sample_name)
  n <- length(homer_results[!file.exists(homer_results)])
  if (n > 0) {
    if (workers == "auto" && n <= 4) {
      workers <- n
    } else if (workers == "auto" && n > 4) {
      workers <- 4
    }
    log_info("Use {workers} workers, each assigned {ncores} cores...")
    plan(multisession, workers = workers)
  }
  # This does run in parallel!
  future_walk(
    .x = seq_along(sample_name),
    .f = function(i) {
      if (!file.exists(homer_results[i])) {
        system2(
          command = find_motifs_genome,
          args = c(
            bed[i],
            genome,
            outdir[i],
            "-len",
            paste(len, collapse = ","),
            "-rna",
            "-p",
            ncores
          ),
          env = c(paste0(
            "PATH=", dirname(find_motifs_genome), ":$PATH"
          ))
        )
      }
    }
  )
  x$homer <- homer_results
  return(x)
}

#### cat_annotate_peaks ----
cat_annotate_peaks <- function(annotate_peaks,
                               peak_file,
                               genome) {
  bed_annotation <-
    str_replace_all(peak_file, ".bed", "_annotation.bed")

  system2(
    command = annotate_peaks,
    args = c(
      peak_file,
      genome,
      ">",
      bed_annotation
    )
  )
  return(bed_annotation)
}


cat_guitar <- function(x, group = "group") {
  bed <- x$bed
  group <- x[[group]]

  library(Guitar)
  # importdifferentdataformatsintoanamedlistobject.
  # Thesegenomicfeaturesareusingmm10genomeassembly
  stBedFiles <-
    list(
      system.file(
        "extdata",
        "m6A_mm10_exomePeak_1000peaks_bed12.bed",
        package = "Guitar"
      ),
      system.file(
        "extdata",
        "m6A_mm10_exomePeak_1000peaks_bed6.bed",
        package = "Guitar"
      )
    )
  # BuildGuitarCoordinates
  txdb_file <- system.file("extdata", "mm10_toy.sqlite",
    package = "Guitar"
  )
  txdb <- loadDb(txdb_file)
  # GuitarPlot
  GuitarPlot(
    txTxdb = txdb,
    stBedFiles = stBedFiles,
    headOrtail = TRUE,
    enableCI = FALSE,
    mapFilterTranscript = TRUE,
    pltTxType = c("mrna"),
    stGroupName = c("BED12", "BED6")
  )
}
