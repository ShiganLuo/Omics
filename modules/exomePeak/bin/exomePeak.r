# Load required libraries
library(argparse)
library(exomePeak)

options(width = 60)

run_peak_calling <- function(gtf, ip_bams, input_bams, outdir = "./") {
  result <- exomepeak(GENE_ANNO_GTF = gtf,
                      IP_BAM = ip_bams,
                      INPUT_BAM = input_bams,
                      OUT_DIR = outdir)
  print(names(result))
  result
}

run_diff_analysis <- function(gtf, ip_bams, input_bams, treated_ip_bams, treated_input_bams, outdir = "./") {
  result <- exomepeak(GENE_ANNO_GTF = gtf,
                      IP_BAM = ip_bams,
                      INPUT_BAM = input_bams,
                      TREATED_IP_BAM = treated_ip_bams,
                      TREATED_INPUT_BAM = treated_input_bams,
                      OUT_DIR = outdir)
  print(names(result))
  print(is.na(result$con_sig_diff_peaks))
  result
}


main <- function() {
  parser <- ArgumentParser(description = "exomePeak peak calling and differential analysis")
  parser$add_argument("--gtf", required = TRUE, help = "Path to the GTF annotation file")
  parser$add_argument("--ip_bams", nargs = "+", required = TRUE, help = "List of IP BAM files")
  parser$add_argument("--input_bams", nargs = "+", required = TRUE, help = "List of Input BAM files")
  parser$add_argument("--treated_ip_bams", nargs = "*", help = "List of treated IP BAM files")
  parser$add_argument("--treated_input_bams", nargs = "*", help = "List of treated Input BAM files")
  parser$add_argument("--outdir", default = "", help = "Output directory for results")

  args <- parser$parse_args()


  has_treated <- length(args$treated_ip_bams) > 0 && length(args$treated_input_bams) > 0
  if (has_treated) {
    diff_result <- run_diff_analysis(gtf = args$gtf,
                                     ip_bams = args$ip_bams,
                                     input_bams = args$input_bams,
                                     treated_ip_bams = args$treated_ip_bams,
                                     treated_input_bams = args$treated_input_bams,
                                     outdir = args$outdir
                                     )
  } else {
    peak_result <- run_peak_calling(gtf = args$gtf,
                                  ip_bams = args$ip_bams,
                                  input_bams = args$input_bams,
                                  outdir = args$outdir)
  }
}

if (interactive() || sys.nframe() == 0) {
  main()
}

