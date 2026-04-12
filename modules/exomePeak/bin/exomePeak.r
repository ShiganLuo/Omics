# Load required libraries
library("argparse")
library("exomePeak")

options(width = 60)

run_peak_calling <- function(gtf, ip_bams, input_bams) {
  result <- exomepeak(GENE_ANNO_GTF = gtf,
                      IP_BAM = ip_bams,
                      INPUT_BAM = input_bams)
  print(names(result))
  result
}

run_diff_analysis <- function(gtf, ip_bams, input_bams, treated_ip_bams, treated_input_bams) {
  result <- exomepeak(GENE_ANNO_GTF = gtf,
                      IP_BAM = ip_bams,
                      INPUT_BAM = input_bams,
                      TREATED_IP_BAM = treated_ip_bams,
                      TREATED_INPUT_BAM = treated_input_bams)
  print(names(result))
  print(is.na(result$con_sig_diff_peaks))
  result
}

write_peak_outputs <- function(result, outprefix) {
  if (is.null(outprefix) || outprefix == "") {
    return(invisible(NULL))
  }
  out_dir <- dirname(outprefix)
  if (out_dir != "." && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  con_peaks <- result$con_peaks
  all_peaks <- result$all_peaks

  con_df <- as.data.frame(mcols(con_peaks))
  all_df <- as.data.frame(mcols(all_peaks))

  write.table(con_df,
              file = paste0(outprefix, "peak_consistent_peaks.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(all_df,
              file = paste0(outprefix, "peak_all_peaks.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  summary_lines <- c(
    paste0("names(result): ", paste(names(result), collapse = ","))
  )
  writeLines(summary_lines, con = paste0(outprefix, "peak_summary.txt"))
}

write_diff_outputs <- function(result, outprefix) {
  if (is.null(outprefix) || outprefix == "") {
    return(invisible(NULL))
  }
  out_dir <- dirname(outprefix)
  if (out_dir != "." && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }


  if (!is.null(result$diff_peaks)) {
    diff_df <- as.data.frame(mcols(result$diff_peaks))
    write.table(diff_df,
                file = paste0(outprefix, "diff_peaks.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    lines <- c(
      "No significant differential peaks identified with the specified thresholds."
    )
    writeLines(lines, con = paste0(outprefix, "diff_peaks.tsv"))
  }

  if (!is.null(result$sig_siff_peaks)) {
    sig_siff_df <- as.data.frame(mcols(result$sig_siff_peaks))
    write.table(sig_siff_df,
                file = paste0(outprefix, "sig_siff_peaks.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    lines <- c(
      "No significant differential peaks identified with the specified thresholds."
    )
    writeLines(lines, con = paste0(outprefix, "sig_siff_peaks.tsv"))
  }

  if (!is.null(result$con_sig_diff_peaks)) {
    con_sig_diff_df <- as.data.frame(mcols(result$con_sig_diff_peaks))
    write.table(con_sig_diff_df,
                file = paste0(outprefix, "con_sig_diff_peaks.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    lines <- c(
      "No significant differential peaks identified among the consistent peaks with the specified thresholds."
    )
    writeLines(lines, con = paste0(outprefix, "con_sig_diff_peaks.tsv"))
  }

}

main <- function() {
  parser <- ArgumentParser(description = "exomePeak peak calling and differential analysis")
  parser$add_argument("--gtf", required = TRUE, help = "Path to the GTF annotation file")
  parser$add_argument("--ip_bams", nargs = "+", required = TRUE, help = "List of IP BAM files")
  parser$add_argument("--input_bams", nargs = "+", required = TRUE, help = "List of Input BAM files")
  parser$add_argument("--treated_ip_bams", nargs = "*", help = "List of treated IP BAM files")
  parser$add_argument("--treated_input_bams", nargs = "*", help = "List of treated Input BAM files")
  parser$add_argument("--outprefix", default = "", help = "Output file prefix (can include directory)")

  args <- parser$parse_args()


  has_treated <- length(args$treated_ip_bams) > 0 && length(args$treated_input_bams) > 0
  if (has_treated) {
    diff_result <- run_diff_analysis(gtf = args$gtf,
                                     ip_bams = args$ip_bams,
                                     input_bams = args$input_bams,
                                     treated_ip_bams = args$treated_ip_bams,
                                     treated_input_bams = args$treated_input_bams)
    write_diff_outputs(diff_result, args$outprefix)
  } else {
    peak_result <- run_peak_calling(gtf = args$gtf,
                                  ip_bams = args$ip_bams,
                                  input_bams = args$input_bams)
    write_peak_outputs(peak_result, args$outprefix)
  }
}

if (interactive() || sys.nframe() == 0) {
  main()
}

