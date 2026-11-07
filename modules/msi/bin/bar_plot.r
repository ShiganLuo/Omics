#' Create a bar plot of hUSI scores and save as JPEG.
#'
#' @param vec Named numeric vector of hUSI scores.
#' @param outpng Path to save the output JPEG file.
#' @param base_width_mm Base device width in millimeters; grows with sample count.
#' @param base_height_mm Base device height in millimeters; grows with label length.
#' @param width_per_sample Extra width per bar to keep labels visible when many samples.
#' @param height_per_label Extra height per character of the longest label.
#' @param xlab_line Line offset for the x-axis title to avoid collisions with tick labels.
bar <- function(vec, outpng,
        base_width_mm = 150,
        base_height_mm = 120,
        width_per_sample = 20,
        height_per_label = 4,
        xlab_line = 3,
        xlab = "sample",
        ylab = "hUSI"
    ) {
    n <- length(vec)
    labels <- names(vec)
    max_label_len <- ifelse(length(labels) > 0, max(nchar(labels)), 5)

    # Scale device to fit labels: width grows with sample count, height with label length.
    width_mm <- max(base_width_mm, width_per_sample * n)
    height_mm <- max(base_height_mm, 80 + height_per_label * max_label_len)

    jpeg(width = width_mm, height = height_mm, units = "mm", res = 300, file = outpng)
    bottom_mar <- max(5, 2 + max_label_len * 0.3)  # add margin if labels are long
    par(mar = c(bottom_mar, 4, 2, 1))
    colors <- grDevices::rainbow(max(3, n))[seq_len(n)]
        barplot(vec, col = colors, xlab = "", ylab = ylab, cex.lab = 1.1, border = NA,
            las = 1, cex.names = min(0.9, 10 / max_label_len))
    mtext(xlab, side = 1, line = xlab_line, cex = 0.9)
    dev.off()
}

if (sys.nframe() == 0) {
    if (!requireNamespace("argparse", quietly = TRUE)) {
        stop("Package 'argparse' is required. Install with install.packages('argparse').")
    }
    suppressMessages(library(argparse))

    parser <- ArgumentParser(description = "Plot scores as a bar chart")
    parser$add_argument("-i","--input", required = TRUE, help = "Input table with sample and score columns")
    parser$add_argument("-o","--output", required = TRUE, help = "Output image path (jpg)")
    parser$add_argument("-s","--sample-col", dest = "sample_col", default = "sample", help = "Column name for sample IDs")
    parser$add_argument("-c","--score-col", dest = "score_col", default = "score", help = "Column name for scores")
    parser$add_argument("-t","--sep", default = "\t", help = "Field separator, default tab")
    parser$add_argument("-x","--xlab", default = "sample", help = "X-axis label")
    parser$add_argument("-y","--ylab", default = "MSI", help = "Y-axis label")

    args <- parser$parse_args()

    if (!file.exists(args$input)) {
        stop("Input file not found: ", args$input)
    }

    df <- read.table(args$input, header = TRUE, sep = args$sep, check.names = FALSE, stringsAsFactors = FALSE)
    if (!all(c(args$sample_col, args$score_col) %in% colnames(df))) {
        stop("Input file must contain columns ", args$sample_col, " and ", args$score_col)
    }

    vec <- df[[args$score_col]]
    names(vec) <- df[[args$sample_col]]
    bar(vec, outpng = args$output, xlab = args$xlab, ylab = args$ylab)
}
