# Valid output image formats (R base-graphics device set). pgf/raw/rgba are
# matplotlib-only and fail fast here (modules.md #9).
VALID_FORMATS <- c("png", "pdf", "svg", "ps", "eps", "tif", "tiff", "jpg", "jpeg")

#' Open a graphics device for one output format (device-per-extension helper).
#'
#' Sizing is in inches for every format; raster devices share `res`, so the
#' physical geometry is identical across devices.
#'
#' @param path Output file path (with extension).
#' @param fmt One of VALID_FORMATS.
#' @param width_in Device width in inches.
#' @param height_in Device height in inches.
#' @param res Resolution for raster devices. Default 300.
open_device <- function(path, fmt, width_in, height_in, res = 300) {
    switch(fmt,
        png = grDevices::png(path, width = width_in, height = height_in,
                             units = "in", res = res),
        pdf = grDevices::pdf(path, width = width_in, height = height_in),
        svg = grDevices::svg(path, width = width_in, height = height_in),
        ps  = grDevices::postscript(path, width = width_in, height = height_in),
        eps = grDevices::postscript(path, width = width_in, height = height_in,
                                    onefile = FALSE, horizontal = FALSE),
        tif =, tiff = grDevices::tiff(path, width = width_in, height = height_in,
                                      units = "in", res = res),
        jpg =, jpeg = grDevices::jpeg(path, width = width_in, height = height_in,
                                      units = "in", res = res),
        stop(sprintf("Unsupported image format(s): %s; valid: %s", fmt,
                     paste(VALID_FORMATS, collapse = " "))))
}

#' Create a bar plot of hUSI scores and save it in each requested format.
#'
#' @param vec Named numeric vector of hUSI scores.
#' @param out_base Output path prefix WITHOUT extension; one file is written
#'   per requested format as `<out_base>.<fmt>`.
#' @param formats Character vector of output formats (png/pdf/svg/ps/eps/tif/
#'   tiff/jpg/jpeg). Default "jpg" (historical single-JPEG behavior).
#' @param base_width_mm Base device width in millimeters; grows with sample count.
#' @param base_height_mm Base device height in millimeters; grows with label length.
#' @param width_per_sample Extra width per bar to keep labels visible when many samples.
#' @param height_per_label Extra height per character of the longest label.
#' @param xlab_line Line offset for the x-axis title to avoid collisions with tick labels.
bar <- function(vec, out_base,
        formats = "jpg",
        base_width_mm = 150,
        base_height_mm = 120,
        width_per_sample = 20,
        height_per_label = 4,
        xlab_line = 3,
        xlab = "sample",
        ylab = "hUSI"
    ) {
    formats <- as.character(formats)
    bad <- setdiff(formats, VALID_FORMATS)
    if (length(bad)) {
        stop(sprintf("Unsupported image format(s): %s; valid: %s",
                     paste(bad, collapse = ", "), paste(VALID_FORMATS, collapse = " ")))
    }

    n <- length(vec)
    labels <- names(vec)
    max_label_len <- ifelse(length(labels) > 0, max(nchar(labels)), 5)

    width_mm <- max(base_width_mm, width_per_sample * n)
    height_mm <- max(base_height_mm, 80 + height_per_label * max_label_len)
    width_in <- width_mm / 25.4
    height_in <- height_mm / 25.4

    for (fmt in formats) {
        outfile <- sprintf("%s.%s", out_base, fmt)
        open_device(outfile, fmt, width_in, height_in, res = 300)
        bottom_mar <- max(5, 2 + max_label_len * 0.3)
        par(mar = c(bottom_mar, 4, 2, 1))
        colors <- grDevices::rainbow(max(3, n))[seq_len(n)]
            barplot(vec, col = colors, xlab = "", ylab = ylab, cex.lab = 1.1, border = NA,
                las = 1, cex.names = min(0.9, 10 / max_label_len))
        mtext(xlab, side = 1, line = xlab_line, cex = 0.9)
        dev.off()
    }
}

if (sys.nframe() == 0) {
    if (!requireNamespace("argparse", quietly = TRUE)) {
        stop("Package 'argparse' is required. Install with install.packages('argparse').")
    }
    suppressMessages(library(argparse))

    parser <- ArgumentParser(description = "Plot scores as a bar chart")
    parser$add_argument("-i","--input", required = TRUE, help = "Input table with sample and score columns")
    parser$add_argument("-o","--output", required = TRUE,
                        help = "Output path prefix without extension (files written as <prefix>.<fmt>); a legacy full filename like out.jpg is accepted, its extension is dropped")
    parser$add_argument("-s","--sample-col", dest = "sample_col", default = "sample", help = "Column name for sample IDs")
    parser$add_argument("-c","--score-col", dest = "score_col", default = "score", help = "Column name for scores")
    parser$add_argument("-t","--sep", default = "\t", help = "Field separator, default tab")
    parser$add_argument("-x","--xlab", default = "sample", help = "X-axis label")
    parser$add_argument("-y","--ylab", default = "MSI", help = "Y-axis label")
    parser$add_argument("-f","--format", action = "append", dest = "formats",
                        help = "Output image format, repeatable (default: jpg); valid: png pdf svg ps eps tif tiff jpg jpeg")

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
    out_base <- sub("\\.(png|pdf|svg|ps|eps|tif|tiff|jpg|jpeg)$", "", args$output)
    formats <- if (is.null(args$formats)) "jpg" else as.character(unlist(args$formats))
    bar(vec, out_base = out_base, formats = formats, xlab = args$xlab, ylab = args$ylab)
}
