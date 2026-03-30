#' @export
log_msg <- function(level = c("INFO","WARN","ERROR"), ..., quit = FALSE) {
  level <- match.arg(level)
  msg <- paste(...)
  prefix <- paste0(
    "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
    "[", level, "] "
  )
  message(prefix, msg)
  if (quit) stop(msg, call. = FALSE)
}

load_script_source <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    script_dir <- if(length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
    rfiles <- list.files(
        path = file.path(script_dir, "src"),
        pattern = "\\.[Rr]$",
        full.names = TRUE,
        recursive = TRUE
    )
    log_script <- file.path(script_dir, "util", "LogUtil.r")
    GeneUtil <- file.path(script_dir, "util", "GeneUtil.r")
    rfiles <- c(log_script, GeneUtil, rfiles)
    lapply(rfiles, source)
    log_msg("INFO", "load", paste(rfiles, collapse = ", "))
}
