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
