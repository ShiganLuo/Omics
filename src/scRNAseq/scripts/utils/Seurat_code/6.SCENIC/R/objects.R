#' @include utilities.R
NULL

# setOldClass("pins_board_folder")
# setClassUnion("pins_board", c("pins_board_folder", NULL))

#' @title Hephaestus class
#' @slot ... Other slots from
#' \code{\link[Seurat:Seurat]{Seurat}}
#' @importClassesFrom Seurat Seurat
#' @exportClass Hephaestus
#' @keywords internal
setClass(
  Class = "Hephaestus",
  contains = "Seurat",
  slots = c(
    feature = "list",
    board = c("data.frame", NULL),
    reference = "list",
    outdir = c("character")
  )
)

#' @title Construct a Hephaestus object
#' @param assays A 'list' or 'SimpleList' of matrix-like elements
#' All elements of the list must have the same dimensions, we also
#' recommend they have names, e.g. list(Abundance=xx1, RareAbundance=xx2).
#' @param colData An optional DataFrame describing the samples.
#' @param ... additional parameters, see also the usage
#' of \code{\link[Seuart]{Seuart}}.
#' @return Hephaestus object
#' @importFrom methods new
#' @importFrom SeuratObject CreateSeuratObject
#' @export
#' @examples
#' set.seed(717)
#' xx <- matrix(abs(round(rnorm(100, sd = 4), 0)), 10)
#' xx <- data.frame(xx)
#' rownames(xx) <- paste0("row", seq_len(10))
#' hephaestus <- create_hephaestus(assays = xx)
#' hephaestus
create_hephaestus <- function(counts = NULL,
                              metadata = NULL,
                              project = "Hephaestus",
                              board = NULL,
                              outdir = "./",
                              ...) {
  metadata <- as.data.frame(metadata)
  # 1. Create from counts
  if (!is.null(counts) && is.null(metadata)) {
    seurat_obj <-
      CreateSeuratObject(
        counts = counts,
        project = project,
        ...
      )
  } else if (!is.null(metadata) && !is.null(counts)) {
    seurat_obj <-
      CreateSeuratObject(
        counts = counts,
        meta.data = metadata,
        project = "Hephaestus",
        ...
      )
  } else if (is.null(counts) && !is.null(metadata)) {
    # Create from metadata
    # check sample,group,r1 and r2 column name
    requred_colnames <- c("sample", "group", "input_r1")
    if (!all(requred_colnames %in% colnames(metadata))) {
      stop("The colnames of the metadata: sample, group, input_r1, input_r2")
    }
    # check r1, r2 file
    hp_check_file(c(metadata$input_r1, metadata$input_r2))
    counts <- data.frame(row.names = metadata$sample)
    counts <- as.data.frame(t(counts))
    # create board
    x <-
      c(
        "input_r1",
        "input_r2",
        "ip_r1",
        "ip_r2",
        "input_bam",
        "ip_bam"
      )
    x <- x[x %in% colnames(metadata)]
    board <- metadata[, x]
    rownames(board) <- metadata$sample
    metadata <- metadata[, !(colnames(metadata) %in% x)]
    # create
    rownames(metadata) <- metadata$sample
    seurat_obj <-
      suppressWarnings(
        CreateSeuratObject(
          counts = counts,
          meta.data = metadata,
          project = "Hephaestus",
          ...
        )
      )
  }
  hephaestus <- new(
    "Hephaestus",
    seurat_obj,
    board = board,
    outdir = outdir,
    feature = list(),
    reference = list()
  )
  return(hephaestus)
}
