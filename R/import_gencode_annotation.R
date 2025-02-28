#' Import GENCODE annotation
#'
#' @param file A GENCODE annotation file in GTF format.
#'
#' @return A tibble with GENCODE annotation.
#' @export
#'
#' @examples
#' ### Do not run ###
#' # gencode_annotation <- import_gencode_annotation(file = "path_to_gencode_annotation.gtf")
#'
import_gencode_annotation <- function(file) {

  if(!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package \"rtracklayer\" must be installed to use this function.")
  }

  if(!grepl(pattern = "^.*\\.gtf$", x = file)) {
    stop("Input file must be in GTF format.")
  }

  output <- rtracklayer::import(file)

  tibble::as_tibble(output)

}


