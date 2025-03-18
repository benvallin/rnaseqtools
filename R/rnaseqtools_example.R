#' Get path to rnaseqtools example
#'
#' rnaseqtools comes bundled with a number of sample files in its inst/extdata directory. This function make them easy to access.
#'
#' @param file name of file. If NULL, the example files will be listed.
#'
#' @export
#'
#' @examples
#' rnaseqtools_example()
#'
#' rnaseqtools_example(file = "meta_info_ex.json")
#'
rnaseqtools_example <- function(file = NULL) {

  if (is.null(file)) {

    dir(system.file("extdata", package = "rnaseqtools"))

  } else {

    system.file("extdata", file, package = "rnaseqtools", mustWork = T)

  }

}
