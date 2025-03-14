#' Get leading edge
#'
#' @param input data.frame or tibble with set_name column and leadingEdge list-column, as produced by rnaseqtools::multi_fgsea.
#' @param set_name character vector of length 1 representing gene set name.
#'
#' @return a character vector of gene IDs representing the leading edge for the queried gene set.
#' @export
#'
#' @examples
#' # Get genes in leading edge
#' leading_edge <- get_leading_edge(input = multi_fgsea_results,
#'                                  set_name = "REACTOME_PHOSPHOLIPID_METABOLISM")
#'
get_leading_edge <- function(input,
                             set_name) {

  if(!is.data.frame(input) ||
     !all(c("set_name", "leadingEdge") %in% colnames(input)) ||
     !is.list(input[["leadingEdge"]])) {
    stop("Input must be a data.frame or tibble with set_name column and leadingEdge list-column.",
         call. = F)
  }

  if(!is.character(set_name) ||
     length(set_name) != 1L) {
    stop("Invalid set_name argument.",
         call. = F)
  }

  if(!set_name %in% input$set_name) {
    stop("Gene set ", set_name, " is not in input.",
         call. = F)
  }

  input[input$set_name == set_name, "leadingEdge"][[1]][[1]]

}
