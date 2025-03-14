#' Get overlapping genes
#'
#' @param input data.frame or tibble with set_name column and overlapGenes list-column, as produced by rnaseqtools::multi_fora.
#' @param set_name character vector of length 1 representing gene set name.
#'
#' @return a character vector of gene IDs representing the overlapping genes for the queried gene set.
#' @export
#'
#' @examples
#' # Get genes in overlap
#' overlap <- get_overlap(input = multi_fora_results,
#'                        set_name = "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
#'
get_overlap <- function(input,
                        set_name) {

  if(!is.data.frame(input) ||
     !all(c("set_name", "overlapGenes") %in% colnames(input)) ||
     !is.list(input[["overlapGenes"]])) {
    stop("Input must be a data.frame or tibble with set_name column and overlapGenes list-column.",
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

  input[input$set_name == set_name, "overlapGenes"][[1]][[1]]

}
