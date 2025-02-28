#' Convert gene identifiers to indices for gene sets - multiple collections
#'
#' @param input tibble with gene_id column and collections list-columns.
#' @param collections character vector of gene collection names; all elements must be column names of input.
#' @param identifiers character vector of gene identifiers.
#' @param gene_id character vector of length 1 representing gene ID; must be a column name of input.
#'
#' @return a list.
#' @export
#'
#' @examples
#' # example
#'
batch_ids2indices <- function(input,
                              collections,
                              identifiers,
                              gene_id = "ensembl_gene_id_version") {

  if(!requireNamespace("limma", quietly = TRUE)) {
    stop("Package \"limma\" must be installed to use this function.")
  }

  if(!tibble::is_tibble(input)) {
    stop("Invalid input argument.")
  }

  if(!is.character(collections) ||
     !all(collections %in% colnames(input)) ||
     !all(unlist(lapply(X = collections,
                        FUN = function(x) { class(input[[x]]) == "list" })))) {
    stop("Invalid collections argument.")
  }

  if(!is.character(gene_id) ||
     length(gene_id) != 1L ||
     !gene_id %in% colnames(input)) {
    stop("Invalid gene_id argument.")
  }

  output <- lapply(X = collections,
                   FUN = function(x) {

                     temp <- tibble::as_tibble(input) %>%
                       dplyr::select(tidyselect::all_of(c(gene_id, x))) %>%
                       tidyr::unnest(tidyselect::all_of(x), keep_empty = F) %>%
                       tidyr::nest(data = -"set_name")

                     temp <- stats::setNames(object = lapply(X = temp$data,
                                                             FUN = function(x) { x[[gene_id]] }),
                                             nm = temp$set_name)

                     temp <- limma::ids2indices(gene.sets = temp,
                                                identifiers = identifiers,
                                                remove.empty = T)

                   })

  stats::setNames(object = output, nm = collections)

}
