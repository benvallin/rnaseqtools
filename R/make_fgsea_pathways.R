#' Make fgsea pathways
#'
#' @param input tibble with gene_id column and collections list-columns.
#' @param collections character vector of gene collection names; all elements must be column names of input.
#' @param gene_id character vector of length 1 representing gene ID; must be a column name of input.
#'
#' @return a list of collection-specific sublists containing gene set-specific character vectors of gene IDs.
#' @export
#'
#' @examples
#' ### Do not run ###
#' # collections <- make_fgsea_pathways(input = gene_metadata, collections = collection_names)
#'
make_fgsea_pathways <- function(input,
                                collections,
                                gene_id = "ensembl_gene_id_version") {

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

                     temp <- tibble::as_tibble(input)

                     temp <- temp[, c(gene_id, x)]

                     temp <- tidyr::unnest(data = temp,
                                           cols = tidyselect::all_of(x),
                                           keep_empty = F)

                     temp <- tidyr::nest(.data = temp,
                                         data = tidyselect::all_of(gene_id))

                     temp <- dplyr::mutate(.data = temp,
                                           data = lapply(X = temp$data,
                                                         FUN = function(x) { x[[gene_id]] }))

                     as.list(stats::setNames(object = temp[["data"]],
                                             nm = temp[["set_name"]]))

                   })

  stats::setNames(object = output, nm = collections)

}
