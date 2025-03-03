#' Make fgsea pathways
#'
#' @param input data.frame or tibble with gene_id column and collections list-columns. Each collection must contain column set_name.
#' @param collections character vector of gene collection names. All elements must be column names of input.
#' @param gene_id character vector of length 1 representing gene ID. Must be a column name of input.
#'
#' @return a list of collection-specific sublists containing gene set-specific character vectors of gene IDs.
#' @export
#'
#' @examples
#' ### Do not run ###
#' # collections <- make_fgsea_pathways(input = msigdb_collection_table,
#' #                                    collections = collection_names,
#' #                                    gene_id = "ensembl_gene_id")
#'
make_fgsea_pathways <- function(input,
                                collections,
                                gene_id = "ensembl_gene_id_version") {

  if(!is.data.frame(input)) {
    stop("Input must be a data.frame or tibble.",
         call. = F)
  }

  if(!is.character(collections) ||
     !all(collections %in% colnames(input)) ||
     !all(unlist(lapply(X = collections,
                        FUN = function(x) { class(input[[x]]) == "list" })))) {
    stop("Invalid collections argument.",
         call. = F)
  }

  if(!is.character(gene_id) ||
     length(gene_id) != 1L ||
     !gene_id %in% colnames(input)) {
    stop("Invalid gene_id argument.",
         call. = F)
  }

  output <- lapply(X = collections,
                   FUN = function(x) {

                     temp <- tibble::as_tibble(input)

                     temp <- temp[, c(gene_id, x)]

                     temp <- tidyr::unnest(data = temp,
                                           cols = tidyselect::all_of(x),
                                           keep_empty = F)

                     if(!"set_name" %in% colnames(temp)) {
                       stop("Column set_name not in ", x, " collection.",
                            call. = F)
                     }

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
