#' Convert gene identifiers to indices for gene sets - multiple collections
#'
#' @param input tibble with gene_id column and collections list-columns.
#' @param collections character vector of gene collection names. All elements must be column names of input.
#' @param identifiers character vector of gene identifiers.
#' @param gene_id character vector of length 1 representing gene ID. Must be a column name of input.
#'
#' @return a list of collection-specific sublists containing gene set-specific integer vectors of gene indices in the vector identifiers.
#' @export
#'
#' @examples
#' ### Do not run ###
#' # collections <- batch_ids2indices(input = msigdb_collection_table,
#' #                                  collections = collection_names,
#' #                                  identifiers = msigdb_collection_table$ensembl_gene_id,
#' #                                  gene_id = "ensembl_gene_id")
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

  if(!is.character(identifiers)) {
    stop("Invalid identifiers argument.")
  }

  if(!is.character(gene_id) ||
     length(gene_id) != 1L ||
     !gene_id %in% colnames(input)) {
    stop("Invalid gene_id argument.")
  }

  identifiers <- unique(identifiers)

  n_identifiers <- length(identifiers)
  n_genes <- length(unique(input[[gene_id]]))
  n_identifiers_in_genes <- sum(identifiers %in% unique(input[[gene_id]]))

  if(n_identifiers_in_genes < n_identifiers) {

    warning("Only ", n_identifiers_in_genes, " / ", n_identifiers,
            " (", round(n_identifiers_in_genes / n_genes * 100, 3), "%) identifiers found in column ", gene_id, " of input.",
            "\nCheck that identifiers and gene_id refer to the same thing!",
            call. = F, immediate. = T)

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
