#' Convert gene identifiers to indices for gene sets - multiple collections
#'
#' @param input data.frame or tibble with gene_id column and collections list-columns.
#' @param collections character vector of gene collection names. All elements must be column names of input.
#' @param identifiers character vector of gene identifiers.
#' @param gene_id character vector of length 1 representing gene ID. Must be a column name of input.
#'
#' @return a list of collection-specific sublists containing gene set-specific integer vectors of gene indices in the vector identifiers.
#' @export
#'
#' @examples
#' # Make MSigDB collection table
#' msigdb_collection_table = get_msigdb_collections()
#'
#' # Define collections of interest
#' collection_names <- c("MSigDB_H", "MSigDB_C2_CP:REACTOME", "MSigDB_C2_CP:KEGG", "MSigDB_C5_GO:BP")
#'
#' # Extract gene identifiers from log2 COUNT+1 matrix
#' identifiers <- rownames(log2_tpm1p)
#'
#' # Get lists of indices for selected collections
#' collections <- multi_ids2indices(input = msigdb_collection_table,
#'                                  collections = collection_names,
#'                                  identifiers = identifiers,
#'                                  gene_id = "ensembl_gene_id")
#'
multi_ids2indices <- function(input,
                              collections,
                              identifiers,
                              gene_id = "ensembl_gene_id") {

  if(!requireNamespace("limma", quietly = TRUE)) {
    stop("Package \"limma\" must be installed to use this function.",
         call. = F)
  }

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

  if(!is.character(identifiers)) {
    stop("Invalid identifiers argument.",
         call. = F)
  }

  if(!is.character(gene_id) ||
     length(gene_id) != 1L ||
     !gene_id %in% colnames(input)) {
    stop("Invalid gene_id argument.",
         call. = F)
  }

  identifiers <- unique(identifiers)

  n_identifiers <- length(identifiers)
  n_identifiers_in_input <- sum(identifiers %in% unique(input[[gene_id]]))

  if(n_identifiers_in_input < n_identifiers) {

    warning("Only ", n_identifiers_in_input, " / ", n_identifiers, " (",
            round(n_identifiers_in_input / n_identifiers * 100, 3), "%) identifiers found in column ",
            gene_id, " of input.\nCheck that identifiers and gene_id refer to the same thing!\n",
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
