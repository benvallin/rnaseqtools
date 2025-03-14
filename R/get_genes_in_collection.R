#' Get genes in collection
#'
#' @param input data.frame or tibble with gene_ids columns and collection list-column. Collection must contain column set_name.
#' @param collection_name character vector of length 1 representing gene collection name. Must be a column name of input.
#' @param set_name character vector of length 1 representing gene set name in collection.
#' @param gene_ids character vector of gene IDs. All elements must be column names of input.
#'
#' @return a tibble with gene IDs for queried set name in collection.
#' @export
#'
#' @examples
#' # Make MSigDB collection table
#' msigdb_collection_table = get_msigdb_collections()
#'
#' # Define collection of interest
#' collection_name <- "MSigDB_H"
#'
#' # Define gene set name of interest
#' set_name <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
#'
#' # Get genes in gene set of interest
#' genes <- get_genes_in_collection(input = msigdb_collection_table,
#'                                  collection = collection_name,
#'                                  set_name = set_name,
#'                                  gene_ids = "ensembl_gene_id")
#'
get_genes_in_collection <- function(input,
                                    collection_name,
                                    set_name,
                                    gene_ids = c("ensembl_gene_id", "gene_symbol")) {


  if(!is.data.frame(input)) {
    stop("Input must be a data.frame or tibble.",
         call. = F)
  }

  if(!is.character(collection_name) ||
     length(collection_name) != 1L ||
     !collection_name %in% colnames(input) ||
     !is.list(input[[collection_name]])) {
    stop("Invalid collection_name argument.",
         call. = F)
  }

  if(!is.character(set_name) ||
     length(set_name) != 1L) {
    stop("Invalid set_name argument.",
         call. = F)
  }

  if(!is.character(gene_ids) ||
     !all(gene_ids %in% colnames(input))) {
    stop("Invalid gene_ids argument.",
         call. = F)
  }


  output <- input[, c(gene_ids, collection_name)]

  output <- tidyr::unnest(data = output,
                          cols = tidyselect::all_of(collection_name))

  if(!"set_name" %in% colnames(output)) {
    stop("Column set_name not in ", collection_name, " collection.",
         call. = F)
  }

  output <- output[output$set_name == set_name,]

  if(dim(output)[[1]] == 0L) {
    warning(set_name, " is not a gene set of collection ", collection_name, ".\nReturning empty table!\n",
            call. = F, immediate. = T)
  }

  output$collection <- collection_name

  output[, c("collection", "set_name", gene_ids)]

}
