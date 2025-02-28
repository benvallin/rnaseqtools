#' Get MSigBD collections
#'
#' @param species character vector of length 1 representing species; must be a species_name returned by msigdbr::msigdbr_species.
#' @param collection_table tibble with columns collection_name, gs_cat and gs_subcat as returned by list_msigdb_collections.
#'
#' @return a tibble with ensembl_gene_id column and collections list-columns.
#' @export
#'
#' @examples
#' ### Do not run ###
#' # collections = get_msigdb_collections()
#'
get_msigdb_collections <- function(species = "Homo sapiens",
                                   collection_table = list_msigdb_collections()) {

  if(!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package \"msigdbr\" must be installed to use this function.")
  }

  available_species <- msigdbr::msigdbr_species()
  available_species <- available_species$species_name

  if(!is.character(species) ||
     length(species) != 1L ||
     !species %in% available_species) {
    stop("Invalid species argument.")
  }

  if(!tibble::is_tibble(collection_table) ||
     !all(c("collection_name", "gs_cat", "gs_subcat") %in% colnames(collection_table))) {
    stop("Invalid collection_table argument.")
  }

  collections <- apply(X = collection_table,
                       MARGIN = 1,
                       FUN = function(x) { as.list(x) },
                       simplify = F)

  output <- lapply(X = collections,
                   FUN = function(x) {

                     collection_name <- x$collection_name

                     category <- x$gs_cat

                     if(is.na(x$gs_subcat)) {
                       subcategory <- NULL
                     } else {
                       subcategory <- x$gs_subcat
                     }

                     temp <- msigdbr::msigdbr(species = species,
                                              category = category,
                                              subcategory = subcategory)

                     temp <- temp[, c("ensembl_gene", "gs_id", "gs_name")]

                     colnames(temp) <- c("ensembl_gene_id", "set_id", "set_name")

                     temp <- stats::na.omit(temp)

                     temp <- tidyr::nest(.data = temp,
                                         .by = "ensembl_gene_id",
                                         .key = collection_name)

                   })

  Reduce(f = function(...) dplyr::full_join(..., by = dplyr::join_by("ensembl_gene_id")),
         x = output)

}
