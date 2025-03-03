#' Get MSigBD collections
#'
#' @param species character vector of length 1 representing species. Must be a species_name returned by msigdbr::msigdbr_species.
#' @param msigdb_collection_summary data.frame or tibble with columns collection_name, gs_cat and gs_subcat as returned by list_msigdb_collections.
#'
#' @return a tibble with ensembl_gene_id column and collections list-columns. Each collection contains columns set_id and set_name.
#' @export
#'
#' @examples
#' msigdb_collection_table = get_msigdb_collections()
#'
get_msigdb_collections <- function(species = "Homo sapiens",
                                   msigdb_collection_summary = list_msigdb_collections()) {

  if(!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package \"msigdbr\" must be installed to use this function.",
         call. = F)
  }

  available_species <- msigdbr::msigdbr_species()
  available_species <- available_species$species_name

  if(!is.character(species) ||
     length(species) != 1L ||
     !species %in% available_species) {
    stop("Invalid species argument.",
         call. = F)
  }

  if(!is.data.frame(msigdb_collection_summary) ||
     !all(c("collection_name", "gs_cat", "gs_subcat") %in% colnames(msigdb_collection_summary))) {
    stop("Invalid msigdb_collection_summary argument.",
         call. = F)
  }

  collections <- apply(X = msigdb_collection_summary,
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
