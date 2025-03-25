#' Get MSigDB collections
#'
#' @param species character vector of length 1 representing species for output genes. Must be a species_name returned by msigdbr::msigdbr_species.
#' @param db_species character vector of length 1 representing species abbreviation for database query. Must be "Hs" or "HS" for human and "Mm" or "MM" for mouse databases.
#' @param msigdb_collection_summary data.frame or tibble with columns collection_id, gs_collection and gs_subcollection as returned by rnaseqtools::list_msigdb_collections.
#'
#' @return a tibble with ensembl_gene_id column and collections list-columns. Each collection contains columns set_id and set_name.
#' @export
#'
#' @examples
#' msigdb_collection_table = get_msigdb_collections()
#'
get_msigdb_collections <- function(species = "Homo sapiens",
                                   db_species = "Hs",
                                   msigdb_collection_summary = list_msigdb_collections(db_species = "Hs")) {

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

  if(!is.character(db_species) ||
     length(db_species) != 1L ||
     !db_species %in% c("HS", "Hs", "MM", "Mm")) {
    stop("Invalid db_species argument.",
         call. = F)
  }

  if(!is.data.frame(msigdb_collection_summary) ||
     !all(c("collection_id", "gs_collection", "gs_subcollection") %in% colnames(msigdb_collection_summary))) {
    stop("Invalid msigdb_collection_summary argument.",
         call. = F)
  }

  collections <- apply(X = msigdb_collection_summary,
                       MARGIN = 1,
                       FUN = function(x) { as.list(x) },
                       simplify = F)

  output <- lapply(X = collections,
                   FUN = function(x) {

                     collection_id <- x$collection_id

                     collection <- x$gs_collection

                     if(x$gs_subcollection == "") {
                       subcollection <- NULL
                     } else {
                       subcollection <- x$gs_subcollection
                     }

                     temp <- msigdbr::msigdbr(species = species,
                                              collection = collection,
                                              subcollection = subcollection)

                     temp <- temp[, c("ensembl_gene", "gs_id", "gs_name")]

                     colnames(temp) <- c("ensembl_gene_id", "set_id", "set_name")

                     temp <- stats::na.omit(temp)

                     temp <- tidyr::nest(.data = temp,
                                         .by = "ensembl_gene_id",
                                         .key = collection_id)

                   })

  Reduce(f = function(...) dplyr::full_join(..., by = dplyr::join_by("ensembl_gene_id")),
         x = output)

}
