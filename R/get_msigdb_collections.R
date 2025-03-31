#' Get MSigDB collections
#'
#' @param species character vector of length 1 representing species for output genes. Must be a species_name returned by msigdbr::msigdbr_species.
#' @param db_species character vector of length 1 representing species abbreviation for database query. Must be "Hs" or "HS" for human and "Mm" or "MM" for mouse databases.
#'
#' @return a tibble with ensembl_gene_id column and collections list-columns. Each collection contains columns set_id and set_name.
#' @export
#'
#' @examples
#' msigdb_collection_table = get_msigdb_collections()
#'
get_msigdb_collections <- function(species = "Homo sapiens",
                                   db_species = "Hs") {

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

  output <- msigdbr::msigdbr(species = species,
                             db_species = db_species,
                             collection = NULL,
                             subcollection = NULL)

  output$collection_id <- ifelse(output$gs_subcollection == "",
                                 paste0("MSigDB_", output$gs_collection),
                                 paste0("MSigDB_", output$gs_collection, "_", output$gs_subcollection))

  output <- output[, c("collection_id",  "gs_id", "gs_name", "db_ensembl_gene")]

  colnames(output) <- c("collection_id", "set_id", "set_name", "ensembl_gene_id")

  output <- tidyr::nest(.data = output,
                        data = tidyselect::all_of(c("set_id", "set_name")))

  output <- tidyr::pivot_wider(data = output,
                               names_from = "collection_id",
                               values_from = "data")

  output[, sort(colnames(output))]

}
