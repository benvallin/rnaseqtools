#' List MSigDB collections
#'
#' @param db_species character vector of length 1 representing species abbreviation for database query. Must be "Hs" or "HS" for human and "Mm" or "MM" for mouse databases.
#'
#' @return a tibble summarizing available MSigDB collections.
#' @export
#'
#' @examples
#' msigdb_collection_summary = list_msigdb_collections()
#'
list_msigdb_collections <- function(db_species = "Hs") {

  if(!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package \"msigdbr\" must be installed to use this function.",
         call. = F)
  }

  if(!requireNamespace("msigdbdf", quietly = TRUE)) {
    stop("Package \"msigdbdf\" must be installed to use this function.",
         call. = F)
  }

  if(!is.character(db_species) ||
     length(db_species) != 1L ||
     !db_species %in% c("HS", "Hs", "MM", "Mm")) {
    stop("Invalid db_species argument.",
         call. = F)
  }

  output <- msigdbr::msigdbr_collections(db_species = db_species)

  output$collection_id <- ifelse(output$gs_subcollection == "",
                                 paste0("MSigDB_", output$gs_collection),
                                 paste0("MSigDB_", output$gs_collection, "_", output$gs_subcollection))

  output[, c("collection_id", setdiff(colnames(output), "collection_id"))]

}
