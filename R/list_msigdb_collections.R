#' List MSigDB collections
#'
#' @return a tibble summarizing available MSigDB collections with columns collection_name, gs_cat and gs_subcat.
#' @export
#'
#' @examples
#' msigdb_collection_summary = list_msigdb_collections()
#'
list_msigdb_collections <- function() {

  if(!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package \"msigdbr\" must be installed to use this function.",
         call. = F)
  }

  output <- msigdbr::msigdbr_collections()

  output <- output[, c("gs_cat", "gs_subcat")]

  output$gs_subcat <- ifelse(output$gs_subcat == "",
                             NA_character_,
                             output$gs_subcat)

  output$collection_name <- ifelse(is.na(output$gs_subcat),
                                   paste0("MSigDB_", output$gs_cat),
                                   paste0("MSigDB_", output$gs_cat, "_", output$gs_subcat))

  output[, c("collection_name", "gs_cat", "gs_subcat")]

}
