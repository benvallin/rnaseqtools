#' List MSigBD collections
#'
#' @return a tibble with columns collection_name, gs_cat and gs_subcat.
#' @export
#'
#' @examples
#' ### Do not run ###
#' # collection_table = list_msigdb_collections()
#'
list_msigdb_collections <- function() {

  if(!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package \"msigdbr\" must be installed to use this function.")
  }

  output <- msigdbr::msigdbr_collections()

  output <- output[, c("gs_cat", "gs_subcat")]

  output$gs_subcat <- ifelse(output$gs_subcat == "",
                             NA_character_,
                             output$gs_subcat)

  output$collection_name <- paste0("MSigDB_",
                                   output$gs_cat, "_", output$gs_subcat)

  output[, c("collection_name", "gs_cat", "gs_subcat")]

}
