#' Run limma gene set tests on multiple collections
#'
#' @param y matrix giving log-expression or log-ratio values with rownames corresponding to gene_id.
#' @param input data.frame or tibble with gene_id column and collections list-columns.
#' @param collections character vector of gene collection names. All elements must be column names of input.
#' @param design design matrix.
#' @param contrast contrast for which the test is required. Can be an integer specifying a column of design, or the name of a column of design, or a numeric contrast vector of length equal to the number of columns of design.
#' @param gene_id character vector of length 1 representing gene ID. Must be a column name of input.
#' @param min_set_size minimal size of a gene set to test. All gene sets below the threshold are excluded.
#' @param max_set_size maximal size of a gene set to test. All gene sets above the threshold are excluded.
#' @param padj_threshold padj threshold. All gene sets with padj equal or above the threshold are filtered out.
#' @param limma_test character vector of length 1 representing limma gene set test.
#' @param ... optional arguments passed to limma test.
#'
#' @return a tibble with results of limma test runs on multiple collections.
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
multi_limma_gsea <- function(y,
                             input,
                             collections,
                             design,
                             contrast,
                             gene_id = "ensembl_gene_id",
                             min_set_size = 1,
                             max_set_size = dim(y)[[1]] - 1,
                             padj_threshold = Inf,
                             limma_test = "fry",
                             ...) {

  if(!requireNamespace("limma", quietly = TRUE)) {
    stop("Package \"limma\" must be installed to use this function.",
         call. = F)
  }

  limma_test <- paste0("limma::", limma_test)

  message("Performing ", limma_test, " test...")

  collections <- multi_ids2indices(input = input,
                                   collections = collections,
                                   identifiers = rownames(y))

  output <- lapply(X = names(collections),
                   FUN = function(x) {

                     sets <- collections[[x]]
                     sets <- sets[purrr::map_lgl(.x = sets, .f = ~ length(.x) >= min_set_size)]
                     sets <- sets[purrr::map_lgl(.x = sets, .f = ~ length(.x) <= max_set_size)]

                     if(length(sets) == 0) {

                       NULL

                     } else {

                       temp <- do.call(what = eval(parse(text = limma_test)),
                                        args = list(y = y,
                                                    index = sets,
                                                    design = design,
                                                    contrast = contrast,
                                                    ...))

                       temp <- tibble::rownames_to_column(.data = temp,
                                                          var = "set_name") %>%
                         tibble::as_tibble() %>%
                         dplyr::mutate(collection = x) %>%
                         dplyr::select("collection", "set_name", tidyselect::everything())

                       if("FDR" %in% colnames(temp)) {

                         temp <- temp %>%
                           dplyr::filter(.data$FDR < padj_threshold) %>%
                           dplyr::arrange(dplyr::desc(.data$Direction), .data$FDR)

                       }
                     }
                   })

  dplyr::bind_rows(output)

}
