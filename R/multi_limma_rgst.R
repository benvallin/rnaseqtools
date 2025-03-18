#' Run limma rotation gene set tests on multiple collections
#'
#' @param y matrix giving log-expression or log-ratio values with rownames corresponding to gene_id.
#' @param input data.frame or tibble with gene_id column and collections list-columns.
#' @param collections character vector of gene collection names. All elements must be column names of input.
#' @param design design matrix.
#' @param contrast contrast for which the test is required. Can be an integer specifying a column of design, or the name of a column of design, or a numeric contrast vector of length equal to the number of columns of design.
#' @param gene_id character vector of length 1 representing gene ID. Must be a column name of input.
#' @param limma_test character vector of length 1 representing limma rotation gene set test. Must be one of "fry", "roast", or "mroast".
#' @param min_set_size minimal size of a gene set to test. All gene sets below the threshold are excluded.
#' @param max_set_size maximal size of a gene set to test. All gene sets above the threshold are excluded.
#' @param padj_threshold padj threshold. All gene sets with padj equal or above the threshold are filtered out.
#' @param ... optional arguments passed to limma test.
#'
#' @return a tibble with results of limma rotation gene set test runs on multiple collections.
#' @export
#'
#' @examples
#' # Make MSigDB collection table
#' msigdb_collection_table = get_msigdb_collections()
#'
#' # Define collections of interest
#' collection_names <- c("MSigDB_H", "MSigDB_C2_CP:REACTOME", "MSigDB_C2_CP:KEGG", "MSigDB_C5_GO:BP")
#'
#' # Define design matrix
#' model_formula <- "~ treatment + donor_id"
#'
#' design <- model.matrix(object = formula(model_formula), data = sample_metadata_ex)
#'
#' design <- design[, which(colSums(design) != 0), drop = FALSE]
#'
#' # Run limma::fry on selected collections
#' fry_results <- multi_limma_rgst(y = log2_tpm1p_ex,
#'                                 input = msigdb_collection_table,
#'                                 collections = collection_names,
#'                                 design = design,
#'                                 contrast = "treatmenttreated",
#'                                 gene_id = "ensembl_gene_id",
#'                                 limma_test = "fry")
#'
multi_limma_rgst <- function(y,
                             input,
                             collections,
                             design = NULL,
                             contrast = ncol(design),
                             gene_id = "ensembl_gene_id",
                             limma_test = "fry",
                             min_set_size = 1,
                             max_set_size = dim(y)[[1]] - 1,
                             padj_threshold = Inf,
                             ...) {

  if(!requireNamespace("limma", quietly = TRUE)) {
    stop("Package \"limma\" must be installed to use this function.",
         call. = F)
  }

  if(!is.matrix(y) || is.null(rownames(y))) {
    stop("Invalid y argument.",
         call. = F)
  }

  available_tests <- c("fry", "roast", "mroast")

  if(!is.character(limma_test) ||
     length(limma_test) != 1L ||
     !limma_test %in% available_tests) {
    stop("Invalid species argument.",
         call. = F)
  }

  if(!is.numeric(min_set_size) ||
     length(min_set_size) != 1L ||
     min_set_size < 1L) {
    stop("Invalid min_set_size argument.",
         call. = F)
  }

  if(!is.numeric(max_set_size) ||
     length(max_set_size) != 1L ||
     max_set_size <= min_set_size) {
    stop("Invalid max_set_size argument.",
         call. = F)
  }

  if(!is.numeric(padj_threshold) ||
     length(padj_threshold) != 1L ||
     padj_threshold < 0) {
    stop("Invalid padj_threshold argument.",
         call. = F)
  }

  limma_test <- paste0("limma::", limma_test)

  message("Performing ", limma_test, " test...")

  collections <- multi_ids2indices(input = input,
                                   collections = collections,
                                   identifiers = rownames(y),
                                   gene_id = gene_id)

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
