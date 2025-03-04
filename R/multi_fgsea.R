#' Run fgsea::fgsea on multiple collections
#'
#' @param input data.frame or tibble with gene_id column and collections list-columns.
#' @param collections character vector of gene collection names. All elements must be column names of input.
#' @param stats named vector of gene-level stats. Names should be the same as in gene_id.
#' @param gene_id character vector of length 1 representing gene ID. Must be a column name of input.
#' @param min_set_size minimal size of a gene set to test. All gene sets below the threshold are excluded.
#' @param max_set_size maximal size of a gene set to test. All gene sets above the threshold are excluded.
#' @param padj_threshold padj threshold. All gene sets with padj equal or above the threshold are filtered out.
#' @param ... optional arguments passed to fgsea::fgsea.
#'
#' @return a tibble with results of fgsea::fgsea runs on multiple collections.
#' @export
#'
#' @examples
#' # Make MSigDB collection table
#' msigdb_collection_table = get_msigdb_collections()
#'
#' # Define collections of interest
#' collection_names <- c("MSigDB_H", "MSigDB_C2_CP:REACTOME", "MSigDB_C2_CP:KEGG", "MSigDB_C5_GO:BP")
#'
#' # Extract gene ranking metric from DESeq2 results
#' stats <- deseq2_results %>%
#'   dplyr::arrange(dplyr::desc(sign_log2fc_times_minus_log10pvalue)) %>%
#'   dplyr::pull(sign_log2fc_times_minus_log10pvalue, ensembl_gene_id)
#'
#' # Run fgsea::fgsea on selected collections
#' fgsea_results <- multi_fgsea(input = msigdb_collection_table,
#'                              collections = collection_names,
#'                              stats = stats,
#'                              gene_id = "ensembl_gene_id")
#'
multi_fgsea <- function(input,
                        collections,
                        stats,
                        gene_id = "ensembl_gene_id",
                        min_set_size = 1,
                        max_set_size = length(stats) - 1,
                        padj_threshold = Inf,
                        ...) {

  if(!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package \"fgsea\" must be installed to use this function.",
         call. = F)
  }

  if(!is.numeric(stats) || is.null(names(stats))) {
    stop("Invalid stats argument.",
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

  stats <- stats::na.omit(stats)

  collections <- make_fgsea_pathways(input = input,
                                     collections = collections,
                                     gene_id = gene_id)

  output <- lapply(X = names(collections),
                   FUN = function(x) {

                     temp <- fgsea::fgsea(pathways = collections[[x]],
                                          stats = stats,
                                          minSize = min_set_size,
                                          maxSize = max_set_size,
                                          ...)

                     temp <- tibble::as_tibble(temp)

                     temp <- dplyr::mutate(.data = temp,
                                           collection = x,
                                           n_leadingEdge = purrr::map_int(.x = .data$leadingEdge,
                                                                          .f = ~ length(.x)),
                                           n_leadingEdge_over_total = paste0(.data$n_leadingEdge, " / ", .data$size),
                                           pct_leadingEdge = (.data$n_leadingEdge / .data$size) * 100)

                     temp <- dplyr::select(.data = temp,
                                           "collection", set_name = "pathway", tidyselect::everything())

                     temp <- dplyr::filter(.data = temp,
                                           .data$padj < padj_threshold)

                     temp <- dplyr::arrange(.data = temp,
                                            dplyr::desc(.data$NES), .data$padj)

                   })

  dplyr::bind_rows(output)

}
