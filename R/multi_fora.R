#' Run fgsea::fora on multiple collections
#'
#' @param input data.frame or tibble with gene_id column and collections list-columns.
#' @param collections character vector of gene collection names. All elements must be column names of input.
#' @param genes set of query genes.
#' @param universe a universe from which genes were selected.
#' @param gene_id character vector of length 1 representing gene ID. Must be a column name of input.
#' @param min_set_size minimal size of a gene set to test. All gene sets below the threshold are excluded.
#' @param max_set_size maximal size of a gene set to test. All gene sets above the threshold are excluded.
#' @param padj_threshold padj threshold. All gene sets with padj equal or above the threshold are filtered out.
#' @param ... optional arguments passed to fgsea::fora.
#'
#' @return a tibble with results of fgsea::fora runs on multiple collections.
#' @export
#'
#' @examples
#' ### Do not run ###
#' # fora_results <- multi_ora(input = msigdb_collection_table,
#' #                           collections = collection_names,
#' #                           genes = genes,
#' #                           universe = universe,
#' #                           gene_id = "ensembl_gene_id")
#'
multi_fora <- function(input,
                       collections,
                       genes,
                       universe,
                       gene_id = "ensembl_gene_id_version",
                       min_set_size = 1,
                       max_set_size = NULL,
                       padj_threshold = Inf,
                       ...) {

  if(!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package \"fgsea\" must be installed to use this function.",
         call. = F)
  }

  if(!is.character(genes)) {
    stop("Invalid genes argument.",
         call. = F)
  }

  if(!is.character(universe)) {
    stop("Invalid universe argument.",
         call. = F)
  }

  if(!is.numeric(min_set_size) ||
     length(min_set_size) != 1L ||
     min_set_size < 1L) {
    stop("Invalid min_set_size argument.",
         call. = F)
  }

  if(!is.null(max_set_size)) {

    if(!is.numeric(max_set_size) ||
       length(max_set_size) != 1L ||
       max_set_size <= min_set_size) {
      stop("Invalid max_set_size argument.",
           call. = F)
    }

  }

  if(!is.numeric(padj_threshold) ||
     length(padj_threshold) != 1L ||
     padj_threshold < 0) {
    stop("Invalid padj_threshold argument.",
         call. = F)
  }

  collections <- make_fgsea_pathways(input = input,
                                     collections = collections,
                                     gene_id = gene_id)

  output <- lapply(X = names(collections),
                   FUN = function(x) {

                     temp_collection <- input %>%
                       dplyr::select(tidyselect::all_of(c(gene_id, x))) %>%
                       tidyr::unnest(tidyselect::all_of(x), keep_empty = F)

                     temp_universe <- temp_collection %>%
                       dplyr::filter(.data[[gene_id]] %in% universe) %>%
                       dplyr::pull(.data[[gene_id]]) %>%
                       unique()

                     temp_genes <- temp_collection %>%
                       dplyr::filter(.data[[gene_id]] %in% genes) %>%
                       dplyr::pull(.data[[gene_id]]) %>%
                       unique()

                     temp_max_set_size <- ifelse(is.null(max_set_size),
                                                 length(temp_universe),
                                                 max_set_size)

                     temp <- fgsea::fora(pathways = collections[[x]],
                                         genes = temp_genes,
                                         universe = temp_universe,
                                         minSize = min_set_size,
                                         maxSize = temp_max_set_size,
                                         ...)

                     temp <- tibble::as_tibble(temp) %>%
                       dplyr::mutate(collection = x,
                                     n_overlapGenes = paste0(.data$overlap, " / ", .data$size),
                                     pct_overlapGenes = (.data$overlap / .data$size) * 100) %>%
                       dplyr::select("collection", set_name = "pathway", tidyselect::everything()) %>%
                       dplyr::filter(.data$padj < padj_threshold) %>%
                       dplyr::arrange(.data$padj)

                   })

  dplyr::bind_rows(output)

}
