#' Run fgsea::fora on multiple collections
#'
#' @param input data.frame or tibble with gene_id column and collections list-columns.
#' @param collections character vector of gene collection names. All elements must be column names of input.
#' @param genes set of query genes. Should be in the gene_id column of input.
#' @param universe a universe from which genes were selected. Should be in the gene_id column of input.
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
#' # Make MSigDB collection table
#' msigdb_collection_table = get_msigdb_collections()
#'
#' # Define collections of interest
#' collection_names <- c("MSigDB_H", "MSigDB_C2_CP:REACTOME", "MSigDB_C2_CP:KEGG", "MSigDB_C5_GO:BP")
#'
#' # Extract genes significantly downregulated (log2fc < 0 and padj < 0.05) from MAST results
#' genes <- mast_results %>%
#'   dplyr::filter(log2fc < 0,
#'                 padj < 0.05) %>%
#'   dplyr::pull(ensembl_gene_id)
#'
#' # Extract corresponding universe from MAST results
#' universe <- mast_results %>%
#'   dplyr::filter(!is.na(log2fc)) %>%
#'   dplyr::pull(ensembl_gene_id)
#'
#' # Run fgsea::fora on selected collections
#' fora_results <- multi_fora(input = msigdb_collection_table,
#'                            collections = collection_names,
#'                            genes = genes,
#'                            universe = universe,
#'                            gene_id = "ensembl_gene_id")
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

  n_genes <- length(unique(genes))
  n_genes_in_input <- sum(genes %in% unique(input[[gene_id]]))

  if(n_genes_in_input < n_genes) {

    warning("Only ", n_genes_in_input, " / ", n_genes, " (",
            round(n_genes_in_input / n_genes * 100, 3), "%) genes found in column ",
            gene_id, " of input.\nCheck that genes and gene_id refer to the same thing!\n",
            call. = F, immediate. = T)

  }

  n_universe <- length(unique(universe))
  n_universe_in_input <- sum(universe %in% unique(input[[gene_id]]))

  if(n_universe_in_input < n_universe) {

    warning("Only ", n_universe_in_input, " / ", n_universe, " (",
            round(n_universe_in_input / n_universe * 100, 3), "%) universe genes found in column ",
            gene_id, " of input.\nCheck that universe and gene_id refer to the same thing!\n",
            call. = F, immediate. = T)

  }

  n_genes_not_in_universe <- length(unique(genes[!genes %in% universe]))

  if(n_genes_not_in_universe > 0L) {

    warning(n_genes_not_in_universe, " / ", n_genes,
            " genes are not in universe.\nMake sure that all genes belong to the universe!\n",
            call. = F, immediate. = T)

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
