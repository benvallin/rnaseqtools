#' Search genes in leading edges
#'
#' @param input data.frame or tibble with leadingEdge list-column and n_leadingEdge column, as produced by rnaseqtools::multi_fgsea.
#' @param gene_ids character vector of gene IDs to search in leading edges. Must be the same type of gene IDs as in the leadingEdge list-column.
#' @param group_name character vector of length 1 representing group name for the searched genes.
#'
#' @return the input table augmented with gene search results.
#' @export
#'
#' @examples
#' # Get leading edge for KEGG_RIBOSOME gene set
#' genes <- get_leading_edge(input = multi_fgsea_results,
#'                           set_name = "KEGG_RIBOSOME")
#'
#' # Search genes representing KEGG_RIBOSOME leading edge in all leading edges
#' search_results <- search_genes_in_leading_edges(input = multi_fgsea_results,
#'                                                 gene_ids = genes,
#'                                                 group_name = "ribosome")
#'
search_genes_in_leading_edges <- function(input,
                                          gene_ids,
                                          group_name = "searched") {

  if(!is.data.frame(input) ||
     !all(c("leadingEdge", "n_leadingEdge") %in% colnames(input)) ||
     !is.list(input[["leadingEdge"]])) {
    stop("Input must be a data.frame or tibble with leadingEdge list-column and n_leadingEdge column.",
         call. = F)
  }

  if(!is.character(gene_ids)) {
    stop("Invalid gene_ids argument.",
         call. = F)
  }

  if(!is.character(group_name) ||
     length(group_name) != 1L) {
    stop("Invalid group_name argument.",
         call. = F)
  }

  output <- input %>%
    dplyr::mutate(leadingEdge_searched = purrr::map(.x = .data$leadingEdge,
                                                    .f = ~ intersect(.x, gene_ids)),
                  n_leadingEdge_searched = purrr::map_int(.x = .data$leadingEdge_searched,
                                                          .f = ~ length(.x)),
                  n_leadingEdge_searched_over_leadingEdge = paste0(.data$n_leadingEdge_searched, " / ", .data$n_leadingEdge),
                  pct_leadingEdge_searched_in_leadingEdge = (.data$n_leadingEdge_searched / .data$n_leadingEdge) * 100) %>%
    dplyr::arrange(dplyr::desc(.data$pct_leadingEdge_searched_in_leadingEdge),
                   .data$collection) %>%
    dplyr::mutate(set_number = seq_along(input[[1]])) %>%
    dplyr::rename_with(.fn = function(x) gsub(pattern = "searched", replacement = group_name, x = x),
                       .cols = c("leadingEdge_searched",
                                 "n_leadingEdge_searched",
                                 "n_leadingEdge_searched_over_leadingEdge",
                                 "pct_leadingEdge_searched_in_leadingEdge")) %>%
    dplyr::select("set_number", tidyselect::everything())

  if(all(purrr::map_lgl(.x = output[[paste0("n_leadingEdge_", group_name)]],
                        .f = ~ .x == 0L))) {

    warning("None of the gene IDs were found in the leading edges!\n",
            "Check that gene_ids and input's leadingEdge are the same type of gene ID.\n",
            call. = F, immediate. = T)

  }

  output

}
