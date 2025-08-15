#' Search genes in overlaps
#'
#' @param input data.frame or tibble with columns collection and overlap and list-column overlapGenes, as produced by rnaseqtools::multi_fora.
#' @param genes character vector of gene IDs to search in overlaps. Must be the same type of gene IDs as in the overlapGenes list-column.
#' @param group_name character vector of length 1 representing group name for the searched genes.
#'
#' @return the input table augmented with gene search results.
#' @export
#'
#' @examples
#' # Get overlap for KEGG_RIBOSOME gene set
#' genes <- get_overlap(input = multi_fora_results_ex,
#'                      set_name = "KEGG_RIBOSOME")
#'
#' # Search genes representing KEGG_RIBOSOME overlap in all overlaps
#' search_results <- search_genes_in_overlaps(input = multi_fora_results_ex,
#'                                            genes = genes,
#'                                            group_name = "ribosome")
#'
search_genes_in_overlaps <- function(input,
                                     genes,
                                     group_name = "searched") {

  if(!is.data.frame(input) ||
     !all(c("collection", "overlapGenes", "overlap") %in% colnames(input)) ||
     !is.list(input[["overlapGenes"]])) {
    stop("Input must be a data.frame or tibble with columns collection and overlap and list-column overlapGenes.",
         call. = F)
  }

  if(!is.character(genes)) {
    stop("Invalid genes argument.",
         call. = F)
  }

  if(!is.character(group_name) ||
     length(group_name) != 1L) {
    stop("Invalid group_name argument.",
         call. = F)
  }

  output <- input %>%
    dplyr::mutate(overlapGenes_searched = purrr::map(.x = .data$overlapGenes,
                                                     .f = ~ intersect(.x, genes)),
                  n_overlapGenes_searched = purrr::map_int(.x = .data$overlapGenes_searched,
                                                           .f = ~ length(.x)),
                  n_overlapGenes_searched_over_overlapGenes = paste0(.data$n_overlapGenes_searched, " / ", .data$overlap),
                  pct_overlapGenes_searched_in_overlapGenes = (.data$n_overlapGenes_searched / .data$overlap) * 100) %>%
    dplyr::arrange(dplyr::desc(.data$pct_overlapGenes_searched_in_overlapGenes),
                   .data$collection) %>%
    dplyr::mutate(set_number = seq_along(input[[1]])) %>%
    dplyr::rename_with(.fn = function(x) gsub(pattern = "searched", replacement = group_name, x = x),
                       .cols = c("overlapGenes_searched",
                                 "n_overlapGenes_searched",
                                 "n_overlapGenes_searched_over_overlapGenes",
                                 "pct_overlapGenes_searched_in_overlapGenes")) %>%
    dplyr::select("set_number", tidyselect::everything())

  if(all(purrr::map_lgl(.x = output[[paste0("n_overlapGenes_", group_name)]],
                        .f = ~ .x == 0L))) {

    warning("None of the gene IDs were found in the overlaps!\n",
            "Check that genes and input's overlapGenes are the same type of gene ID.\n",
            call. = F, immediate. = T)

  }

  output

}
