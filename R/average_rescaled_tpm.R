#' Average rescaled TPMs
#'
#' @param input numeric matrix of rescaled TPMs with gene IDs as rownames and sample IDs as colnames. Typically produced by rnaseqtools::rescale_tpm.
#' @param grp_df data.frame or tibble with column sample_id and one or more grouping columns.
#' @param sample_id character vector of length 1 representing sample ID. Must be a column name of grp_df. Sample IDs should be of the same type in colnames of input and column sample_id of group_df.
#' @param gene_id character vector of length 1 representing gene ID. Used to label row names of input.
#'
#' @return a matrix with rescaled TPMs averaged by groups.
#' @export
#'
#' @examples
#' # Define genes of interest
#' genes <- calcium_genes_ex %>%
#'   dplyr::filter(protein_complex == "SERCA") %>%
#'   dplyr::pull(ensembl_gene_id)
#'
#' # Rescale TPMs for selected genes
#' resc_tpm <- rescale_tpm(input = bulk_tpm_ex, genes = genes)
#'
#' # Average rescaled TPMs by disease status using a grp_df providing sample metadata
#' mean_resc_tpm <- average_rescaled_tpm(input = resc_tpm,
#'                                       grp_df = bulk_sample_metadata_ex,
#'                                       sample_id = donor_id,
#'                                       gene_id = ensembl_gene_id)
#'
average_rescaled_tpm <- function(input,
                                 grp_df,
                                 sample_id = "sample_name",
                                 gene_id = "ensembl_gene_id") {

  sample_id <- rlang::enexpr(sample_id)
  sample_id <- ifelse(is.symbol(sample_id), deparse(sample_id), eval(sample_id))

  gene_id <- rlang::enexpr(gene_id)
  gene_id <- ifelse(is.symbol(gene_id), deparse(gene_id), eval(gene_id))

  if(!is.matrix(input) ||
     !is.numeric(input) ||
     is.null(rownames(input)) ||
     is.null(colnames(input))) {
    stop("input must be a numerical matrix with gene IDs as rownames and sample IDs as colnames.",
         call. = F)
  }

  if(!all(round(colSums(input)) == 100L)) {
    stop("input is not a matrix of rescaled TPMs: its column sums are not all equal to 100.",
         call. = F)
  }

  if(!is.character(sample_id) ||
     length(sample_id) != 1L) {
    stop("Invalid gene_id argument.",
         call. = F)
  }

  if(!is.data.frame(grp_df) ||
     length(colnames(grp_df)) < 2L ||
     !sample_id %in% colnames(grp_df)) {
    stop("grp_df must be a data.frame / tibble containing at least 2 columns: <sample_id> and one or more grouping columns.",
         call. = F)
  }

  sample_ids <- intersect(colnames(input),
                          grp_df[[sample_id]])

  if(length(sample_ids) == 0L) {
    stop("None of the sample IDs in column <sample_id> of group_df are present in colnames of input.",
         call. = F)
  }

  grp_vars <- setdiff(colnames(grp_df), sample_id)

  if(!all(purrr::map_lgl(.x = grp_vars,
                         .f = ~ is.character(grp_df[[.x]]) || is.factor(grp_df[[.x]])))) {
    stop("All grouping variables in grp_df should be of class character or factor.",
         call. = F)
  }

  if(!is.character(gene_id) ||
     length(gene_id) != 1L ||
     is.na(gene_id)) {
    stop("Invalid gene_id argument.",
         call. = F)
  }

  message("Grouping by ", paste0(grp_vars, collapse = ", "), " and ", gene_id, " for averaging rescaled TPMs...\n")

  tibble::as_tibble(input,
                    rownames = gene_id) %>%
    dplyr::select(gene_id, sample_ids) %>%
    tidyr::pivot_longer(cols = -gene_id,
                        names_to = sample_id,
                        values_to = "rescaled_tpm") %>%
    dplyr::left_join(grp_df,
                     by = sample_id) %>%
    dplyr::group_by(dplyr::across(.cols = c(grp_vars, gene_id))) %>%
    dplyr::mutate(mean_rescaled_tpm = mean(.data$rescaled_tpm, na.rm = T),
                  sd_rescaled_tpm = stats::sd(.data$rescaled_tpm, na.rm = T),
                  sem_rescaled_tpm = .data$sd_rescaled_tpm / length(.data[[sample_id]])) %>%
    dplyr::ungroup() %>%
    dplyr::select(colnames(grp_df),
                  tidyselect::everything())

}
