#' Average summed rescaled TPMs
#'
#' @param input numeric matrix of rescaled TPMs with gene IDs as rownames and sample IDs as colnames. Typically produced by rnaseqtools::rescale_tpm.
#' @param grp_df data.frame or tibble with column <sample_id> and one or more grouping columns.
#' @param gs_df data.frame or tibble with two columns: <gene_id> and a column representing gene sets.
#' @param sample_id character vector of length 1 or name representing sample ID. Must be a column name of grp_df. Sample IDs should be of the same type in colnames of input and column <sample_id> of grp_df.
#' @param gene_id character vector of length 1 or name representing gene ID. Gene IDs should be of the same type in rownames of input and column <gene_id> of gs_df.
#'
#' @return a matrix with rescaled TPMs summed by sample ID and gene set and then averaged by groups.
#' @export
#'
#' @examples
#' # Define gs_df
#' gs_df <- calcium_genes_ex %>%
#'   dplyr::select(ensembl_gene_id, protein_complex) %>%
#'   dplyr::filter(!is.na(protein_complex))
#'
#' # Pull genes of interest from gs_df
#' genes <- gs_df %>%
#'   dplyr::pull(ensembl_gene_id)
#'
#' # Rescale TPMs for selected genes
#' resc_tpm <- rescale_tpm(input = bulk_tpm_ex, genes = genes)
#'
#' # Sum rescaled TPMs by donor ID and protein complex using a gs_df providing gene metadata
#' # and then average summed rescaled TPMs by disease status using a grp_df providing sample metadata
#' mean_sum_resc_tpm <- average_sum_rescaled_tpm(input = resc_tpm,
#'                                               grp_df = bulk_sample_metadata_ex,
#'                                               gs_df = gs_df,
#'                                               sample_id = donor_id,
#'                                               gene_id = ensembl_gene_id)
#'
average_sum_rescaled_tpm <- function(input,
                                     grp_df,
                                     gs_df,
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

  if(length(unique(rownames(input))) != length(rownames(input))) {
    stop("Gene IDs in rownames of input are not all unique.",
         call. = F)
  }

  if(length(unique(colnames(input))) != length(colnames(input))) {
    stop("Sample IDs in colnames of input are not all unique.",
         call. = F)
  }

  if(!all(round(colSums(input)) == 100L)) {
    stop("input is not a matrix of rescaled TPMs: its column sums are not all equal to 100.",
         call. = F)
  }

  if(!is.character(sample_id) ||
     length(sample_id) != 1L) {
    stop("Invalid sample_id argument.",
         call. = F)
  }

  if(!is.data.frame(grp_df) ||
     length(colnames(grp_df)) < 2L ||
     !sample_id %in% colnames(grp_df)) {
    stop("grp_df must be a data.frame / tibble containing at least 2 columns: <sample_id> and one or more grouping columns.",
         call. = F)
  }

  if(length(unique(grp_df[[sample_id]])) != length(grp_df[[sample_id]])) {
    stop("Sample IDs in column <sample_id> of grp_df are not all unique.",
         call. = F)
  }

  sample_ids <- intersect(colnames(input),
                          grp_df[[sample_id]])

  if(length(sample_ids) == 0L) {
    stop("None of the sample IDs in column <sample_id> of grp_df are present in colnames of input.",
         call. = F)
  }

  if(!is.character(gene_id) ||
     length(gene_id) != 1L ||
     is.na(gene_id)) {
    stop("Invalid gene_id argument.",
         call. = F)
  }

  if(!is.data.frame(gs_df) ||
     length(colnames(gs_df)) != 2L ||
     !gene_id %in% colnames(gs_df)) {
    stop("gs_df must be a data.frame containing only 2 columns: <gene_id> and a column representing gene sets.",
         call. = F)
  }

  if(length(unique(gs_df[[gene_id]])) != length(gs_df[[gene_id]])) {
    stop("Gene IDs in column <gene_id> of gs_df are not all unique.",
         call. = F)
  }

  gene_ids <- intersect(rownames(input),
                        gs_df[[gene_id]])

  if(length(gene_ids) == 0L) {
    stop("None of the gene IDs in column <gene_id> of gs_df are present in rownames of input.",
         call. = F)
  }

  grp_vars <- setdiff(colnames(grp_df), sample_id)

  if(!all(purrr::map_lgl(.x = grp_vars,
                         .f = ~ is.character(grp_df[[.x]]) || is.factor(grp_df[[.x]])))) {
    stop("All grouping variables in grp_df should be of class character or factor.",
         call. = F)
  }

  gs_var <- setdiff(colnames(gs_df), gene_id)

  if(!is.character(gs_df[[gs_var]]) &&
     !is.factor(gs_df[[gs_var]])) {
    stop("The grouping variable ", gs_var, " in gs_df should be of class character or factor.",
         call. = F)
  }

  if(gs_var %in% colnames(grp_df)) {
    stop("Column ", gs_var, " is present both in gs_df and grp_df.",
         call. = F)
  }

  if(!all(colnames(input) %in% grp_df[[sample_id]])) {
    warning("The sample IDs in colnames of input are not all present in the column <sample_id> of grp_df. Those samples will be ommitted.\n",
            call. = F, immediate. = T)
  }

  if(!all(grp_df[[sample_id]] %in% colnames(input))) {
    warning("The sample IDs in the column <sample_id> of grp_df are not all present in colnames of input. Those samples will be ommitted.\n",
            call. = F, immediate. = T)
  }

  if(!all(rownames(input) %in% gs_df[[gene_id]])) {
    warning("The gene IDs in rownames of input are not all present in the column <gene_id> of gs_df. Those genes will be omitted.\n",
            call. = F, immediate. = T)
  }

  if(!all(gs_df[[gene_id]] %in% rownames(input))) {
    warning("The gene IDs in the column <gene_id> of gs_df are not all present in rownames of input. Those genes will be omitted.\n",
            call. = F, immediate. = T)
  }

  if(any(is.na(gs_df[[gs_var]]))) {
    warning("Missing ", gs_var, " value for one or more gene IDs in gs_df. Those genes will be omitted.\n",
            call. = F, immediate. = T)
  }

  message("Grouping by ", sample_id, " and ", gs_var, " for summing rescaled TPMs...\n")

  message("Grouping by ", paste0(grp_vars, collapse = ", "), " and ", gs_var, " for averaging summed rescaled TPMs...\n")

  grp_var_nm <- paste0(grp_vars, collapse = "_")

  grp_df[[grp_var_nm]] <- Reduce(f = function(...) paste(..., sep = "_"),
                                 x = grp_df[, grp_vars])

  grp_df <- grp_df[, c(sample_id, grp_var_nm)]

  tibble::as_tibble(input,
                    rownames = gene_id) %>%
    dplyr::filter(.data[[gene_id]] %in% gene_ids) %>%
    dplyr::select(gene_id, sample_ids) %>%
    tidyr::pivot_longer(cols = -gene_id,
                        names_to = sample_id,
                        values_to = "rescaled_tpm") %>%
    dplyr::left_join(grp_df,
                     by = sample_id) %>%
    dplyr::left_join(gs_df,
                     by = gene_id) %>%
    dplyr::group_by(dplyr::across(.cols = c(sample_id, gs_var))) %>%
    dplyr::mutate(sum_rescaled_tpm = sum(.data$rescaled_tpm, na.rm = T)) %>%
    dplyr::group_by(dplyr::across(.cols = c(grp_var_nm, gs_var))) %>%
    dplyr::mutate(mean_sum_rescaled_tpm = mean(.data$sum_rescaled_tpm, na.rm = T),
                  sd_sum_rescaled_tpm = stats::sd(.data$sum_rescaled_tpm, na.rm = T),
                  sem_sum_rescaled_tpm = .data$sd_sum_rescaled_tpm / length(.data[[sample_id]])) %>%
    dplyr::ungroup() %>%
    tidyr::nest(data = c(gene_id, "rescaled_tpm")) %>%
    dplyr::select(colnames(grp_df),
                  tidyselect::everything())

}
