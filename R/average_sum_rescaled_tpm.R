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

  if(!is.data.frame(gs_df) ||
     length(colnames(gs_df)) != 2L ||
     !gene_id %in% colnames(gs_df)) {
    stop("gs_df must be a data.frame containing only 2 columns: <gene_id> and a column representing gene sets.",
         call. = F)
  }

  gs_var <- setdiff(colnames(gs_df), gene_id)

  if(gs_var %in% colnames(grp_df)) {
    stop("Column ", gs_var, " is present both in gs_df and grp_df.",
         call. = F)
  }

  if(!is.character(gs_df[[gs_var]]) &&
     !is.factor(gs_df[[gs_var]])) {
    stop("The grouping variable ", gs_var, " in gs_df should be of class character or factor.",
         call. = F)
  }

  if(!all(rownames(input) %in% gs_df[[gene_id]])) {
    stop("The gene IDs in rownames of input are not all present in the column <gene_id> of gs_df.",
         call. = F)
  }

  if(!is.character(gene_id) ||
     length(gene_id) != 1L ||
     is.na(gene_id)) {
    stop("Invalid gene_id argument.",
         call. = F)
  }

  # message("Grouping samples by ", paste0(grp_vars, collapse = ", "), "...")

  tibble::as_tibble(input,
                    rownames = gene_id) %>%
    dplyr::select(gene_id, sample_ids) %>%
    tidyr::pivot_longer(cols = -gene_id,
                        names_to = sample_id,
                        values_to = "rescaled_tpm") %>%
    dplyr::left_join(grp_df,
                     by = sample_id) %>%
    dplyr::left_join(gs_df,
                     by = gene_id) %>%
    dplyr::filter(!is.na(gs_var)) %>%
    dplyr::group_by(dplyr::across(.cols = c(grp_vars, sample_id, gs_var))) %>%
    dplyr::mutate(sum_rescaled_tpm = sum(.data$rescaled_tpm, na.rm = T)) %>%
    dplyr::group_by(across(.cols = c(grp_vars, gs_var))) %>%
    dplyr::mutate(mean_sum_rescaled_tpm = mean(.data$sum_rescaled_tpm, na.rm = T),
                  sd_sum_rescaled_tpm = stats::sd(.data$sum_rescaled_tpm, na.rm = T),
                  sem_sum_rescaled_tpm = .data$sd_sum_rescaled_tpm / length(.data[[sample_id]])) %>%
    dplyr::ungroup() %>%
    dplyr::select(colnames(grp_df),
                  tidyselect::everything(),
                  -c(gene_id, rescaled_tpm)) %>%
    unique()

}
