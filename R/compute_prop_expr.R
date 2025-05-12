#' Compute proportion of expressing cells
#'
#' @param input data.frame or tibble with character column <gene_id> and double columns representing sample-specific gene counts.
#' @param gene_id character vector of length 1 or name representing gene ID. Must be a column name of input.
#' @param sample_metadata data.frame or tibble with character / factor columns <sample_id_var> and <group_id_vars>.
#' @param sample_id_var character vector of length 1 or name representing sample ID. Must be a column name of sample_metadata.
#' @param group_id_vars character vector or vector of names representing grouping variables. Must be column name(s) of sample_metadata.
#' @param min_cnt numeric vector of length 1 representing count threshold above which a gene should be considered as expressed.
#'
#' @return a tibble with number of cells, number of expressing cells and proportion of expressing cells per biological group ID.
#' @export
#'
#' @examples
#' # Convert log2(TPM+1) matrix to TPM tibble with gene IDs as column
#' tpm_df <- tibble::as_tibble(x = (2^sc_log2_tpm1p_ex)-1,
#'                             rownames = "ensembl_gene_id")
#'
#' # Compute proportion of expressing cells per disease status
#' propr_expr <- compute_prop_expr(input = tpm_df,
#'                                 gene_id = ensembl_gene_id,
#'                                 sample_metadata = sc_sample_metadata_ex,
#'                                 sample_id_var = barcode,
#'                                 group_id_vars = treatment,
#'                                 min_cnt = 10)
#'
compute_prop_expr <- function(input,
                              gene_id = "ensembl_gene_id",
                              sample_metadata,
                              sample_id_var,
                              group_id_vars,
                              min_cnt = 0) {

  gene_id <- rlang::enexpr(gene_id)
  gene_id <- ifelse(is.symbol(gene_id), deparse(gene_id), eval(gene_id))

  sample_id_var <- rlang::enexpr(sample_id_var)
  sample_id_var <- ifelse(is.symbol(sample_id_var), deparse(sample_id_var), eval(sample_id_var))

  group_id_vars <- rlang::enexpr(group_id_vars)

  if(is.call(group_id_vars)) {

    for(i in seq_along(group_id_vars)[2:length(group_id_vars)]) {
      group_id_vars[[i]] <- ifelse(is.symbol(group_id_vars[[i]]),
                                   deparse(group_id_vars[[i]]),
                                   eval(group_id_vars[[i]]))
    }

    group_id_vars <- eval(group_id_vars)

  } else {

    group_id_vars <- ifelse(is.symbol(group_id_vars), deparse(group_id_vars), eval(group_id_vars))

  }

  if(!is.character(gene_id) ||
     length(gene_id) != 1L) {
    stop("Invalid gene_id argument.",
         call. = F)
  }

  if(!is.character(sample_id_var) ||
     length(sample_id_var) != 1L) {
    stop("Invalid sample_id_var argument.",
         call. = F)
  }

  if(!is.character(group_id_vars)) {
    stop("Invalid group_id_vars argument.",
         call. = F)
  }

  if(!is.data.frame(input) ||
     !gene_id %in% colnames(input) ||
     !is.character(input[[gene_id]]) ||
     !all(purrr::map_lgl(.x = setdiff(colnames(input), gene_id),
                         .f = ~ is.double(input[[.x]])))) {
    stop("Input must be a data.frame or tibble with only character column <gene_id> and double columns representing sample-specific gene counts.",
         call. = F)
  }

  if(!is.data.frame(sample_metadata) ||
     !all(c(sample_id_var, group_id_vars) %in% colnames(sample_metadata)) ||
     (!is.character(sample_metadata[[sample_id_var]]) && !is.factor(sample_metadata[[sample_id_var]])) ||
     any(purrr::map_lgl(.x = group_id_vars,
                        .f = ~ (!is.character(sample_metadata[[.x]]) && !is.factor(sample_metadata[[.x]]))))) {

    stop("sample_metadata must be a data.frame or tibble with character / factor columns <sample_id_var> and <group_id_vars>.",
         call. = F)
  }

  if(!all(setdiff(colnames(input), gene_id) %in% unique(as.character(sample_metadata[[sample_id_var]])))) {
    stop("The names of double columns in input must all be present in column <sample_id_var> of sample_metadata.",
         call. = F)
  }

  if(!is.numeric(min_cnt) ||
     length(min_cnt) != 1L) {
    stop("Invalid min_cnt argument.",
         call. = F)
  }

  output <- tibble::column_to_rownames(.data = input,
                                       var = gene_id)

  output <- output > min_cnt

  output <- tibble::rownames_to_column(.data = as.data.frame(output),
                                       var = gene_id)

  output <- tidyr::pivot_longer(data = output,
                                cols = setdiff(colnames(output), gene_id),
                                names_to = sample_id_var,
                                values_to = "cnt")

  output <- dplyr::left_join(x = output,
                             y = sample_metadata[, c(sample_id_var, group_id_vars)],
                             by = sample_id_var)

  output <- tidyr::nest(.data = output,
                        sample_data = tidyselect::all_of(c(sample_id_var, "cnt")))

  output$n <- purrr::map_int(.x = output$sample_data,
                             .f = ~ dim(.x)[[1]])

  output$n_expressing <- purrr::map_int(.x = output$sample_data,
                                        .f = ~ sum(.x$cnt))

  output$pct_expressing <- output$n_expressing / output$n * 100

  output[, c(gene_id, group_id_vars, "n", "n_expressing", "pct_expressing")]

}
