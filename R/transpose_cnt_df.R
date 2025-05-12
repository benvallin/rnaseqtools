#' Transpose count data.frame
#'
#' @param input data.frame or tibble with character column <gene_id> or <sample_id> and double columns representing sample-specific gene counts.
#' @param gene_id character vector of length 1 or name representing gene ID.
#' @param sample_id character vector of length 1 or name representing sample ID.
#'
#' @return a transposed tibble of gene counts.
#' @export
#'
#' @examples
#' # Convert log2(TPM+1) matrix to tibble with gene IDs as column
#' log2_tpm1p_df <- tibble::as_tibble(x = sc_log2_tpm1p_ex,
#'                                    rownames = "ensembl_gene_id")
#'
#' # Transpose log2(TPM+1) tibble with barcode as column
#' log2_tpm1p_df_t <- transpose_cnt_df(input = log2_tpm1p_df,
#'                                     gene_id = "ensembl_gene_id",
#'                                     sample_id = "barcode")
#'
transpose_cnt_df <- function(input,
                             gene_id = "ensembl_gene_id",
                             sample_id = "barcode") {

  gene_id <- rlang::enexpr(gene_id)
  gene_id <- ifelse(is.symbol(gene_id), deparse(gene_id), eval(gene_id))

  sample_id <- rlang::enexpr(sample_id)
  sample_id <- ifelse(is.symbol(sample_id), deparse(sample_id), eval(sample_id))

  if(!is.character(gene_id) ||
     length(gene_id) != 1L) {
    stop("Invalid gene_id argument.",
         call. = F)
  }

  if(!is.character(sample_id) ||
     length(sample_id) != 1L) {
    stop("Invalid sample_id argument.",
         call. = F)
  }

  if(!is.data.frame(input) ||
     (!gene_id %in% colnames(input) && !sample_id %in% colnames(input)) ||
     (gene_id %in% colnames(input) && sample_id %in% colnames(input))) {
    stop("Input must be a data.frame or tibble with character column <gene_id> or <sample_id> and double columns representing sample-specific gene counts.",
         call. = F)
  }

  if(gene_id %in% colnames(input)) {

    input_col <- gene_id
    output_col <- sample_id

  } else {

    input_col <- sample_id
    output_col <- gene_id

  }

  if(!is.character(input[[input_col]]) ||
     !all(purrr::map_lgl(.x = setdiff(colnames(input), input_col),
                         .f = ~ is.double(input[[.x]])))) {
    stop("Invalid column type in input.",
         call. = F)
  }

  output <- tibble::column_to_rownames(.data = input,
                                       var = input_col)

  output <- t(output)

  tibble::as_tibble(x = tibble::rownames_to_column(.data = as.data.frame(output),
                                                   var = output_col))

}
