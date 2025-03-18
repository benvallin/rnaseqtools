#' Tidy DESeq2 results
#'
#' @param input DESeqResults object produced by DESeq2::results.
#' @param gene_id character vector of length 1 representing gene ID. Used to create a column of gene IDs in output.
#'
#' @return a tibble with DESeq2 results.
#' @export
#'
#' @examples
#' # Extract DESeq2 log fold change results from DESeqResults object
#' deseq2_results <- tidy_DESeqResults(input = DESeqResults_ex,
#'                                     gene_id = "ensembl_gene_id")
#'
tidy_DESeqResults <- function(input,
                              gene_id = "ensembl_gene_id") {

  if(!requireNamespace("DESeq2", quietly = TRUE)) {

    stop("Package \"DESeq2\" must be installed to use this function.",
         call. = F)

  }

  if(!methods::is(input, "DESeqResults") ||
     !all(c("log2FoldChange", "pvalue") %in% names(input))) {
    stop("Input must be a DESeqResults object produced by DESeq2::results.",
         call. = F)
  }

  if(!is.character(gene_id) ||
     length(gene_id) != 1L) {
    stop("Invalid gene_id argument.",
         call. = F)
  }

  output <- tibble::as_tibble(input, rownames = gene_id)

  output[["sign_log2fc_times_minus_log10pvalue"]] = sign(output[["log2FoldChange"]]) * -log(x = output[["pvalue"]],
                                                                                            base = 10)

  output %>%
    dplyr::arrange(dplyr::desc(.data$sign_log2fc_times_minus_log10pvalue))

}
