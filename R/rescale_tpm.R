#' Rescale TPMs
#'
#' @param input numeric matrix of TPMs with gene IDs as rownames and sample IDs as colnames.
#' @param genes character vector representing gene IDs to subset for TPM rescaling. Must be of the same gene ID type as rownames(input).
#'
#' @return a matrix with rescaled TPMs for the selected genes.
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
rescale_tpm <- function(input, genes) {

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

  if(!is.character(genes)) {
    stop("genes must be a character vector",
         call. = F)
  }

  genes <- unique(genes)
  n_genes <- length(genes)
  n_genes_in_input <- sum(genes %in% rownames(input), na.rm = TRUE)

  if(n_genes_in_input == 0L) {
    stop("None of the genes are present in rownames of input.\nCheck that <genes> and rownames(input) refer to the same thing!\n",
         call. = F)
  }

  if(n_genes_in_input < n_genes) {
    warning("Only ", n_genes_in_input, " / ", n_genes, " (",
            round(n_genes_in_input / n_genes * 100, 3),
            "%) genes found in rownames of input.\nCheck that <genes> and rownames(input) refer to the same thing!\n",
            call. = F, immediate. = T)
  }

  output <- input[rownames(input) %in% genes, , drop = F]

  output <- sweep(x = output,
                  MARGIN = 2,
                  STATS = colSums(output),
                  FUN = "/") * 100

}

