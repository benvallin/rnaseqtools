#' Make gene metadata
#'
#' @param input A GENCODE annotation tibble produced by import_gencode_annotation.
#'
#' @return A tibble with gene metadata.
#' @export
#'
#' @examples
#' ### Do not run ###
#' # gene_metadata <- make_gene_metadata(input = gencode_annotation)
#'
make_gene_metadata <- function(input) {

  output <- unique(input[!is.na(input$transcript_id),
                         c("gene_id", "gene_name", "gene_type", "seqnames")])

  colnames(output) <- c("ensembl_gene_id_version", "gene_symbol", "gene_type", "chr_name")

  output$ensembl_gene_id <- gsub(pattern = "\\..*$",
                                 replacement = "",
                                 x = output$ensembl_gene_id_version)

  output[, c("ensembl_gene_id_version", "ensembl_gene_id", "gene_symbol", "gene_type", "chr_name")]

}


