#' Make transcript metadata
#'
#' @param input A GENCODE annotation tibble produced by import_gencode_annotation.
#'
#' @return A tibble with transcript metadata.
#' @export
#'
#' @examples
#' ### Do not run ###
#' # transcript_metadata <- make_transcript_metadata(input = gencode_annotation)
#'
make_transcript_metadata <- function(input) {

  output <- unique(input[!is.na(input$transcript_id), c("transcript_id", "gene_id")])

  colnames(output) <- c("ensembl_transcript_id_version", "ensembl_gene_id_version")

  output

}


