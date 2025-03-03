#' Make transcript metadata
#'
#' @param input A GENCODE annotation data.frame or tibble as produced by import_gencode_annotation.
#'
#' @return A tibble with transcript metadata.
#' @export
#'
#' @examples
#' ### Do not run ###
#' # transcript_metadata <- make_transcript_metadata(input = gencode_annotation)
#'
make_transcript_metadata <- function(input) {

  keep <- c("transcript_id", "gene_id", "gene_name")

  if(!is.data.frame(input) || !all(keep %in% colnames(input))) {
    stop("Input must be a tibble with columns \"transcript_id\", \"gene_id\" and \"gene_name\".",
         call. = F)
  }

  output <- unique(input[!is.na(input$transcript_id), keep])

  colnames(output) <- c("ensembl_transcript_id_version", "ensembl_gene_id_version", "gene_symbol")

  output

}


