#' Make transcript metadata
#'
#' @param input A GENCODE annotation data.frame or tibble as produced by import_gencode_annotation.
#'
#' @return A tibble with transcript metadata.
#' @export
#'
#' @examples
#' # Download GENCODE annotation file
#' persistent_data_dir_path <- tools::R_user_dir(package = "rnaseqtools", which = "data")
#'
#' download_gencode_annotation(species = "human",
#'                             release = 46,
#'                             out_dir_path = persistent_data_dir_path)
#'
#' # Import GENCODE annotation from downloaded file
#' gencode_file_path <- paste0(persistent_data_dir_path,
#'                             "/",
#'                             "gencode.v46.primary_assembly.annotation.gtf.gz")
#'
#' gencode_annotation <- import_gencode_annotation(file = gencode_file_path)
#'
#' # Make transcript metadata from GENCODE annotation
#' transcript_metadata <- make_transcript_metadata(input = gencode_annotation)
#'
#' # Clean up
#' file.remove(gencode_file_path)
#'
#' rm(persistent_data_dir_path, gencode_file_path, gencode_annotation)
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


