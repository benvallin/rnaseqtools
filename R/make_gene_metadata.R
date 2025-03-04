#' Make gene metadata
#'
#' @param input A GENCODE annotation data.frame or tibble as produced by import_gencode_annotation.
#'
#' @return A tibble with gene metadata.
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
#' # Make gene metadata from GENCODE annotation
#' gene_metadata <- make_gene_metadata(input = gencode_annotation)
#'
#' # Clean up
#' file.remove(gencode_file_path)
#'
#' rm(persistent_data_dir_path, gencode_file_path, gencode_annotation)
#'
make_gene_metadata <- function(input) {

  keep <- c("transcript_id", "gene_id", "gene_name", "gene_type", "seqnames")

  if(!is.data.frame(input) || !all(keep %in% colnames(input))) {
    stop("Input must be a tibble with columns \"transcript_id\", \"gene_id\", \"gene_name\", \"gene_type\" and \"seqnames\".",
         call. = F)
  }

  keep <- setdiff(keep, "transcript_id")

  output <- unique(input[!is.na(input$transcript_id), keep])

  colnames(output) <- c("ensembl_gene_id_version", "gene_symbol", "gene_type", "chr_name")

  output$ensembl_gene_id <- gsub(pattern = "\\..*$",
                                 replacement = "",
                                 x = output$ensembl_gene_id_version)

  output[, c("ensembl_gene_id_version", "ensembl_gene_id", "gene_symbol", "gene_type", "chr_name")]

}


