#' Import GENCODE annotation
#'
#' @param file A GENCODE annotation file in GTF format.
#'
#' @return A tibble with GENCODE annotation.
#' @export
#'
#' @examples
#' # Download GENCODE annotation file
#' out_dir_path <- tools::R_user_dir(package = "rnaseqtools", which = "data")
#'
#' download_gencode_annotation(species = "human",
#'                             release = 46,
#'                             out_dir_path = out_dir_path)
#'
#' # Import GENCODE annotation from downloaded file
#' gencode_file_path <- paste0(out_dir_path, "/gencode.v46.primary_assembly.annotation.gtf.gz")
#'
#' gencode_annotation <- import_gencode_annotation(file = gencode_file_path)
#'
#' # Clean up
#' file.remove(gencode_file_path)
#'
#' rm(out_dir_path, gencode_file_path)
#'
import_gencode_annotation <- function(file) {

  if(!requireNamespace("rtracklayer", quietly = TRUE)) {

    stop("Package \"rtracklayer\" must be installed to use this function.",
         call. = F)

  }

  if(!grepl(pattern = "^.*\\.gtf$", x = file) &&
     !grepl(pattern = "^.*\\.gtf\\.gz$", x = file)) {

    stop("Input file must be in GTF format.",
         call. = F)

  }

  output <- rtracklayer::import(file)

  tibble::as_tibble(output)

}


