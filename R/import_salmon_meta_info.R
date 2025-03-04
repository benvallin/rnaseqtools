#' Import Salmon's metadata
#'
#' @param file A meta_info.json file produced by Salmon software.
#'
#' @return A tibble with key metrics extracted from input meta_info.json file.
#' @export
#'
#' @examples
#' # Import Salmon's run metadata
#' meta_info <- import_salmon_meta_info(file = rnaseqtools_example(file = "meta_info.json"))
#'
import_salmon_meta_info <- function(file) {

  if(!requireNamespace("rjson", quietly = TRUE)) {
    stop("Package \"rjson\" must be installed to use this function.",
         call. = F)
  }

  if(!grepl(pattern = "^.*\\.json$", x = file)) {
    stop("Input file must be in JSON format.",
         call. = F)
  }

  output <- rjson::fromJSON(file = file)

  keep <- c("library_types",
            "frag_dist_length", "frag_length_mean", "frag_length_sd",
            "seq_bias_correct", "gc_bias_correct", "keep_duplicates",
            "num_valid_targets", "num_decoy_targets",
            "num_processed", "num_mapped",
            "num_decoy_fragments", "num_dovetail_fragments",
            "num_fragments_filtered_vm", "num_alignments_below_threshold_for_mapped_fragments_vm",
            "percent_mapped")

  if(!all(keep %in% names(output))) {
    stop("Missing Salmon's metrics in input file.",
         call. = F)
  }

  output <- tibble::as_tibble(output[keep])

  output$file <- file

  output[, c("file", keep)]

}
