#' MsaAAMultipleAlignment to tibble
#'
#' @param input MsaAAMultipleAlignment object produced by msa::msa.
#'
#' @return a tibble with columns name and sequence.
#' @export
#'
#' @examples
#' # Extract multiple sequence alignment from MsaAAMultipleAlignment object
#' msa_df <- msa_to_tibble(input = msa_ex)
#'
msa_to_tibble <- function(input) {

  if(!requireNamespace("Biostrings", quietly = TRUE)) {

    stop("Package \"Biostrings\" must be installed to use this function.",
         call. = F)

  }

  if(!methods::is(input, "MsaAAMultipleAlignment") ||
     !"unmasked" %in% methods::slotNames(input)) {
    stop("Input must be a MsaAAMultipleAlignment object produced by msa::msa.",
         call. = F)
  }

  temp <- purrr::map_chr(.x = Biostrings::unmasked(input),
                         .f = as.character)

  tibble::tibble(name = names(temp),
                 sequence = unname(temp))

}


