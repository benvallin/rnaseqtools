#' Split txi by feature type
#'
#' @param txi txi object produced by tximport::tximport.
#' @param split_df data.frame or tibble with 2 columns representing feature IDs. Each column represents a separate feature type, and the features in a given row are considered representatives of the same feature.
#'
#' @return a list with distinct abundance, counts and length vectors for the 2 feature types.
#' @export
#'
#' @examples
#'
#' # Split txi by splicing status (spliced / unspliced)
#' txi_split <- split_splici_txi(txi = txi_ex,
#'                               split_df = split_df_ex)
#'
split_splici_txi <- function(txi, split_df) {

  if(!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package \"SummarizedExperiment\" must be installed to use this function.",
         call. = F)
  }

  if(!requireNamespace("tximeta", quietly = TRUE)) {
    stop("Package \"tximeta\" must be installed to use this function.",
         call. = F)
  }

  if(!is.list(txi) ||
     length(txi) != 4L ||
     !all(c("abundance", "counts", "length", "countsFromAbundance") %in% names(txi))) {
    stop("Invalid txi argument.",
         call. = F)
  }

  if(!is.data.frame(split_df) ||
     dim(split_df)[[2]] != 2L) {
    stop("Invalid split_df argument.",
         call. = F)
  }

  output <- lapply(X = c("abundance", "counts", "length"),
                   FUN = function(layer) {

                     se <- SummarizedExperiment::SummarizedExperiment(assays = stats::setNames(object = list(txi[[layer]]),
                                                                                               nm = layer))

                     txis <- tximeta::splitSE(se = se,
                                              splitDf = split_df,
                                              assayName = layer)

                     stats::setNames(object = list(SummarizedExperiment::assay(txis,
                                                                               colnames(split_df)[[1]]),
                                                   SummarizedExperiment::assay(txis,
                                                                               colnames(split_df)[[2]])),
                                     nm = c(colnames(split_df)[[1]],
                                            colnames(split_df)[[2]]))

                   })

  output <- stats::setNames(object = c(output, txi$countsFromAbundance),
                            nm = c("abundance", "counts", "length", "countsFromAbundance"))

}
