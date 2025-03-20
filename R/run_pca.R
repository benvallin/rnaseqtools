#' Run PCA
#'
#' @param input numeric matrix or data.frame which provides the data for the principal components analysis.
#' @param pcs principal components to return.
#' @param ntop number of top variable features to use for PCA.
#' @param center logical value indicating whether the variables should be shifted to be zero centered.
#' @param scale. logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place.
#' @param rownames character vector of length 1 representing row names.
#'
#' @return a tibble with PCA results.
#' @export
#'
#' @examples
#'
#' # Run PCA
#' pca_results <- run_pca(input = log2_tpm1p_ex,
#'                        pcs = 1:10)
#'
run_pca <- function(input,
                    pcs = 1:2,
                    ntop = NULL,
                    center = T,
                    scale. = T,
                    rownames = "barcode") {

  if(is.null(ntop)) {

    ntop <- dim(input)[[1]]

  }

  var <- matrixStats::rowVars(x = input)

  top_var <- order(var, decreasing = T)[seq_len(min(ntop, length(var)))]

  top_mtx <- t(input[top_var,])

  pca <- stats::prcomp(x = top_mtx,
                       center = center,
                       scale. = scale.)

  pct_var <- (pca$sdev^2/sum(pca$sdev^2))*100

  pct_var <- stats::setNames(object = pct_var[pcs],
                             nm = paste0("PC", pcs))

  pcs <- names(pct_var)

  df <- tibble::as_tibble(x = pca$x[,pcs],
                          rownames = rownames)

  attr(df, "pct_var") <- pct_var

  df

}
