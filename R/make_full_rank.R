#' Make matrix full rank
#'
#' @param model_matrix model matrix.
#'
#' @return the input matrix if already full rank.
#' A subset matrix if the input matrix is not full rank due to column(s) full of zeros.
#' An error if the input matrix is not full rank due to linearly dependent columns.
#' @export
#'
#' @examples
#'
#' # Add genotype information to sample metadata
#' sample_metadata <- sc_sample_metadata_ex
#'
#' sample_metadata$genotype <- ifelse(sample_metadata$donor_id == "donor1", "CTRL", "PD")
#'
#' # Define model formula
#' model_formula <- "~ donor_id + treatment + genotype"
#'
#' ### Case when matrix is not full rank due to columns full of zeros
#'
#' # Define model matrix (excluding donor2)
#' model_mtx <- model.matrix(object = formula(model_formula),
#'                           data = sample_metadata[sample_metadata$donor_id != "donor2",])
#'
#' # Attempt to make model matrix full rank
#' message("Error captured:\n",
#'         attr(x = try(make_full_rank(model_matrix = model_mtx), silent = TRUE),
#'              which = "condition")$message)
#'
#' ### Case when matrix is not full rank due to linearly dependent columns
#'
#' # Define model matrix (using full sample metadata)
#' model_mtx <- model.matrix(object = formula(model_formula),
#'                           data = sample_metadata)
#'
#' # Attempt to make model matrix full rank
#' message("Error captured:\n",
#'         attr(x = try(make_full_rank(model_matrix = model_mtx), silent = TRUE),
#'              which = "condition")$message)
#'
make_full_rank <- function(model_matrix) {

  if(!is.matrix(model_matrix) || !is.numeric(model_matrix)) {

    stop("Input model_matrix is not a numeric matrix.",
         call. = F)

  }

  # if(!identical(sort(unique(as.vector(model_matrix))), c(0, 1))) {
  #
  #   stop("Input model_matrix should only contain zeros and ones.",
  #        call. = F)
  #
  # }

  is_full_rank <- function(mtx) { !(qr(mtx)$rank < ncol(mtx)) }

  if (!is_full_rank(mtx = model_matrix)) {

    full_zero_column <- any(colSums(model_matrix) == 0)

    if(full_zero_column) {

      message(
        "Model matrix is not full rank.",
        "\nCause: levels or combinations of levels without any samples have resulted in column(s) full of zeros in the model matrix.",
        "\nPossible fix: discarding full zero column(s) and checking subset model matrix again..."
      )

      output <- model_matrix[, which(colSums(model_matrix) != 0), drop = F]

      if(!is_full_rank(mtx = output)) {

        stop(
          "\nModel matrix still not full rank after discarding full zero column(s).",
          "\nCause: one or more variables or interaction terms in the design formula are linear combinations of the others and must be removed.",
          call. = F
        )

      } else {

        message("Worked! Returning subset model matrix.")

      }

    } else {

      stop(
        "Model matrix is not full rank.",
        "\nCause: one or more variables or interaction terms in the design formula are linear combinations of the others and must be removed.",
        call. = F
      )

    }

  } else {

    message("Model matrix is already full rank.")

    output <- model_matrix

  }

  output

}
