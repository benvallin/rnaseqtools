#' Plot PCA
#'
#' @param input data.frame or tibble with columns containing PC values as produced by rnaseqtools::run_pca, and optionally columns containing fill and shape categorical variables.
#' @param pcs numeric vector of length 2 representing the principal components to plot on x and y axes.
#' @param fill character vector of length 1 representing the categorical variable to use for color-coding datapoints. Should be a column name of input.
#' @param fill_lab character vector of length 1 representing the fill label to use in plot legend.
#' @param fill_val character vector of colors to use for plotting.
#' @param shape character vector of length 1 representing the categorical variable to use for shape-coding datapoints. Should be a column name of input.
#' @param shape_lab character vector of length 1 representing the shape label to use in plot legend.
#' @param shape_val numeric vector of shapes to use for plotting.
#'
#' @return a ggplot object with selected PCs on x and y axes, and optionally color- and shape-coding of datapoints by categorical variables.
#' @export
#'
#' @examples
#' # Run PCA
#' pca_results <- run_pca(input = log2_tpm1p_ex,
#'                        pcs = 1:10)
#'
#' # Join sample metadata to PCA results
#' pca_results <- pca_results %>%
#'   dplyr::left_join(sample_metadata_ex,
#'                    by = dplyr::join_by(barcode))
#'
#' # Plot PC1 against PC2
#' plot_pca(input = pca_results,
#'          pcs = c(1,2),
#'          fill = "donor_id",
#'          fill_lab = NULL,
#'          fill_val = c("#1f77b4","#e377c2", "#279e68"),
#'          shape = "treatment",
#'          shape_lab = NULL,
#'          shape_val = c(21, 23))
#'
plot_pca <- function(input,
                     pcs = c(1,2),
                     fill = NULL,
                     fill_lab = NULL,
                     fill_val = c("#1f77b4", "#ff7f0e", "#279e68", "#8c564b","#e377c2","#aec7e8","#b5bd61"),
                     shape = NULL,
                     shape_lab = NULL,
                     shape_val = seq(21, 25, 1)) {

  if(!is.data.frame(input)) {
    stop("Input must be a data.frame or tibble with columns containing PC values and optionally fill and shape categorical variables.",
         call. = F)
  }

  if(!is.numeric(pcs) ||
     length(pcs) != 2L) {
    stop("Invalid pcs argument.",
         call. = F)
  }

  x <- paste0("PC", pcs[[1]])
  y <- paste0("PC", pcs[[2]])

  if(!all(c(x, y) %in% colnames(input)) ||
     !is.double(input[[x]]) ||
     !is.double(input[[y]])) {
    stop("Columns ", x, " and/or ", y, " not present in input or not of type double.",
         call. = F)
  }

  if((!is.null(fill) && !is.character(fill)) ||
     (is.character(fill) && (length(fill) != 1L ||
                             is.na(fill)))) {
    stop("Invalid fill argument.",
         call. = F)
  }

  if(!is.null(fill)) {

    if(!fill %in% colnames(input) ||
       ((!is.character(input[[fill]]) && !is.factor(input[[fill]])))) {
      stop("Column ", fill, " not present in input or not of type character / factor.",
           call. = F)
    }

    if((!is.null(fill_lab) && !is.character(fill_lab)) ||
       (is.character(fill_lab) && (length(fill_lab) != 1L ||
                                   is.na(fill_lab)))) {
      stop("Invalid fill_lab argument.",
           call. = F)
    }

    if(!is.character(fill_val)) {
      stop("Invalid fill_val argument.",
           call. = F)
    }

    if(length(fill_val) < length(unique(input[[fill]]))) {
      stop("Insufficient fill_val values: provided ", length(fill_val),
           ", need ", length(unique(input[[fill]])), ".",
           call. = F)
    }

  }

  if((!is.null(shape) && !is.character(shape)) ||
     (is.character(shape) && (length(shape) != 1L ||
                              is.na(shape)))) {
    stop("Invalid shape argument.",
         call. = F)
  }

  if(!is.null(shape)) {

    if((!shape %in% colnames(input) ||
        ((!is.character(input[[shape]]) && !is.factor(input[[shape]]))))) {
      stop("Column ", shape, " not present in input or not of type character / factor.",
           call. = F)
    }

    if((!is.null(shape_lab) && !is.character(shape_lab)) ||
       (is.character(shape_lab) && (length(shape_lab) != 1L ||
                                    is.na(shape_lab)))) {
      stop("Invalid shape_lab argument.",
           call. = F)
    }

    if(is.null(fill)) {

      allowed_shapes <- seq(0, 25, 1)
      n_allowed_shapes <- length(allowed_shapes)

      if(length(unique(input[[shape]])) > n_allowed_shapes) {
        stop("Variable ", shape, " has too many levels to shape-code datapoints. Maximum number of levels allowed is ", n_allowed_shapes, ".",
             call. = F)
      }

      if(!all(shape_val %in% allowed_shapes)) {
        stop("Invalid shape_val argument. Provide shape_val values between 0 and 25.",
             call. = F)
      }

    }

    if(!is.null(fill)) {

      allowed_shapes <- seq(21, 25, 1)
      n_allowed_shapes <- length(allowed_shapes)

      if(length(unique(input[[shape]])) > n_allowed_shapes) {
        stop("Variable ", shape, " has too many levels to shape-code datapoints with fill aesthetic. Maximum number of levels allowed is ", n_allowed_shapes, ".",
             call. = F)
      }

      if(!all(shape_val %in% allowed_shapes)) {
        stop("Invalid shape_val argument. Provide shape_val values between 21 and 25 to use with fill aesthetic.",
             call. = F)
      }

    }

    if(length(shape_val) < length(unique(input[[shape]]))) {
      stop("Insufficient shape_val values: provided ", length(shape_val),
           ", need ", length(unique(input[[shape]])), ".",
           call. = F)
    }

  }

  p <- ggplot2::ggplot(data = input,
                       mapping = ggplot2::aes(x = .data[[x]],
                                              y = .data[[y]])) +
    ggpubr::theme_pubr(legend = "right") +
    ggplot2::labs(x = paste0(x, ": ", round(attr(input, "pct_var")[[pcs[[1]]]], digits = 1), "% variance"),
                  y = paste0(y, ": ", round(attr(input, "pct_var")[[pcs[[2]]]], digits = 1), "% variance"))

  if(is.null(shape) && is.null(fill)) {

    p +
      ggplot2::geom_point(size = 5, color = "black", shape = 21L, fill = "#aec7e8")

  } else if(!is.null(shape) && !is.null(fill)) {

    p +
      ggplot2::geom_point(mapping = ggplot2::aes(shape = .data[[shape]],
                                                 fill = .data[[fill]]),
                          size = 5, color = "black") +
      ggplot2::scale_shape_manual(values = shape_val) +
      ggplot2::scale_fill_manual(values = fill_val) +
      ggplot2::labs(shape = shape_lab,
                    fill = fill_lab) +
      ggplot2::guides(shape = ggplot2::guide_legend(order = 1),
                      fill = ggplot2::guide_legend(order = 2,
                                                   override.aes = list(shape = 22)))

  } else if(!is.null(shape) && is.null(fill)) {

    p +
      ggplot2::geom_point(mapping = ggplot2::aes(shape = .data[[shape]]),
                          size = 5, color = "black", fill = "#aec7e8") +
      ggplot2::scale_shape_manual(values = shape_val) +
      ggplot2::labs(shape = shape_lab)

  } else {

    p +
      ggplot2::geom_point(mapping = ggplot2::aes(fill = .data[[fill]]),
                          size = 5, color = "black", shape = 21L) +
      ggplot2::scale_fill_manual(values = fill_val) +

      ggplot2::labs(fill = fill_lab) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 22)))

  }

}
