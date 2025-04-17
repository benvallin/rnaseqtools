#' Plot MSA
#'
#' @param input data.frame or tibble with character columns name and sequence, as produced by rnaseqtools::msa_to_tibble.
#'
#' @return a ggplot object with multiple sequence alignment.
#' @export
#'
#' @examples
#' # Extract multiple sequence alignment from MsaAAMultipleAlignment object
#' msa_df <- msa_to_tibble(input = msa_ex)
#'
#' # Plot multiple sequence alignment
#' plot_msa(input = msa_df)
#'
plot_msa <- function(input) {

  if(!requireNamespace("viridis", quietly = TRUE)) {

    stop("Package \"viridis\" must be installed to use this function.",
         call. = F)

  }

  if(!is.data.frame(input) ||
     !all(c("name", "sequence") %in% colnames(input)) ||
     !is.character(input$name) ||
     !is.character(input$sequence)) {
    stop("Input must be a data.frame or tibble with character columns name and sequence, as produced by rnaseqtools::msa_to_tibble.",
         call. = F)
  }

  df <- input

  df$name <- factor(x = df$name, levels = rev(df$name))

  df$sequence <- strsplit(x = df$sequence, split = "")

  df$position <- lapply(X = df$sequence,
                        FUN = seq_along)

  temp <- do.call(what = cbind, args = df$sequence)

  npos <- dim(temp)[[1]]
  nseq <- dim(temp)[[2]]

  frequency <- vector(mode = "list", length = npos)

  for(i in seq(npos)) {

    vals <- temp[i,]

    frequency[[i]] <- purrr::map_dbl(.x = vals,
                                     .f = function(x) {

                                       if(x == "-") {
                                         NA_real_
                                       } else {
                                         sum(vals == x) / nseq
                                       }

                                     })

  }

  frequency <- do.call(what = cbind, args = frequency)

  df$frequency <- asplit(x = frequency, MARGIN = 1)

  df <- tidyr::unnest(data = df,
                      cols = c("sequence", "position", "frequency"))

  df$name_pos <- as.integer(df$name)

  df$position <- as.factor(df$position)

  ggplot2::ggplot(data = df,
                  mapping = ggplot2::aes(x = .data$position)) +
    ggplot2::geom_col(mapping = ggplot2::aes(y = .data$name_pos,
                                             fill = .data$frequency),
                      color = "white",
                      position = ggplot2::position_identity()) +
    ggplot2::geom_text(mapping = ggplot2::aes(y = .data$name_pos-0.5,
                                              label = .data$sequence),
                       color = "white") +
    ggplot2::scale_y_continuous(breaks = seq(max(df$name_pos))-0.5,
                                labels = levels(df$name),
                                expand = c(0.001, 0.001)) +
    viridis::scale_fill_viridis(option = "magma",
                                direction = -1,
                                end = 0.8,
                                na.value = "#edede9") +
    ggpubr::theme_pubr(legend = "right") +
    ggplot2::theme(plot.margin = ggplot2::margin(),
                   legend.justification = c(1, 1),
                   axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 5)) +
    ggplot2::labs(x = NULL,
                  y = NULL) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(theme = ggplot2::theme(legend.title = ggplot2::element_text(margin = ggplot2::margin(b = 15)),
                                                                          legend.frame = ggplot2::element_rect(linewidth = 0.2, color = "black"),
                                                                          legend.ticks = ggplot2::element_blank())))

}
