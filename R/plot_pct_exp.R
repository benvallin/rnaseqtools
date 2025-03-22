#' Plot percentage expression
#'
#' @param input data.frame or tibble with character columns <gene_id> and <pretty_gene_id>, and double column log2FoldChange.
#' @param genes character vector of gene IDs for which percentage expression in test vs reference conditions should be displayed. Must be of the gene ID type specified in gene_id.
#' @param gene_id character vector of length 1 representing gene ID. Must be a column name of input.
#' @param pretty_gene_id character vector of length 1 representing alternative gene ID to use for display. Must be a column name of input.
#' @param ref_cond character vector of length 1 representing the reference biological condition.
#' @param test_cond character vector of length 1 representing the test biological condition.
#' @param title character vector of length 1 representing plot title.
#' @param fill_colors character vector of length 2 representing the colors to use for reference and test conditions, respectively.
#'
#' @return a ggplot object with percentage expression in test relative to reference conditions for the selected genes.
#' @export
#'
#' @examples
#' # Extract the gene IDs of top 25 downregulated genes with smallest padj
#' genes <- deseq2_results_ex %>%
#'   dplyr::filter(log2FoldChange < 0) %>%
#'   dplyr::arrange(padj) %>%
#'   dplyr::pull(ensembl_gene_id) %>%
#'   head(25)
#'
#' # Plot their percentage expression in test relative to reference conditions
#' plot_pct_exp(input = deseq2_results_ex,
#'              genes = genes,
#'              gene_id = "ensembl_gene_id",
#'              pretty_gene_id = "gene_symbol",
#'              ref_cond = "untreated",
#'              test_cond = "treated",
#'              title = NULL,
#'              fill_colors = c("#5d76cb", "#ca3767"))
#'
plot_pct_exp <- function(input,
                         genes,
                         gene_id = "ensembl_gene_id",
                         pretty_gene_id = "gene_symbol",
                         ref_cond = "ref cond.",
                         test_cond = "test cond.",
                         title = NULL,
                         fill_colors = c("#ca3767", "#5d76cb")) {

  if(!is.character(gene_id) ||
     length(gene_id) != 1L ||
     !gene_id %in% colnames(input)) {
    stop("Invalid gene_id argument.",
         call. = F)
  }

  if(!is.character(pretty_gene_id) ||
     length(pretty_gene_id) != 1L ||
     !pretty_gene_id %in% colnames(input)) {
    stop("Invalid pretty_gene_id argument.",
         call. = F)
  }

  if(!is.data.frame(input) ||
     !"log2FoldChange" %in% colnames(input) ||
     !is.character(input[[gene_id]]) ||
     !is.character(input[[pretty_gene_id]]) ||
     !is.double(input$log2FoldChange)) {
    stop("Input must be a data.frame or tibble with character columns <gene_id> and <pretty_gene_id>, and double column log2FoldChange.",
         call. = F)
  }

  if(!is.character(genes) ||
     (is.character(genes) && (length(genes) == 0L ||
                              any(is.na(genes))))) {
    stop("Invalid genes argument.",
         call. = F)
  }

  if(!is.character(test_cond) ||
     length(test_cond) != 1L) {
    stop("Invalid test_cond argument.",
         call. = F)
  }

  if(!is.character(ref_cond) ||
     length(ref_cond) != 1L) {
    stop("Invalid ref_cond argument.",
         call. = F)
  }

  if((!is.null(title) && !is.character(title)) ||
     (is.character(title) && (length(title) != 1L ||
                              is.na(title)))) {
    stop("Invalid title argument.",
         call. = F)
  }

  if(!is.character(fill_colors) ||
     length(fill_colors) != 2L) {
    stop("Invalid fill_colors argument.",
         call. = F)
  }

  df <- input

  df <- df[df[[gene_id]] %in% genes & !is.na(df$log2FoldChange),]


  df[[pretty_gene_id]] <- forcats::fct_relevel(df[[pretty_gene_id]],
                                               dplyr::arrange(.data = df,
                                                              .data$log2FoldChange) %>%
                                                 dplyr::pull(.data[[pretty_gene_id]]))

  df$reference <- 100L

  df$test <- (2^(df$log2FoldChange))*100

  df <- df[, c(gene_id, pretty_gene_id, "reference", "test")]

  df <- df %>%
    tidyr::pivot_longer(cols = c("reference", "test"),
                        names_to = "condition",
                        values_to = "pct_exp") %>%
    dplyr::mutate(condition = dplyr::case_when(.data$condition == "reference" ~ ref_cond,
                                               .data$condition == "test" ~ test_cond,
                                               TRUE ~ NA_character_) %>%
                    forcats::fct_relevel(c(ref_cond, test_cond)))

  p <- ggplot2::ggplot()

  if(mean(df$pct_exp) < 100) {

    p <- p +
      ggplot2::geom_col(data = df[df$condition == ref_cond,],
                        mapping = ggplot2::aes(x = .data[[pretty_gene_id]],
                                               y = .data$pct_exp,
                                               fill = .data$condition),
                        color = "black", position = "dodge", width = 0.9, alpha = 0.9, linewidth = 0.8) +
      ggplot2::geom_col(data = df[df$condition == test_cond,],
                        mapping = ggplot2::aes(x = .data[[pretty_gene_id]],
                                               y = .data$pct_exp,
                                               fill = .data$condition),
                        color = "black", position = "dodge", width = 0.9, alpha = 0.9, linewidth = 0.8) +
      ggplot2::scale_y_continuous(limits = c(0L, 100L), n.breaks = 10L)

  } else {

    p <- p +
      ggplot2::geom_col(data = df[df$condition == test_cond,],
                        mapping = ggplot2::aes(x = .data[[pretty_gene_id]],
                                               y = .data$pct_exp,
                                               fill = .data$condition),
                        color = "black", position = "dodge", width = 0.9, alpha = 0.9, linewidth = 0.8) +
      ggplot2::geom_col(data = df[df$condition == ref_cond,],
                        mapping = ggplot2::aes(x = .data[[pretty_gene_id]],
                                               y = .data$pct_exp,
                                               fill = .data$condition),
                        color = "black", position = "dodge", width = 0.9, alpha = 0.9, linewidth = 0.8) +
      ggplot2::scale_y_continuous(limits = c(0L, ceiling(max(df$pct_exp))), n.breaks = 12L)

  }

  p +
    ggpubr::theme_pubr(legend = "bottom") +
    ggplot2::scale_fill_manual(limits = c(ref_cond, test_cond),
                               values = fill_colors) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12),
                   axis.title = ggplot2::element_text(size = 14),
                   legend.text = ggplot2::element_text(size = 14),
                   legend.title = ggplot2::element_blank()) +
    ggpubr::rotate_x_text(angle = 30) +
    ggplot2::labs(x = NULL,
                  y = "% expression",
                  title = title)

}
