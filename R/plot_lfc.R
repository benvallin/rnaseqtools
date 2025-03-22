#' Plot log2 fold changes
#'
#' @param input data.frame or tibble with character columns <gene_id> and <pretty_gene_id>, and double columns log2FoldChange, lfcSE and padj.
#' @param genes character vector of gene IDs for which log2 fold changes should be displayed. Must be of the gene ID type specified in gene_id.
#' @param gene_id character vector of length 1 representing gene ID. Must be a column name of input.
#' @param pretty_gene_id character vector of length 1 representing alternative gene ID to use for display. Must be a column name of input.
#' @param ref_cond character vector of length 1 representing the reference biological condition.
#' @param test_cond character vector of length 1 representing the test biological condition.
#' @param key_genes character vector of gene IDs representing key genes among genes. Must be of the gene ID type specified in gene_id.
#' @param key_genes_name character vector of length 1 representing the name of the key genes group.
#' @param title character vector of length 1 representing plot title.
#' @param desc_lfc logical vector of length 1 indicating if log2 fold changes should be ordered in descending order.
#' @param xintercepts numeric vector of x positions where vertical lines should be drawn. Passed to ggplot2::geom_vline.
#' @param padj_lab logical vector of length 1 indicating if padj values should be displayed.
#' @param lab_nudge numerical vector of length 1 representing the amount of horizontal nudge to use for padj labels.
#'
#' @return a ggplot object with log2 fold changes for the selected genes in test vs reference conditions.
#' @export
#'
#' @examples
#' # Extract the gene IDs of top 25 upregulated genes with smallest padj
#' genes <- deseq2_results_ex %>%
#'   dplyr::filter(log2FoldChange > 0) %>%
#'   dplyr::arrange(padj) %>%
#'   dplyr::pull(ensembl_gene_id) %>%
#'   head(25)
#'
#' # Plot their log2 fold changes
#' plot_lfc(input = deseq2_results_ex,
#'          genes = genes,
#'          gene_id = "ensembl_gene_id",
#'          pretty_gene_id = "gene_symbol",
#'          ref_cond = "untreated",
#'          test_cond = "treated",
#'          desc_lfc = FALSE,
#'          xintercepts = seq(0.5, 2, 0.5),
#'          padj_lab = TRUE,
#'          lab_nudge = 0.02)
#'
plot_lfc <- function(input,
                     genes,
                     gene_id = "ensembl_gene_id",
                     pretty_gene_id = "gene_symbol",
                     ref_cond = "ref cond.",
                     test_cond = "test cond.",
                     key_genes = NULL,
                     key_genes_name = "leading edge",
                     title = NULL,
                     desc_lfc = TRUE,
                     xintercepts = NULL,
                     padj_lab = TRUE,
                     lab_nudge = 0.02) {

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
     !all(c("log2FoldChange", "lfcSE", "padj") %in% colnames(input)) ||
     !is.character(input[[gene_id]]) ||
     !is.character(input[[pretty_gene_id]]) ||
     !is.double(input$log2FoldChange) ||
     !is.double(input$lfcSE) ||
     !is.double(input$padj)) {
    stop("Input must be a data.frame or tibble with character columns <gene_id> and <pretty_gene_id>, and double columns log2FoldChange, lfcSE and padj.",
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

  if(!is.null(key_genes) &&
     !is.character(key_genes)) {
    stop("Invalid key_genes argument.",
         call. = F)
  }

  if(!is.character(key_genes_name) ||
     length(key_genes_name) != 1L) {
    stop("Invalid key_genes_name argument.",
         call. = F)
  }

  if((!is.null(title) && !is.character(title)) ||
     (is.character(title) && (length(title) != 1L ||
                              is.na(title)))) {
    stop("Invalid title argument.",
         call. = F)
  }

  if(!is.logical(desc_lfc) ||
     length(desc_lfc) != 1L) {
    stop("Invalid desc_lfc argument.",
         call. = F)
  }

  if((!is.null(xintercepts) && !is.numeric(xintercepts)) ||
     (is.numeric(xintercepts) && (length(xintercepts) == 0L ||
                                  any(is.na(xintercepts))))) {
    stop("Invalid xintercepts argument.",
         call. = F)
  }

  if(!is.logical(padj_lab) ||
     length(padj_lab) != 1L) {
    stop("Invalid padj_lab argument.",
         call. = F)
  }

  if(!is.numeric(lab_nudge) ||
     length(lab_nudge) != 1L) {
    stop("Invalid lab_nudge argument.",
         call. = F)
  }

  df <- input

  df <- df[df[[gene_id]] %in% genes & !is.na(df$log2FoldChange),]

  df <- df %>%
    dplyr::mutate(padj_cat = dplyr::case_when(.data$padj < 0.01 ~ "padj < 0.01",
                                              .data$padj < 0.05 ~ "padj < 0.05",
                                              .data$padj < 0.1 ~ "padj < 0.1",
                                              .data$padj >= 0.1 ~ "padj \u2265 0.1",
                                              is.na(.data$padj) ~ "padj not computed"),
                  padj_color = dplyr::case_when(.data$padj < 0.01 ~ "#e0007f",
                                                .data$padj < 0.05 ~ "#b892ff",
                                                .data$padj < 0.1 ~ "#29a655",
                                                .data$padj >= 0.1 ~ "#ff9100",
                                                is.na(.data$padj) ~ "grey"))

  p <- ggplot2::ggplot(data = df,
                       mapping = ggplot2::aes(x = .data$log2FoldChange,
                                              y = forcats::fct_reorder(.data[[pretty_gene_id]],
                                                                       .data$log2FoldChange,
                                                                       .desc = desc_lfc),
                                              fill = .data$padj_cat)) +
    ggplot2::geom_hline(yintercept = df[[pretty_gene_id]],
                        color = "grey",
                        alpha = 0.2) +
    ggplot2::geom_vline(xintercept = 0,
                        linetype = "dashed",
                        color = "#c9184a",
                        linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = xintercepts,
                        linetype = "dashed",
                        color = "orange",
                        linewidth = 0.5) +
    ggpubr::theme_pubr(legend = "bottom") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 12),
                   axis.text.y = ggplot2::element_text(size = 12),
                   axis.title.x = ggplot2::element_text(size = 14),
                   legend.text = ggplot2::element_text(size = 14),
                   strip.text.x = ggplot2::element_text(size = 14)) +
    ggplot2::scale_fill_manual(values = df %>%
                                 dplyr::arrange(.data$padj) %>%
                                 dplyr::pull(.data$padj_color) %>%
                                 unique()) +
    ggplot2::geom_linerange(mapping = ggplot2::aes(xmin = .data$log2FoldChange - .data$lfcSE,
                                                   xmax = .data$log2FoldChange + .data$lfcSE)) +
    ggplot2::labs(title = title,
                  x = bquote(log[2] * " fold change (" * .(test_cond) * " vs " * .(ref_cond) * ")"),
                  y = NULL,
                  fill = NULL,
                  shape = NULL)

  if(!is.null(xintercepts)) {

    p <- p +
      ggplot2::scale_x_continuous(breaks = c(0, xintercepts))

  }

  if(padj_lab) {

    p <- p +
      ggplot2::geom_text(mapping = ggplot2::aes(x = .data$log2FoldChange + .data$lfcSE + lab_nudge,
                                                label = ifelse(is.na(.data$padj),
                                                               "",
                                                               formatC(.data$padj, format = "e", digits = 1))),
                         hjust = 0)

  }

  if(is.null(key_genes)) {

    p +
      ggplot2::geom_point(shape = 21,
                          color = "black",
                          size = 5,
                          stroke = 1)

  } else {

    df$in_key_genes <- ifelse(df[[gene_id]] %in% key_genes,
                              paste0("in ", key_genes_name),
                              paste0("not in ", key_genes_name))

    if(length(unique(df$in_key_genes)) == 2) {

      shapes <- c(23, 21)

    } else if(unique(df$in_key_genes) == paste0("in ", key_genes_name)) {

      shapes <- 23

    } else {

      shapes <- 21

    }

    p +
      ggplot2::geom_point(data = df,
                          ggplot2::aes(shape = .data$in_key_genes),
                          color = "black",
                          size = 5,
                          stroke = 1) +
      ggplot2::scale_shape_manual(values = shapes) +
      ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 22, size = 5)))

  }
}
