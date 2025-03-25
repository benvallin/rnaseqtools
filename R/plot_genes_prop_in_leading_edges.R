#' Plot genes proportions in leading edges
#'
#' @param input data.frame or tibble with columns collection and n_leadingEdge and list-column leadingEdge, as produced by rnaseqtools::multi_fgsea.
#' @param genes character vector of gene IDs to search in leading edges. Must be the same type of gene IDs as in the leadingEdge list-column.
#' @param x_lab character vector of length 1 representing the label of x axis.
#' @param gene_group_name character vector of length 1 representing the name of the genes group. Will be used in y axis label.
#' @param collection_levels character vector with collection levels in desired order.
#' @param collection_colors character vector with collection level colors in desired order.
#'
#' @return a ggplot object with genes proportions in leading edges for the selected genes.
#' @export
#'
#' @examples
#' # Make MSigDB collection table
#' msigdb_collection_table = get_msigdb_collections()
#'
#' # Define collection of interest
#' collection_name <- "MSigDB_C2_CP:KEGG_LEGACY"
#'
#' # Define gene set name of interest
#' set_name <- "KEGG_RIBOSOME"
#'
#' # Get genes in gene set of interest
#' genes <- get_genes_in_collection(input = msigdb_collection_table,
#'                                  collection = collection_name,
#'                                  set_name = set_name,
#'                                  gene_ids = "ensembl_gene_id")
#'
#' genes <- genes$ensembl_gene_id
#'
#' # Plot genes proportions in leading edges for the selected genes
#' plot_genes_prop_in_leading_edges(input = multi_fgsea_results_ex,
#'                                  genes = genes,
#'                                  x_lab = "gene set",
#'                                  gene_group_name = set_name,
#'                                  collection_levels = c("MSigDB_C2_CP:REACTOME",
#'                                                        "MSigDB_C2_CP:KEGG_LEGACY",
#'                                                        "MSigDB_C5_GO:BP"),
#'                                  collection_colors = c("#5d76cb",
#'                                                        "#29a655",
#'                                                        "#ca3767"))
#'
plot_genes_prop_in_leading_edges <- function(input,
                                             genes,
                                             x_lab = "gene set",
                                             gene_group_name = "searched",
                                             collection_levels = NULL,
                                             collection_colors = NULL) {

  if(!is.character(x_lab) ||
     length(x_lab) != 1L) {
    stop("Invalid x_lab argument.",
         call. = F)
  }

  if(!is.character(gene_group_name) ||
     length(gene_group_name) != 1L) {
    stop("Invalid gene_group_name argument.",
         call. = F)
  }

  df <- search_genes_in_leading_edges(input = input,
                                      genes = genes,
                                      group_name = "searched")

  if(!is.null(collection_levels) &&
     !is.null(collection_colors)) {

    df <- set_level_color(input = df,
                          variable = "collection",
                          levels = collection_levels,
                          colors = collection_colors)

  }

  x_pos_quantiles <- c(0.25*dim(df)[[1]],
                       0.5*dim(df)[[1]],
                       0.75*dim(df)[[1]])

  p <- ggplot2::ggplot(data = df,
                       mapping = ggplot2::aes(x = .data$set_number,
                                              y = .data$pct_leadingEdge_searched_in_leadingEdge)) +
    ggplot2::geom_vline(xintercept = df$set_number,
                        color = "grey",
                        alpha = 0.2) +
    ggplot2::geom_vline(xintercept = x_pos_quantiles,
                        linetype = "dashed",
                        color = "orange",
                        linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = c(0, 25, 50, 75, 100),
                        linetype = "dotted",
                        color = "#c9184a",
                        linewidth = 0.5) +
    ggplot2::geom_label(data = tibble::tibble(x = x_pos_quantiles,
                                              y = 105,
                                              label = c("Q1", "Q2", "Q3")),
                        ggplot2::aes(x = .data$x,
                                     y = .data$y,
                                     label = .data$label)) +
    ggplot2::geom_point(ggplot2::aes(fill = .data$collection),
                        shape = 21,
                        color = "black",
                        size = 3,
                        stroke = 1,
                        alpha = 0.8) +
    ggpubr::theme_pubr(legend = "bottom") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 14),
                   axis.title = ggplot2::element_text(size = 16),
                   legend.text = ggplot2::element_text(size = 16)) +
    ggplot2::scale_y_continuous(limits = c(0L, 105L),
                                breaks = seq(0, 100, 10),
                                labels = seq(0, 100, 10)) +
    ggplot2::labs(x = x_lab,
                  y = paste0("Proportion of leading edge corresponding to ", gene_group_name, " genes (%)"),
                  fill = NULL)

  if(!is.null(collection_levels) &&
     !is.null(collection_colors)) {

    p +
      ggplot2::scale_color_manual(values = levels(df$collection_color), aesthetics = "fill")

  } else {

    p

  }

}

