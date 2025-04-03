#' Plot averaged rescaled TPMs
#'
#' @param input tibble produced by rnaseqtools::average_rescaled_tpm or rnaseqtools::average_sum_rescaled_tpm.
#' @param feature_type character vector of length 1 or name representing feature type. Must be one of genes or genesets.
#' @param feature_id character vector of length 1 or name representing feature ID. Must be a column name of input.
#' @param grp_id character vector of length 1 or name representing sample groups to display on x axis. Must be a column name of input.
#' @param error_bar character vector of length 1 or name representing metric to use for drawing error bars. Must be one of sd or sem.
#' @param y_label character vector of length 1 representing x axis label.
#' @param title character vector of length 1 representing plot title.
#' @param colors character vector of colors to use for column fill based on <feature_id>.
#'
#' @return a ggplot object with averaged (summed) rescaled TPMs split by sample group. Labels express feature abundance as a percentage of all features for each group.
#' @export
#'
#' @examples
#' # Define genes of interest
#' genes <- calcium_genes_ex %>%
#'   dplyr::filter(protein_complex == "SERCA") %>%
#'   dplyr::pull(ensembl_feature_id)
#'
#' # Rescale TPMs for selected genes
#' resc_tpm <- rescale_tpm(input = bulk_tpm_ex, genes = genes)
#'
#' # Average rescaled TPMs by disease status using a grp_df providing sample metadata
#' mean_resc_tpm <- average_rescaled_tpm(input = resc_tpm,
#'                                       grp_df = bulk_sample_metadata_ex,
#'                                       sample_id = donor_id,
#'                                       feature_id = ensembl_feature_id)
#'
#' # Add gene symbols to averaged rescaled TPMs
#' mean_resc_tpm <- mean_resc_tpm %>%
#'   dplyr::mutate(disease_status = forcats::fct_relevel(disease_status, c("healthy", "diseased"))) %>%
#'   dplyr::left_join(gene_metadata_ex[, c("ensembl_feature_id", "gene_symbol")],
#'                    by = dplyr::join_by(ensembl_feature_id))
#'
#' # Plot averaged rescaled TPMs by disease status
#' plot_rescaled_tpm(input = mean_resc_tpm,
#'                   feature_type = "genes",
#'                   feature_id = gene_symbol,
#'                   grp_id = disease_status,
#'                   error_bar = sd,
#'                   title = "SERCA protein complex")
#'
plot_rescaled_tpm <- function(input,
                              feature_type = "genes",
                              feature_id,
                              grp_id,
                              error_bar = "sd",
                              y_label = NULL,
                              title = NULL,
                              colors = c("#5d76cb", "#ff9e00", "#ca3767", "#38a3a5", "#f49cbb", "#858ae3", "#ffd500","#2ba84a", "#fc6dab", "#7637da"))
{

  feature_type <- substitute(feature_type)
  feature_type <- if(is.symbol(feature_type)) deparse(feature_type) else eval(feature_type)

  feature_id <- substitute(feature_id)
  feature_id <- if(is.symbol(feature_id)) deparse(feature_id) else eval(feature_id)

  grp_id <- substitute(grp_id)
  grp_id <- if(is.symbol(grp_id)) deparse(grp_id) else eval(grp_id)

  error_bar <- substitute(error_bar)
  error_bar <- if(is.symbol(error_bar)) deparse(error_bar) else eval(error_bar)

  if(!is.character(feature_type) ||
     length(feature_type) != 1L ||
     !feature_type %in% c("genes", "genesets")) {
    stop("feature_type must be one of genes or genesets.",
         call. = F)
  }

  if(feature_type == "genes") {

    rescaled_tpm <- "rescaled_tpm"
    mean_rescaled_tpm <- "mean_rescaled_tpm"

  } else {

    rescaled_tpm <- "sum_rescaled_tpm"
    mean_rescaled_tpm <- "mean_sum_rescaled_tpm"

  }

  if(!is.character(feature_id) ||
     length(feature_id) != 1L) {
    stop("Invalid feature_id argument.",
         call. = F)
  }

  if(!is.character(grp_id) ||
     length(grp_id) != 1L) {
    stop("Invalid grp_id argument.",
         call. = F)
  }

  if(!is.character(error_bar) ||
     length(error_bar) != 1L ||
     !error_bar %in% c("sd", "sem")) {
    stop("error_bar must be one of sd or sem.",
         call. = F)
  }

  error_bar_var <- paste0(error_bar, "_", rescaled_tpm)

  if(!is.data.frame(input) ||
     !all(c(feature_id, grp_id, rescaled_tpm, mean_rescaled_tpm, error_bar_var) %in% colnames(input)) ||
     (!is.character(input[[feature_id]]) && !is.factor(input[[feature_id]])) ||
     (!is.character(input[[grp_id]]) && !is.factor(input[[grp_id]])) ||
     !is.numeric(input[[rescaled_tpm]]) ||
     !is.numeric(input[[mean_rescaled_tpm]]) ||
     !is.numeric(input[[error_bar_var]])) {
    stop("input must be a data.frame or tibble with character / factor columns <feature_id> and <grp_id> and numeric columns ",
         rescaled_tpm, ", ", mean_rescaled_tpm, " and ", paste0("<error_bar>_", rescaled_tpm), ".",
         call. = F)
  }

  temp <- unique(input[, c(feature_id, grp_id)])

  if(!all(purrr::map_lgl(.x = split(x = temp, f = temp[[grp_id]]),
                         .f = ~ length(unique(.x[[feature_id]])) == length(.x[[feature_id]])))) {
    stop("Feature IDs in column <feature_id> of input are not all unique for every group.",
         call. = F)
  }

  temp <- unique(input[, c(feature_id, grp_id, mean_rescaled_tpm)]) %>%
    dplyr::group_by(dplyr::across(.cols = c(feature_id, grp_id))) %>%
    dplyr::add_count() %>%
    dplyr::pull(.data$n) %>%
    unique()

  if(length(temp) != 1L ||
     temp != 1L) {
    stop("There should be a single ", mean_rescaled_tpm, " value for each combination of <feature_id> and <grp_id>.",
         call. = F)
  }

  if(is.null(y_label)) {

    y_label <- paste0("% of gene set",
                      ifelse(feature_type == "genes", "", " universe"))

  } else {

    if(!is.character(y_label) ||
       length(y_label) != 1L) {
      stop("Invalid y_label argument.",
           call. = F)
    }

  }

  if((!is.null(title) && !is.character(title)) ||
     (is.character(title) && (length(title) != 1L ||
                              is.na(title)))) {
    stop("Invalid title argument.",
         call. = F)
  }

  if(!is.null(colors)) {

    if(!is.character(colors)) {
      stop("Invalid colors argument.",
           call. = F)
    }

    n_features <- length(unique(input[[feature_id]]))

    colors <- unique(colors)

    n_colors <- length(colors)

    if(n_colors < n_features) {
      stop("Need ", n_features-n_colors, " more colors for those ", n_features, " ", feature_type, ".",
           call. = F)
    }

    colors <- factor(x = colors, levels = colors)

    colors <- colors[1:n_features]

    colors <- as.character(sort(x = colors, decreasing = T))

  }

  df <- unique(input[, c(feature_id, grp_id, rescaled_tpm, mean_rescaled_tpm, error_bar_var)]) %>%
    dplyr::group_by(.data[[feature_id]]) %>%
    dplyr::mutate(overall_mean_rescaled_tpm = mean(.data[[rescaled_tpm]], na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$overall_mean_rescaled_tpm) %>%
    dplyr::mutate(!!feature_id := factor(.data[[feature_id]],
                                         levels = unique(.data[[feature_id]]))) %>%
    dplyr::arrange(dplyr::desc(.data[[feature_id]])) %>%
    dplyr::select(grp_id, feature_id, mean_rescaled_tpm, error_bar_var) %>%
    unique() %>%
    dplyr::group_by(.data[[grp_id]]) %>%
    dplyr::mutate(stacked_y = cumsum(.data[[mean_rescaled_tpm]])) %>%
    dplyr::ungroup()

  p <- ggplot2::ggplot(data = df,
                       mapping = ggplot2::aes(x = .data[[grp_id]],
                                              y = .data[[mean_rescaled_tpm]])) +
    ggplot2::geom_col(mapping = ggplot2::aes(fill = .data[[feature_id]]),
                      color = "black") +
    ggplot2::geom_errorbar(mapping = ggplot2::aes(ymin = .data$stacked_y - .data[[error_bar_var]],
                                                  ymax = .data$stacked_y + .data[[error_bar_var]]),
                           position = "identity",
                           color = "black", linewidth = 0.5, width = 0.2) +
    ggrepel::geom_label_repel(mapping = ggplot2::aes(label = paste0(round(.data[[mean_rescaled_tpm]], 2), " %"),
                                                     color = .data[[feature_id]]),
                              position = ggplot2::position_stack(vjust = 0.5),
                              fill = "white", fontface = "bold", size = 4,
                              max.overlaps = Inf,
                              point.size = 0,
                              show.legend = F) +
    ggplot2::scale_y_continuous(breaks = seq(0, 100, 5)) +
    ggpubr::theme_pubr(legend = "right") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 16),
                   axis.text.x = ggplot2::element_text(size = 14),
                   axis.text.y = ggplot2::element_text(size = 12),
                   legend.text = ggplot2::element_text(size = 14),
                   panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                   plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
                   legend.background = ggplot2::element_rect(fill = "transparent", colour = NA)) +
    ggplot2::labs(x = NULL,
                  y = y_label,
                  title = title,
                  fill = NULL,
                  color = NULL)

  if(!is.null(colors)) {

    p <- p +
      ggplot2::scale_fill_manual(values = colors,
                                 aesthetics = c("fill", "color"))

  }

  p

}
