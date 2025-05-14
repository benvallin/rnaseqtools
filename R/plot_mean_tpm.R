#' Plot mean TPMs
#'
#' @param input data.frame or tibble with character columns <gene_id> and <pretty_gene_id> and double columns mean_<cnt_nm> and sem_<cnt_nm>.
#' @param genes character vector of gene IDs for which mean TPM should be be displayed. Must be of the gene ID type specified in gene_id.
#' @param gene_id character vector of length 1 or name representing gene ID. Must be a column name of input.
#' @param pretty_gene_id character vector of length 1 or name representing alternative gene ID to use for display. Must be a column name of input.
#' @param facet_var character vector of length 1 or name representing categorical variable to facet_wrap and fill-code by. Must be a column name of input.
#' @param color_var character vector of length 1 or name representing categorical variable to outline color-code by. Must be a column name of input.
#' @param expr_quartile logical vector of length 1 indicating if expression quartile labels should be displayed.
#' @param ref_gene character vector of length 1 representing the reference gene to use in labels. Ignored if expr_quartile is FALSE. mean_<cnt_nm>_pct_<ref_gene> must be a column name of input.
#' @param fill_colors character vector of colors to use for column fill based on facet_var. If facet_var is NULL, the first element of fill_colors is used.
#' @param outline_colors character vector of colors to use for column outline based on color_var. If color_var is NULL, outline_colors is ignored and black outlines are drawn.
#' @param title character vector of length 1 representing plot title.
#' @param x_label character vector of length 1 representing x axis label.
#' @param cnt_nm character vector of length 1 representing the count name used in input.
#'
#' @return a ggplot object with mean TPMs for the selected genes and optionally expression quartile labels.
#' @export
#'
#' @examples
#' # Convert TPM matrix to tibble with gene IDs as column
#' tpm_df <- tibble::as_tibble(x = bulk_tpm_ex, rownames = "ensembl_gene_id")
#'
#' # Define reference genes
#' ref_gene_ids <- gene_metadata_ex %>%
#'   dplyr::filter(gene_symbol %in% c("ACTB", "GAPDH")) %>%
#'   dplyr::pull(ensembl_gene_id, gene_symbol)
#'
#' # Compute mean TPMs per disease status
#' mean_tpm <- compute_mean_tpm(input = tpm_df,
#'                              gene_id = "ensembl_gene_id",
#'                              sample_metadata = bulk_sample_metadata_ex,
#'                              sample_id_var = "donor_id",
#'                              group_id_var = "disease_status",
#'                              ref_gene_ids = ref_gene_ids,
#'                              gene_groups = NULL,
#'                              cnt_nm = "tpm")
#'
#' # Add gene symbols to mean TPM table
#' mean_tpm <- mean_tpm %>%
#'   dplyr::mutate(disease_status = forcats::fct_relevel(disease_status, c("healthy", "diseased"))) %>%
#'   dplyr::left_join(gene_metadata_ex[, c("ensembl_gene_id", "gene_symbol")],
#'                    by = dplyr::join_by(ensembl_gene_id))
#'
#' # Define genes to display
#' genes <- gene_metadata_ex[match(x = c("SNCA", "TH", "SNAP25", "MAP2", "AQP4"),
#'                                 table = gene_metadata_ex$gene_symbol),
#'                           "ensembl_gene_id"][[1]]
#'
#' # Plot mean TPM for genes of interest
#' plot_mean_tpm(input = mean_tpm,
#'               genes = genes,
#'               gene_id = ensembl_gene_id,
#'               pretty_gene_id = gene_symbol,
#'               facet_var = disease_status,
#'               expr_quartile = TRUE,
#'               ref_gene = "ACTB",
#'               fill_colors = c("#0077b6", "#a4133c"))
#'
plot_mean_tpm <- function(input,
                          genes,
                          gene_id = "ensembl_gene_id",
                          pretty_gene_id = NULL,
                          facet_var = NULL,
                          color_var = NULL,
                          expr_quartile = FALSE,
                          ref_gene = NULL,
                          fill_colors = c("#0077b6", "#a4133c", "#d3d3d3"),
                          outline_colors = c("#5d76cb", "#ff9e00", "#29a655", "#ca3767", "#38a3a5"),
                          title = NULL,
                          x_label = "mean TPM",
                          cnt_nm = "tpm") {

  if(!is.data.frame(input)) {
    stop("Input must be a data.frame.",
         call. = F)
  }

  gene_id <- substitute(gene_id)
  gene_id <- if(is.symbol(gene_id)) deparse(gene_id) else eval(gene_id)

  if(!is.character(gene_id) ||
     length(gene_id) != 1L) {
    stop("Invalid gene_id argument.",
         call. = F)
  }

  pretty_gene_id <- substitute(pretty_gene_id)
  pretty_gene_id <- if(is.symbol(pretty_gene_id)) deparse(pretty_gene_id) else pretty_gene_id

  if(is.null(pretty_gene_id)) {

    pretty_gene_id <- gene_id

  } else {

    if(!is.character(pretty_gene_id) ||
       length(pretty_gene_id) != 1L) {
      stop("Invalid pretty_gene_id argument.",
           call. = F)
    }

  }

  if(!all(c(gene_id, pretty_gene_id) %in% colnames(input)) ||
     !is.character(input[[gene_id]]) ||
     !is.character(input[[pretty_gene_id]])) {
    stop("Columns <gene_id> and / or <pretty_gene_id> are missing or not of class character in input.",
         call. = F)
  }

  if(!is.character(genes) ||
     (is.character(genes) && (length(genes) == 0L ||
                              any(is.na(genes))))) {
    stop("Invalid genes argument.",
         call. = F)
  }

  n_genes <- length(unique(genes))
  n_genes_in_input <- sum(unique(genes) %in% input[[gene_id]], na.rm = TRUE)

  if(n_genes_in_input == 0L) {
    stop("None of the genes are present in column ", gene_id,
         " of input.\nCheck that <genes> and <gene_id> refer to the same thing!\n",
         call. = F)
  }

  if(n_genes_in_input < n_genes) {

    warning("Only ", n_genes_in_input, " / ", n_genes, " (",
            round(n_genes_in_input / n_genes * 100, 3), "%) genes found in column ",
            gene_id, " of input.\nCheck that <genes> and <gene_id> refer to the same thing!\n",
            call. = F, immediate. = T)

  }

  if(!is.character(cnt_nm) ||
     length(cnt_nm) != 1L) {
    stop("Invalid cnt_nm argument.",
         call. = F)
  }

  if(cnt_nm != "tpm") {

    input <- dplyr::rename_with(.data = input,
                                .fn = function(x) gsub(pattern = cnt_nm, replacement = "tpm", x = x),
                                .cols = c(paste0("mean_", cnt_nm),
                                          paste0("sd_", cnt_nm),
                                          paste0("sem_", cnt_nm),
                                          grep(pattern = paste0("^mean_", cnt_nm, "_pct.*$"),
                                               x = colnames(input),
                                               value = TRUE),
                                          grep(pattern = paste0("^mean_", cnt_nm, "_quartile.*$"),
                                               x = colnames(input),
                                               value = TRUE)))

  }

  if(!all(c("mean_tpm", "sem_tpm") %in% colnames(input)) ||
     !is.double(input$mean_tpm) ||
     !is.double(input$sem_tpm)) {
    stop("Columns mean_<cnt_nm> and / or sem_<cnt_nm> are missing or not of class numeric in input.",
         call. = F)
  }

  if(!is.character(fill_colors)) {
    stop("Invalid fill_colors argument.",
         call. = F)
  }

  facet_var <- substitute(facet_var)
  facet_var <- if(is.symbol(facet_var)) deparse(facet_var) else eval(facet_var)

  if(!is.null(facet_var)) {

    if(!is.character(facet_var) ||
       (is.character(facet_var) && (length(facet_var) != 1L ||
                                    is.na(facet_var)))) {
      stop("Invalid facet_var argument.",
           call. = F)
    }

    if(!facet_var %in% colnames(input) ||
       (!is.character(input[[facet_var]]) && !is.factor(input[[facet_var]]))) {
      stop("Columns <facet_var> is missing or not of class character or factor in input.",
           call. = F)
    }

    if(length(fill_colors) < length(unique(input[[facet_var]]))) {
      stop("Insufficient fill_colors values: provided ", length(fill_colors),
           ", need ", length(unique(input[[facet_var]])), ".",
           call. = F)
    }

  }

  if(is.null(facet_var) &&
     (length(input[[gene_id]]) > length(unique(input[[gene_id]])))) {
    stop("Multiple entries for one or more gene ID(s). If that is expected, use the <facet_var> argument.",
         call. = F)
  }

  color_var <- substitute(color_var)
  color_var <- if(is.symbol(color_var)) deparse(color_var) else eval(color_var)

  if(!is.null(color_var)) {

    if(!is.character(color_var) ||
       (is.character(color_var) && (length(color_var) != 1L ||
                                    is.na(color_var)))) {
      stop("Invalid color_var argument.",
           call. = F)
    }

    if(!color_var %in% colnames(input) ||
       (!is.character(input[[color_var]]) && !is.factor(input[[color_var]]))) {
      stop("Columns <color_var> is missing or not of class character or factor in input.",
           call. = F)
    }

    if(!is.character(outline_colors)) {
      stop("Invalid outline_colors argument.",
           call. = F)
    }

    if(length(outline_colors) < length(unique(input[[color_var]]))) {
      stop("Insufficient outline_colors values: provided ", length(outline_colors),
           ", need ", length(unique(input[[color_var]])), ".",
           call. = F)
    }

  }

  if(!is.logical(expr_quartile) ||
     length(expr_quartile) != 1L) {
    stop("Invalid expr_quartile argument.",
         call. = F)
  }

  if(expr_quartile) {

    if(!"mean_tpm_quartile" %in% colnames(input) ||
       !is.character(input$mean_tpm_quartile)) {
      stop("Columns mean_<cnt_nm>_quartile is missing or not of class character in input.",
           call. = F)
    }

    if(!is.null(ref_gene)) {

      if(!is.character(ref_gene) ||
         (is.character(ref_gene) && (length(ref_gene) != 1L ||
                                     is.na(ref_gene)))) {
        stop("Invalid ref_gene argument.",
             call. = F)
      }

      if(!paste0("mean_tpm_pct_", ref_gene) %in% colnames(input) ||
         !is.double(input[[paste0("mean_tpm_pct_", ref_gene)]])) {
        stop("Columns mean_<cnt_nm>_pct_<ref_gene> is missing or not of class numeric in input.",
             call. = F)
      }

    }

  }

  if((!is.null(title) && !is.character(title)) ||
     (is.character(title) && (length(title) != 1L ||
                              is.na(title)))) {
    stop("Invalid title argument.",
         call. = F)
  }

  if(!is.character(x_label) ||
     length(x_label) != 1L) {
    stop("Invalid x_label argument.",
         call. = F)
  }

  input <- input[input[[gene_id]] %in% genes,]

  gene_lvls <- unique(input[match(x = genes, table = input[[gene_id]]), pretty_gene_id][[1]])

  input[[pretty_gene_id]] <- factor(x = input[[pretty_gene_id]], levels = gene_lvls)

  aesthetics <- quote(ggplot2::aes(x = .data$mean_tpm,
                                   y = .data[[pretty_gene_id]]))

  col_layer <- quote(ggplot2::geom_col(width = 0.9, linewidth = 1))

  if(!is.null(facet_var)) {

    aesthetics[["fill"]] <- quote(.data[[facet_var]])

  } else {

    col_layer[["fill"]] <- quote(fill_colors[[1]])

  }

  if(!is.null(color_var)) {

    aesthetics[["color"]] <- quote(.data[[color_var]])

  } else {

    col_layer[["color"]] <- quote("black")

  }

  p <- ggplot2::ggplot(data = input,
                       mapping = eval(aesthetics)) +
    eval(col_layer) +
    ggplot2::geom_linerange(mapping = ggplot2::aes(xmin = .data$mean_tpm - .data$sem_tpm,
                                                   xmax = .data$mean_tpm + .data$sem_tpm))

  if(!is.null(facet_var)) {

    p <- p +
      ggplot2::facet_wrap(~ .data[[facet_var]]) +
      ggplot2::scale_fill_manual(values = fill_colors)

  }

  if(expr_quartile && !is.null(ref_gene)) {

    lab <- quote(paste0(.data$mean_tpm_quartile, " - % ", ref_gene, ": ",
                        round(.data[[paste0("mean_tpm_pct_", ref_gene)]], 3)))

  } else if(expr_quartile && is.null(ref_gene)) {

    lab <- quote(.data$mean_tpm_quartile)
  }

  if(exists("lab")) {

    p <- p +
      ggplot2::geom_label(ggplot2::aes(x = 1.2, label = eval(lab)),
                          color = "black",
                          hjust = "left",
                          size = 4,
                          fill = "#e5e6e4")

  }

  p +
    ggpubr::theme_pubr(legend = "bottom") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 18),
                   axis.text.x = ggplot2::element_text(size = 14),
                   axis.text.y = ggplot2::element_text(size = 14),
                   axis.title.x = ggplot2::element_text(size = 16),
                   strip.text.x = ggplot2::element_text(size = 14),
                   legend.text = ggplot2::element_text(size = 14)) +
    ggplot2::scale_x_log10(n.breaks = 8) +
    ggplot2::scale_color_manual(values = outline_colors) +
    ggplot2::labs(x = x_label,
                  y = NULL,
                  color = NULL,
                  title = title) +
    ggplot2::guides(fill = "none",
                    color = ggplot2::guide_legend(override.aes = list(fill = "white")))

}
