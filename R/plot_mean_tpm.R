plot_mean_tpm <- function(input,
                          genes,
                          gene_id = "ensembl_gene_id",
                          pretty_gene_id = "gene_symbol",
                          facet_var = NULL,
                          color_var = NULL,
                          expr_quartile = FALSE,
                          ref_gene = NULL,
                          fill_colors = c("#0077b6", "#a4133c", "#d3d3d3"),
                          outline_colors = c("#5d76cb", "#ff9e00", "#29a655", "#ca3767", "#38a3a5"),
                          title = NULL,
                          x_label = "mean TPM",
                          cnt_nm = "tpm") {

  gene_id <- substitute(gene_id)
  gene_id <- ifelse(is.symbol(gene_id), deparse(gene_id), gene_id)

  pretty_gene_id <- substitute(pretty_gene_id)
  pretty_gene_id <- ifelse(is.symbol(pretty_gene_id), deparse(pretty_gene_id), pretty_gene_id)

  if(is.null(facet_var) &&
     (length(input[[gene_id]]) > length(unique(input[[gene_id]])))) {
    stop("Multiple entries for one or more gene ID(s). If that is expected, use the <facet_var> argument.",
         call. = F)
  }

  df <- input

  df <- df[df[[gene_id]] %in% genes,]

  if(cnt_nm != "tpm") {

    df <- dplyr::rename_with(.data = df,
                             .fn = function(x) gsub(pattern = cnt_nm, replacement = "tpm", x = x),
                             .cols = c(paste0("mean_", cnt_nm),
                                       paste0("sd_", cnt_nm),
                                       paste0("sem_", cnt_nm),
                                       grep(pattern = paste0("^mean_", cnt_nm, "_pct.*$"),
                                            x = colnames(df),
                                            value = TRUE),
                                       grep(pattern = paste0("^mean_", cnt_nm, "_quartile.*$"),
                                            x = colnames(df),
                                            value = TRUE)))

  }

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

  p <- ggplot2::ggplot(data = df,
                       mapping = eval(aesthetics)) +
    eval(col_layer) +
    ggplot2::geom_linerange(mapping = ggplot2::aes(xmin = .data$mean_tpm - .data$sem_tpm,
                                                   xmax = .data$mean_tpm + .data$sem_tpm))

  if(!is.null(facet_var)) {

    p <- p +
      ggplot2::scale_fill_manual(values = fill_colors) +
      ggplot2::facet_wrap(~ .data[[facet_var]])

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
