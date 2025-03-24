#' Compute mean TPMs
#'
#' @param input data.frame or tibble with character columns <gene_id> and double columns representing sample-specific gene counts.
#' @param gene_id character vector of length 1 representing gene ID. Must be a column name of input.
#' @param sample_metadata data.frame or tibble with character / factor columns <sample_id_var> and <group_id_var>.
#' @param sample_id_var character vector of length 1 representing sample ID. Must be a column name of sample_metadata.
#' @param group_id_var character vector of length 1 representing biological group ID. Must be a column name of sample_metadata.
#' @param ref_gene_ids named character vector of gene IDs to use as reference. Values must be of the gene ID type specified in gene_id and names will be used in column names of output.
#' @param gene_groups list of named character vectors representing gene groups to use as reference.
#' @param cnt_nm ss
#'
#' @return ss
#' @export
#'
#' @examples
compute_mean_tpm <- function(input,
                             gene_id = "ensembl_gene_id_version",
                             sample_metadata,
                             sample_id_var,
                             group_id_var,
                             ref_gene_ids = NULL,
                             gene_groups = NULL,
                             cnt_nm = "tpm") {


  if(!is.character(gene_id) ||
     length(gene_id) != 1L) {
    stop("Invalid gene_id argument.",
         call. = F)
  }

  if(!is.character(sample_id_var) ||
     length(sample_id_var) != 1L) {
    stop("Invalid sample_id_var argument.",
         call. = F)
  }

  if(!is.character(group_id_var) ||
     length(group_id_var) != 1L) {
    stop("Invalid group_id_var argument.",
         call. = F)
  }

  if(!is.character(cnt_nm) ||
     length(cnt_nm) != 1L) {
    stop("Invalid cnt_nm argument.",
         call. = F)
  }

  if(!is.data.frame(input) ||
     !gene_id %in% colnames(input) ||
     !is.character(input[[gene_id]]) ||
     !all(purrr::map_lgl(.x = setdiff(colnames(input), gene_id),
                         .f = ~ is.double(input[[.x]])))) {
    stop("Input must be a data.frame or tibble with only character column <gene_id> and double columns representing sample-specific gene counts.",
         call. = F)
  }

  if(!is.data.frame(sample_metadata) ||
     !all(c(sample_id_var, group_id_var) %in% colnames(sample_metadata)) ||
     (!is.character(sample_metadata[[sample_id_var]]) && !is.factor(sample_metadata[[sample_id_var]])) ||
     (!is.character(sample_metadata[[group_id_var]]) && !is.factor(sample_metadata[[group_id_var]]))) {

    stop("sample_metadata must be a data.frame or tibble with character / factor columns <sample_id_var> and <group_id_var>.",
         call. = F)
  }

  if(!all(setdiff(colnames(input), gene_id) %in% unique(as.character(sample_metadata[[sample_id_var]])))) {
    stop("The names of double columns in input must all be present in column <sample_id_var> of sample_metadata.",
         call. = F)
  }

  if((!is.null(ref_gene_ids) && !is.character(ref_gene_ids)) ||
     (is.character(ref_gene_ids) && (any(is.na(ref_gene_ids)) ||
                                     is.null(names(ref_gene_ids))))) {
    stop("Invalid ref_gene_ids argument.",
         call. = F)
  }

  if((!is.null(gene_groups) && !is.list(gene_groups)) ||
     (is.list(gene_groups) && (!all(purrr::map_lgl(.x = gene_groups, .f = is.character)) ||
                               any(is.na(gene_groups)) ||
                               is.null(names(gene_groups))))) {
    stop("Invalid gene_groups argument.",
         call. = F)
  }

  output <- tidyr::pivot_longer(data = input,
                                cols = setdiff(colnames(input), gene_id),
                                names_to = sample_id_var,
                                values_to = "tpm")

  output <- dplyr::left_join(x = output,
                             y = sample_metadata[, c(sample_id_var, group_id_var)],
                             by = sample_id_var)

  output <- tidyr::nest(.data = output,
                        sample_data = tidyselect::all_of(c(sample_id_var, "tpm")))

  output$mean_tpm <- purrr::map_dbl(.x = output$sample_data,
                                    .f = ~ mean(.x$tpm, na.rm = T))

  output$sd_tpm <- purrr::map_dbl(.x = output$sample_data,
                                  .f = ~ sd(.x$tpm, na.rm = T))

  output$sem_tpm <- purrr::map2_dbl(.x = output$sd_tpm,
                                    .y = output$sample_data,
                                    .f = ~ .x / dim(.y)[[1]])

  output <- tidyr::nest(.data = output,
                        .by = tidyselect::all_of(group_id_var),
                        .key = "group_data")

  output$group_data <- purrr::map(.x = output$group_data,
                                  .f = function(df) {

                                    df$mean_tpm_pct_total <- (df$mean_tpm / sum(df$mean_tpm)) * 100

                                    df

                                  }
  )

  quantiles <- purrr::map(.x = output$group_data,
                          .f = ~ quantile(x = .x$mean_tpm[.x$mean_tpm > 0],
                                          probs = seq(0, 1, 0.25)))

  output$group_data <- purrr::map2(.x = output$group_data,
                                   .y = quantiles,
                                   .f = function(df, q) {

                                     df$mean_tpm_quartile <- dplyr::case_when(
                                       df$mean_tpm > q[["75%"]] ~ "Q4",
                                       df$mean_tpm > q[["50%"]] ~ "Q3",
                                       df$mean_tpm > q[["25%"]] ~ "Q2",
                                       df$mean_tpm > q[["0%"]] ~ "Q1",
                                       T ~ "not detected"
                                     )

                                     df

                                   }
  )

  if(!is.null(ref_gene_ids)) {

    for(i in seq_along(ref_gene_ids)) {

      nm <- paste0("mean_tpm_pct_", names(ref_gene_ids)[[i]])

      output$group_data <- purrr::map(.x = output$group_data,
                                      .f = function(df) {

                                        val <- df$mean_tpm / df[df[[gene_id]] == ref_gene_ids[[i]],
                                                                "mean_tpm"][[1]] * 100

                                        df[[nm]] <- val

                                        df

                                      }
      )
    }
  }

  if(!is.null(gene_groups)) {

    for(i in seq_along(gene_groups)) {

      nms <- paste0(c("mean_tpm_pct_", "mean_tpm_quartile_"),
                    names(gene_groups)[[i]])

      gene_group <- stats::na.omit(gene_groups[[i]])

      output$group_data <- purrr::map(.x = output$group_data,
                                      .f = function(df) {

                                        df[[nms[[1]]]] <- (df$mean_tpm / sum(df[df[[gene_id]] %in% gene_group,
                                                                                "mean_tpm"][[1]])) * 100

                                        df
                                      }
      )

      quantiles <- purrr::map(.x = output$group_data,
                              .f = ~ quantile(x = .x[(.x$mean_tpm > 0 & .x[[gene_id]] %in% gene_group),
                                                     "mean_tpm"][[1]],
                                              probs = seq(0, 1, 0.25)))

      output$group_data <- purrr::map2(.x = output$group_data,
                                       .y = quantiles,
                                       .f = function(df, q) {

                                         df[[nms[[2]]]] <- dplyr::case_when(
                                           df$mean_tpm > q[["75%"]] ~ "Q4",
                                           df$mean_tpm > q[["50%"]] ~ "Q3",
                                           df$mean_tpm > q[["25%"]] ~ "Q2",
                                           df$mean_tpm > q[["0%"]] ~ "Q1",
                                           T ~ "not detected"
                                         )

                                         df
                                       }
      )
    }
  }

  output <- tidyr::unnest(data = output,
                          cols = tidyselect::all_of("group_data"))

  output <- tidyr::unnest(data = output,
                          cols = tidyselect::all_of("sample_data"))

  if(cnt_nm != "tpm") {

    output <- dplyr::rename_with(.data = output,
                                 .fn = function(x) gsub(pattern = "tpm", replacement = cnt_nm, x = x),
                                 .cols = c("tpm",
                                           "mean_tpm",
                                           "sd_tpm",
                                           "sem_tpm",
                                           grep(pattern = "^mean_tpm_pct.*$",
                                                x = colnames(output),
                                                value = TRUE),
                                           grep(pattern = "^mean_tpm_quartile.*$",
                                                x = colnames(output),
                                                value = TRUE)))

  }

  dplyr::select(.data = output,
                group_id_var, sample_id_var, gene_id, tidyselect::everything())

}
