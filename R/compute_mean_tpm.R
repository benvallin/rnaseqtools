compute_mean_tpm <- function(input,
                             gene_id = "ensembl_gene_id_version",
                             sample_metadata,
                             sample_id_var,
                             group_id_var,
                             ref_gene_ids,
                             ref_pretty_gene_ids,
                             cnt_nm = "tpm") {

  output <- input

  if(is.matrix(input)) {

    output <- tibble::as_tibble(input,
                                rownames = gene_id)

  }

  output <- tidyr::pivot_longer(data = output,
                                cols = setdiff(colnames(output), gene_id),
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

  output$quantiles <- purrr::map(.x = output$group_data,
                                 .f = ~ quantile(x = .x$mean_tpm[.x$mean_tpm > 0],
                                                 probs = seq(0, 1, 0.25)))

  output$group_data <- purrr::map2(.x = output$group_data,
                                   .y = output$quantiles,
                                   .f = function(df, q) {

                                     df$mean_tpm_pct_total <- (df$mean_tpm / sum(df$mean_tpm)) * 100

                                     df <- dplyr::mutate(
                                       .data = df,
                                       mean_tpm_quartile = dplyr::case_when(
                                         .data$mean_tpm > q[["75%"]] ~ "Q4",
                                         .data$mean_tpm > q[["50%"]] ~ "Q3",
                                         .data$mean_tpm > q[["25%"]] ~ "Q2",
                                         .data$mean_tpm > q[["0%"]] ~ "Q1",
                                         T ~ "not detected"
                                       )
                                     )
                                     df
                                   }
  )

  # if(!is.null(ref_gene_ids) && !is.null(ref_pretty_gene_ids)) {
  #
  #   for(i in seq_along(ref_gene_ids)) {
  #
  #     nm <- paste0("mean_tpm_pct_", ref_pretty_gene_ids[[i]])
  #
  #     output$group_data <- purrr::map(.x = output$group_data,
  #                                     .f = function(df) {
  #
  #                                       val <- df$mean_tpm / df[df[[gene_id]] == ref_gene_ids[[i]],
  #                                                               "mean_tpm"][[1]] * 100
  #
  #                                       df[[nm]] <- val
  #
  #                                       df
  #
  #                                     }
  #     )
  #   }
  # }

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

  output <- output[, c(group_id_var, "group_data")]

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
                                 "mean_tpm_quartile",
                                 grep(pattern = "^mean_tpm_pct_.*$",
                                      x = colnames(output),
                                      value = TRUE)))

  }

  dplyr::select(.data = output,
                group_id_var, sample_id_var, gene_id, tidyselect::everything())

}
