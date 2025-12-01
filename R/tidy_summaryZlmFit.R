#' Tidy MAST zlm summary
#'
#' @param input summaryZlmFit object produced by MAST::summary(object = zlm, logFC = T, doLRT = T) where zlm is a ZlmFit object. See ?MAST::`summary,ZlmFit-method`.
#' @param coef character vector of length 1 representing the coefficient of interest for LRT.
#' @param coef_scale character vector of length 1 representing coefficient scale. Either 'log2' or 'log'.
#' @param level level of confidence coefficient, should be the same as in the MAST::summary call.
#' @param gene_id character vector of length 1 representing gene ID. Used to rename the primerid column in output.
#'
#' @return a tibble with MAST zlm summary results.
#' @export
#'
#' @examples
#' # Extract MAST log fold change results from summaryZlmFit object
#' mast_results <- tidy_summaryZlmFit(input = summaryZlmFit_ex,
#'                                    coef = "treatmenttreated",
#'                                    coef_scale = "log2",
#'                                    level = 0.95,
#'                                    gene_id = "ensembl_gene_id")
#'
tidy_summaryZlmFit <- function(input, coef, coef_scale = "log2", level = 0.95, gene_id = "ensembl_gene_id") {

  if(!(methods::is(input, "summaryZlmFit") && methods::is(input, "list")) ||
     !all(c("primerid", "coef", "ci.hi", "ci.lo", "Pr(>Chisq)") %in% colnames(input$datatable))) {
    stop("Input must be a summaryZlmFit object produced by MAST::summary(object = zlm, logFC = T, doLRT = T).",
         call. = F)
  }

  if(!is.character(coef) ||
     length(coef) != 1L ||
     !coef %in% input$datatable$contrast) {
    stop("Invalid coef argument.",
         call. = F)
  }

  if(!is.character(coef_scale) ||
     length(coef_scale) != 1L ||
     !coef_scale %in% c("log2", "log")) {
    stop("Invalid coef_scale argument.",
         call. = F)
  }

  if(!is.numeric(level) ||
     length(level) != 1L ||
     level < 0 ||
     level > 1) {
    stop("Invalid level argument.",
         call. = F)
  }

  if(!is.character(gene_id) ||
     length(gene_id) != 1L) {
    stop("Invalid gene_id argument.",
         call. = F)
  }

  output <- tibble::as_tibble(input$datatable)

  output <- dplyr::full_join(x = output[output$contrast == coef & output$component == "H",
                                        c("primerid", "Pr(>Chisq)")],
                             y = output[output$contrast == coef & output$component == "logFC",
                                        c("primerid", "coef", "ci.hi", "ci.lo")],
                             by = "primerid")

  if(coef_scale == "log2") {

    output[["log2FoldChange"]] <- output[["coef"]]

    ci_high_log2 <- output[["ci.hi"]]

    ci_low_log2 <- output[["ci.lo"]]

  } else {

    output[["log2FoldChange"]] <- output[["coef"]] / log(x = 2, base = exp(1))

    ci_high_log2 <- output[["ci.hi"]] / log(x = 2, base = exp(1))

    ci_low_log2 <- output[["ci.lo"]] / log(x = 2, base = exp(1))

  }

  output[["lfcSE"]] <- (ci_high_log2 - ci_low_log2) / (stats::qnorm(1-((1-level)/2))*2)

  output[["pvalue"]] <- output[["Pr(>Chisq)"]]

  output[["padj"]] <- stats::p.adjust(p = output[["pvalue"]],
                                      method = "BH")

  output[["sign_log2fc_times_minus_log10pvalue"]] = sign(output[["log2FoldChange"]]) * -log(x = output[["pvalue"]],
                                                                                            base = 10)

  output[[gene_id]] <- output[["primerid"]]

  output[, c(gene_id, "coef", "ci.hi", "ci.lo",
             "log2FoldChange", "lfcSE", "pvalue", "padj", "sign_log2fc_times_minus_log10pvalue")] %>%
    dplyr::arrange(dplyr::desc(.data$sign_log2fc_times_minus_log10pvalue))

}
