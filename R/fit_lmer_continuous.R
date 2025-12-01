#' Linear mixed effects model for gene expression data
#'
#' @details For genes detected in every (or almost every) cell, MAST zlm(method = "glmer") does not compute fold change estimates because the discrete component of the two-part hurdle model cannot be fit.
#' fit_lmer_continuous() provides an alternative fold change estimation for those genes by fitting a linear mixed effects model for the continuous component only.
#'
#' @param sca SingleCellAssay object.
#' @param assay_name character vector of length 1 representing the assay to use; must be an assay name of <sca>.
#' @param genes character vector of gene IDs to process; must be rownames of <sca>.
#' @param gene_id character vector of length 1 representing gene ID. Used to name the gene ID column of output.
#' @param model_formula character vector of length 1 representing the RHS of the model formula. Must start with '~' and contain '|' for random effect.
#' @param coef_var character vector of length 1 representing the coefficient of interest. Must be a column name of colData(<sca>).
#' @param coef_test_lvl character vector of length 1 representing the test level for the coefficient of interest. Must be a value of <coef_var>.
#' @param cnt_scale character vector of length 1 representing count scale. Either 'log2' or 'log'.
#' @param parallel logical vector of length 1 indicating if parallelisation should be enabled. If TRUE, use n detected cores - 1. Not implemented for windows.
#'
#' @return a tibble with lmerTest::lmer results.
#'
#' @export
#'
#' @examples
#' # Construct SingleCellAssay object from sc_log2_tpm1p_ex and sc_sample_metadata_ex
#' sca <- MAST::FromMatrix(exprsArray = sc_log2_tpm1p_ex,
#'                         cData = as.data.frame(sc_sample_metadata_ex))
#'
#' # Subset MAST glmer test results to only the genes with missing estimates
#' mast_zlm_glmer_failed <- mast_results_ex[is.na(mast_results_ex$log2FoldChange),]
#'
#' # Fit mixed effect model for genes with missing estimates
#' results <- fit_lmer_continuous(sca = sca,
#'                                genes = mast_zlm_glmer_failed$ensembl_gene_id,
#'                                gene_id = "ensembl_gene_id",
#'                                model_formula = "~ treatment + (1 | donor_id)",
#'                                coef_var = "treatment",
#'                                coef_test_lvl = "treated",
#'                                parallel = FALSE)
#'
fit_lmer_continuous <- function(sca,
                                assay_name = "et",
                                genes,
                                gene_id,
                                model_formula,
                                coef_var,
                                coef_test_lvl,
                                cnt_scale = "log2",
                                parallel = FALSE) {

  if(!requireNamespace("SummarizedExperiment", quietly = TRUE)) {

    stop("Package \"SummarizedExperiment\" must be installed to use this function.",
         call. = F)

  }

  if(!requireNamespace("lmerTest", quietly = TRUE)) {

    stop("Package \"lmerTest\" must be installed to use this function.",
         call. = F)

  }

  if(!requireNamespace("broom.mixed", quietly = TRUE)) {

    stop("Package \"broom.mixed\" must be installed to use this function.",
         call. = F)

  }

  if(!inherits(sca, "SingleCellAssay")) {

    stop("sca must be a SingleCellAssay object.",
         call. = F)

  }

  if(!is.character(assay_name) ||
     length(assay_name) != 1L ||
     !assay_name %in% SummarizedExperiment::assayNames(sca)) {

    stop("Invalid assay_name argument.",
         call. = F)

  }

  if(!all(genes %in% rownames(sca))) {

    stop("Not all <genes> are rownames of <sca>.",
         call. = F)

  }

  if(!is.character(gene_id) ||
     length(gene_id) != 1L) {

    stop("Invalid gene_id argument.",
         call. = F)

  }

  if(!is.character(model_formula) ||
     length(model_formula) != 1L ||
     !grepl(pattern = "^~.*\\|.+$", x = model_formula)) {

    stop("Invalid model_formula argument.",
         call. = F)

  }

  if(!is.character(coef_var) ||
     length(coef_var) != 1L ||
     !coef_var %in% colnames(SummarizedExperiment::colData(sca))) {

    stop("Invalid coef_var argument or <coef_var> not a column name of colData(sca).",
         call. = F)

  }

  if(!is.character(coef_test_lvl) ||
     length(coef_test_lvl) != 1L ||
     !coef_test_lvl %in% SummarizedExperiment::colData(sca)[[coef_var]]) {

    stop("Invalid coef_test_lvl argument or <coef_test_lvl> not a value of colData(sca)[[<coef_var>]].",
         call. = F)

  }

  if(!is.character(cnt_scale) ||
     length(cnt_scale) != 1L ||
     !cnt_scale %in% c("log2", "log")) {

    stop("Invalid cnt_scale argument.",
         call. = F)

  }

  if(!is.logical(parallel) ||
     length(parallel) != 1L) {

    stop("Invalid parallel argument.",
         call. = F)

  }

  cnt <- SummarizedExperiment::assay(sca, assay_name)

  metadata <- tibble::as_tibble(SummarizedExperiment::colData(sca),
                                rownames = "barcode_id")

  fit_fun <- function(gene) {

    gene_df <- cnt[gene,]

    gene_df <- tibble::tibble(barcode_id = names(gene_df),
                              y = as.numeric(gene_df))

    gene_df <- gene_df %>%
      dplyr::left_join(metadata,
                       by = "barcode_id")

    fit <- lmerTest::lmer(formula = stats::formula(paste0("y ", model_formula)),
                          data = gene_df,
                          REML = FALSE)

    res <- broom.mixed::tidy(x = fit, effects = "fixed")

    res <- res[res$term == paste0(coef_var, coef_test_lvl),]

    tibble::tibble(gene_id = gene,
                   log2FoldChange = res$estimate,
                   lfcSE = res$std.error,
                   pvalue = res$p.value)

  }

  if(parallel) {

    windows <- .Platform$OS.type == "windows"

    if(windows) {

      message("Parallelisation not implemented for windows... Setting parallel = FALSE.")

      parallel <- FALSE

    }

  }

  if(parallel) {

    output <- parallel::mclapply(X = genes,
                                 FUN = fit_fun,
                                 mc.cores = parallel::detectCores() - 1)

  } else {

    output <- lapply(X = genes,
                     FUN = fit_fun)

  }

  output <- dplyr::bind_rows(output)

  if(cnt_scale == "log") {

    output$log2FoldChange <- output$log2FoldChange / log(x = 2, base = exp(1))

    output$lfcSE <- output$lfcSE / log(x = 2, base = exp(1))

  }

  output$method <- "lmer_continuous"

  colnames(output) <- c(gene_id, "log2FoldChange", "lfcSE", "pvalue", "method")

  output

}
