#' Gene metadata
#'
#' Gene metadata derived from GENCODE annotation - Human release 46.
#' It can be produced by successive calls to download_gencode_annotation, import_gencode_annotation and make_gene_metadata.
#'
#' @format A tibble with 63140 rows and 5 variables:
#' \describe{
#'   \item{ensembl_gene_id_version}{ensembl gene ID with version number (e.g.: ENSG00000145335.17)}
#'   \item{ensembl_gene_id}{ensembl gene ID (e.g.: ENSG00000145335)}
#'   \item{gene_symbol}{gene symbol (e.g.: SNCA)}
#'   \item{gene_type}{gene type (e.g.: protein_coding)}
#'   \item{chr_name}{chromosome name (e.g.: chr4))}
#' }
"gene_metadata"

#' DESeq2 results
#'
#' Example results table produced by DESeq2::results. The table was augmented with columns for alternative gene IDs and gene ranking metric.
#'
#' @format A tibble with 8524 rows and 10 variables:
#' \describe{
#'   \item{ensembl_gene_id_version}{ensembl gene ID with version number (e.g.: ENSG00000145335.17)}
#'   \item{ensembl_gene_id}{ensembl gene ID (e.g.: ENSG00000145335)}
#'   \item{gene_symbol}{gene symbol (e.g.: SNCA)}
#'   \item{baseMean}{see ?DESeq2::results}
#'   \item{log2FoldChange}{see ?DESeq2::results}
#'   \item{lfcSE}{see ?DESeq2::results}
#'   \item{stat}{see ?DESeq2::results}
#'   \item{pvalue}{see ?DESeq2::results}
#'   \item{padj}{see ?DESeq2::results}
#'   \item{sign_log2fc_times_minus_log10pvalue}{gene ranking metric computed as sign(log2FoldChange) * -log(x = pvalue, base = 10)}
#' }
"deseq2_results"

#' MAST results
#'
#' Example results table derived from MAST::zlm using the workflow described at \url{https://bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html}.
#' The table was augmented with columns for alternative gene IDs, adjusted p value, log2 fold change and gene ranking metric.
#'
#' @format A tibble with 8524 rows and 10 variables:
#' \describe{
#'   \item{ensembl_gene_id_version}{ensembl gene ID with version number (e.g.: ENSG00000145335.17)}
#'   \item{ensembl_gene_id}{ensembl gene ID (e.g.: ENSG00000145335)}
#'   \item{gene_symbol}{gene symbol (e.g.: SNCA)}
#'   \item{Pr(>Chisq)}{see ?DESeq2::results}
#'   \item{coef}{see ?DESeq2::results}
#'   \item{ci.hi}{see ?DESeq2::results}
#'   \item{ci.lo}{see ?DESeq2::results}
#'   \item{padj}{adjusted p value computed as p.adjust(p = `Pr(>Chisq)`, method = "BH")}
#'   \item{log2fc}{log2 fold change computed as coef / log(x = 2, base = exp(1))}
#'   \item{sign_log2fc_times_minus_log10pvalue}{gene ranking metric computed as sign(log2fc) * -log(x = `Pr(>Chisq)`, base = 10)}
#' }
"mast_results"
