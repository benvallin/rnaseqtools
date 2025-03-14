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

#' log2(TPM+1) matrix
#'
#' Example log2(TPM+1) matrix with column and row names representing cell barcodes and ensembl gene IDs, respectively.
#' The matrix has 1000 rows and 547 columns.
#'
"log2_tpm1p"

#' MAST results
#'
#' Example results table derived from MAST::zlm using the workflow described at \url{https://bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html}.
#' The table was augmented with columns for alternative gene IDs, adjusted p value, log2 fold change and gene ranking metric.
#'
#' @format A tibble with 1000 rows and 10 variables:
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

#' DESeq2 results
#'
#' Example results table produced by DESeq2::results. The table was augmented with columns for alternative gene IDs and gene ranking metric.
#'
#' @format A tibble with 1000 rows and 10 variables:
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

#' Sample metadata
#'
#' Example sample metadata table matching the log2(TPM+1) matrix.
#'
#' @format A tibble with 547 rows and 3 variables:
#' \describe{
#'   \item{barcode}{cell barcodes matching the column names of the log2(TPM+1) matrix}
#'   \item{donor_id}{donor ID (donor1, donor2 or donor3)}
#'   \item{treatment}{treatment status (untreated or treated)}
#' }
"sample_metadata"

#' Multi fgsea results
#'
#' Example results table produced by rnaseqtools::multi_fgsea.
#'
#' @format A tibble with 32 rows and 12 variables:
#' \describe{
#'   \item{collection}{collection name (e.g.: MSigDB_C2_CP:REACTOME)}
#'   \item{set_name}{gene set name (e.g.: REACTOME_PHOSPHOLIPID_METABOLISM)}
#'   \item{pval}{see ?fgsea::fgsea}
#'   \item{padj}{see ?fgsea::fgsea}
#'   \item{log2err}{see ?fgsea::fgsea}
#'   \item{ES}{see ?fgsea::fgsea}
#'   \item{NES}{see ?fgsea::fgsea}
#'   \item{size}{see ?fgsea::fgsea}
#'   \item{leadingEdge}{see ?fgsea::fgsea}
#'   \item{n_leadingEdge}{number of genes in leading edge)}
#'   \item{n_leadingEdge_over_total}{number of genes in leading edge over number of genes in gene set}
#'   \item{pct_leadingEdge}{percentage of gene set being leading edge}
#' }
"multi_fgsea_results"

#' Multi fora results
#'
#' Example results table produced by rnaseqtools::multi_fora.
#'
#' @format A tibble with 49 rows and 9 variables:
#' \describe{
#'   \item{collection}{collection name (e.g.: MSigDB_H)}
#'   \item{set_name}{gene set name (e.g.: HALLMARK_OXIDATIVE_PHOSPHORYLATION)}
#'   \item{pval}{see ?fgsea::fora}
#'   \item{padj}{see ?fgsea::fora}
#'   \item{overlap}{see ?fgsea::fora}
#'   \item{size}{see ?fgsea::fora}
#'   \item{overlapGenes}{vector with overlapping genes}
#'   \item{n_overlapGenes}{number of overlapping genes over number of genes in gene set}
#'   \item{pct_overlapGenes}{percentage of gene set being overlapping genes}
#' }
"multi_fora_results"
