#' Gene metadata
#'
#' Example gene metadata derived from GENCODE annotation - Human release 46.
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
"gene_metadata_ex"

#' Bulk RNA-seq sample metadata
#'
#' Example sample metadata table matching the TPM matrix.
#'
#' @format A tibble with 6 rows and 2 variables:
#' \describe{
#'   \item{donor_id}{donor IDs matching the column names of the TPM matrix (donor1 to donor6)}
#'   \item{treatment}{disease status (healthy or diseased)}
#' }
"bulk_sample_metadata_ex"

#' Bulk-seq TPM matrix
#'
#' Example TPM matrix with column and row names representing donor IDs and ensembl gene IDs, respectively.
#' The matrix has 10073 rows and 6 columns.
#'
"bulk_tpm_ex"

#' scRNA-seq sample metadata
#'
#' Example sample metadata table matching the log2(TPM+1) matrix.
#'
#' @format A tibble with 547 rows and 3 variables:
#' \describe{
#'   \item{barcode}{cell barcodes matching the column names of the log2(TPM+1) matrix}
#'   \item{donor_id}{donor ID (donor1, donor2 or donor3)}
#'   \item{treatment}{treatment status (untreated or treated)}
#' }
"sc_sample_metadata_ex"

#' scRNA-seq log2(TPM+1) matrix
#'
#' Example log2(TPM+1) matrix with column and row names representing cell barcodes and ensembl gene IDs, respectively.
#' The matrix has 500 rows and 547 columns.
#'
"sc_log2_tpm1p_ex"

#' MAST summaryZlmFit
#'
#' Example summaryZlmFit object produced by MAST::summary(object = zlm, logFC = T, doLRT = T) using the workflow described at:
#' \url{https://bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html}.
#' The object contains a data.table with 1000 rows and 7 columns.
#'
"summaryZlmFit_ex"

#' MAST results
#'
#' Example MAST results table produced by rnaseqtools::tidy_summaryZlmFit.
#' The table was augmented with columns for alternative gene IDs.
#'
#' @format A tibble with 500 rows and 11 variables:
#' \describe{
#'   \item{ensembl_gene_id_version}{ensembl gene ID with version number (e.g.: ENSG00000145335.17)}
#'   \item{ensembl_gene_id}{ensembl gene ID (e.g.: ENSG00000145335)}
#'   \item{gene_symbol}{gene symbol (e.g.: SNCA)}
#'   \item{Pr(>Chisq)}{see ?MAST::`summary,ZlmFit-method`}
#'   \item{coef}{see ?MAST::`summary,ZlmFit-method`}
#'   \item{ci.hi}{see ?MAST::`summary,ZlmFit-method`}
#'   \item{ci.lo}{see ?MAST::`summary,ZlmFit-method`}
#'   \item{log2FoldChange}{log2 fold change computed as coef / log(x = 2, base = exp(1))}
#'   \item{lfcSE}{log2 fold change standard error computed from upper and lower bounds of confidence interval and level of confidence coefficient}
#'   \item{padj}{adjusted p value computed as p.adjust(p = `Pr(>Chisq)`, method = "BH")}
#'   \item{sign_log2fc_times_minus_log10pvalue}{gene ranking metric computed as sign(log2FoldChange) * -log(x = `Pr(>Chisq)`, base = 10)}
#' }
"mast_results_ex"

#' DESeq2 DESeqResults
#'
#' Example DESeqResults object produced by DESeq2::results.
#'
"DESeqResults_ex"

#' DESeq2 results
#'
#' Example DESeq2 results table produced by rnaseqtools::tidy_DESeqResults.
#' The table was augmented with columns for alternative gene IDs.
#'
#' @format A tibble with 500 rows and 10 variables:
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
"deseq2_results_ex"

#' Multi fgsea results
#'
#' Example results table produced by rnaseqtools::multi_fgsea.
#'
#' @format A tibble with 12 columns:
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
"multi_fgsea_results_ex"

#' Multi fora results
#'
#' Example results table produced by rnaseqtools::multi_fora.
#'
#' @format A tibble with 9 columns:
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
"multi_fora_results_ex"

#' Splici txi
#'
#' Example txi object produced by tximport::tximport and to supply to rnaseqtools::split_splici_txi.
#' It contains gene-level estimates of counts, abundances and lengths for spliced and unspliced transcripts.
#' To generate reference files for spliced and unspliced abundance estimation with alignment-free methods, see the worflow described at:
#' \url{https://bioconductor.org/packages/release/bioc/vignettes/eisaR/inst/doc/rna-velocity.html}.
#'
"txi_ex"

#' Split dataframe
#'
#' Example split_df to supply to rnaseqtools::split_splici_txi.
#' To generate reference files for spliced and unspliced abundance estimation with alignment-free methods, see the worflow described at:
#' \url{https://bioconductor.org/packages/release/bioc/vignettes/eisaR/inst/doc/rna-velocity.html}.
#'
#' @format A data.frame with 2 columns:
#' \describe{
#'   \item{spliced}{gene ID for spliced transcript (e.g.: ENSG00000005486.17 )}
#'   \item{unspliced}{gene ID for unspliced transcript (e.g.: ENSG00000005486.17-I)}
#' }
"split_df_ex"
