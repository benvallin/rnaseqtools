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
