% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_genes_in_collection.R
\name{get_genes_in_collection}
\alias{get_genes_in_collection}
\title{Get genes in collection}
\usage{
get_genes_in_collection(
  input,
  collection_name,
  set_name,
  gene_ids = c("ensembl_gene_id", "gene_symbol")
)
}
\arguments{
\item{input}{data.frame or tibble with gene_ids columns and collection list-column. Collection must contain column set_name.}

\item{collection_name}{character vector of length 1 representing gene collection name. Must be a column name of input.}

\item{set_name}{character vector of length 1 representing gene set name in collection.}

\item{gene_ids}{character vector of gene IDs. All elements must be column names of input.}
}
\value{
a tibble with gene IDs for queried set name in collection.
}
\description{
Get genes in collection
}
\examples{
# Make MSigDB collection table
msigdb_collection_table = get_msigdb_collections()

# Define collection of interest
collection_name <- "MSigDB_H"

# Define gene set name of interest
set_name <- "HALLMARK_OXIDATIVE_PHOSPHORYLATION"

# Get genes in gene set of interest
genes <- get_genes_in_collection(input = msigdb_collection_table,
                                 collection = collection_name,
                                 set_name = set_name,
                                 gene_ids = "ensembl_gene_id")

}
