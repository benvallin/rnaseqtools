% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_leading_edge.R
\name{get_leading_edge}
\alias{get_leading_edge}
\title{Get leading edge}
\usage{
get_leading_edge(input, set_name)
}
\arguments{
\item{input}{data.frame or tibble with set_name column and leadingEdge list-column, as produced by rnaseqtools::multi_fgsea.}

\item{set_name}{character vector of length 1 representing gene set name.}
}
\value{
a character vector of gene IDs representing the leading edge for the queried gene set.
}
\description{
Get leading edge
}
\examples{
# Get genes in leading edge
leading_edge <- get_leading_edge(input = multi_fgsea_results_ex,
                                 set_name = "REACTOME_TRANSLATION")

}
