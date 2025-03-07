% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_transcript_metadata.R
\name{make_transcript_metadata}
\alias{make_transcript_metadata}
\title{Make transcript metadata}
\usage{
make_transcript_metadata(input)
}
\arguments{
\item{input}{A GENCODE annotation data.frame or tibble as produced by import_gencode_annotation.}
}
\value{
A tibble with transcript metadata.
}
\description{
Make transcript metadata
}
\examples{
# Download GENCODE annotation file
out_dir_path <- tools::R_user_dir(package = "rnaseqtools", which = "data")

download_gencode_annotation(species = "human",
                            release = 46,
                            out_dir_path = out_dir_path)

# Import GENCODE annotation from downloaded file
gencode_file_path <- paste0(out_dir_path,
                            "/",
                            "gencode.v46.primary_assembly.annotation.gtf.gz")

gencode_annotation <- import_gencode_annotation(file = gencode_file_path)

# Make transcript metadata from GENCODE annotation
transcript_metadata <- make_transcript_metadata(input = gencode_annotation)

# Clean up
file.remove(gencode_file_path)

rm(out_dir_path, gencode_file_path, gencode_annotation)

}
