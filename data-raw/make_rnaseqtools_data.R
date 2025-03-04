# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to output directory for raw data
raw_data_dir_path <- "inst/extdata/"

# Path to output directory for persistent user data
persistent_data_dir_path <- tools::R_user_dir(package = "rnaseqtools", which = "data")

# Path to Salmon's meta_info file (json)
salmon_meta_info_path <- "/scratch/ben/rnaseq/seq_data/ben/coculture_div67/gencode_v46/04_salmon_quant/Sample_2_S36_L003_CGTAGAACCTTGAGTT/aux_info/meta_info.json"

# Set up ------------------------------------------------------------------

# Append trailing "/" to directory paths if missing
for(x in ls(pattern = "^.*dir_path$")) {
  value <- ifelse(grepl(pattern = "^.*/$", x = eval(as.symbol(x))),
                  eval(as.symbol(x)),
                  paste0(eval(as.symbol(x)), "/"))
  assign(x = x, value = value)
}

# Create output directory for raw data if necessary
if(!dir.exists(raw_data_dir_path)) { dir.create(raw_data_dir_path, recursive = T) }

# Import required libraries
library(devtools)
library(tidyverse)

# Load all rnaseqtools functions
load_all()

# GENCODE-related datasets ------------------------------------------------

# Download GENCODE annotation file
download_gencode_annotation(species = "human",
                            release = 46,
                            out_dir_path = persistent_data_dir_path)

# Import GENCODE annotation from downloaded file
gencode_file_path <- paste0(persistent_data_dir_path, "/gencode.v46.primary_assembly.annotation.gtf.gz")

gencode_annotation <- import_gencode_annotation(file = gencode_file_path)

# Make gene metadata from GENCODE annotation
gene_metadata <- make_gene_metadata(input = gencode_annotation)

use_data(gene_metadata, overwrite = T)

# Clean up
file.remove(gencode_file_path)
rm(persistent_data_dir_path, gencode_file_path, gencode_annotation)

# Salmon-related datasets -------------------------------------------------

# Create example meta_info.json file
file.copy(from = salmon_meta_info_path,
          to = paste0(raw_data_dir_path, "meta_info.json"),
          overwrite = T)



