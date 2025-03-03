## code to prepare `DATASET` dataset goes here

library(devtools)

load_all()

# Download GENCODE annotation file
out_dir_path <- tools::R_user_dir(package = "rnaseqtools", which = "data")

download_gencode_annotation(species = "human",
                            release = 46,
                            out_dir_path = out_dir_path)

# Import GENCODE annotation from downloaded file
gencode_file_path <- paste0(out_dir_path, "/gencode.v46.primary_assembly.annotation.gtf.gz")

gencode_annotation <- import_gencode_annotation(file = gencode_file_path)

# Make gene metadata from GENCODE annotation
gene_metadata <- make_gene_metadata(input = gencode_annotation)

use_data(gene_metadata, overwrite = T)

# Clean up
file.remove(gencode_file_path)
rm(out_dir_path, gencode_file_path, gencode_annotation)





