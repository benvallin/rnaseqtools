# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to output directory for raw data
raw_data_dir_path <- "inst/extdata/"

# Path to output directory for persistent user data
persistent_data_dir_path <- tools::R_user_dir(package = "rnaseqtools", which = "data")

# Path to Salmon's meta_info file (json)
salmon_meta_info_path <- "/scratch/ben/rnaseq/seq_data/ben/coculture_div67/gencode_v46/04_salmon_quant/Sample_2_S36_L003_CGTAGAACCTTGAGTT/aux_info/meta_info.json"

# Path to DESeq2 results file (rds)
deseq2_res_path <- "/scratch/ben/rnaseq/data/ben/coculture_div67/output/dge_gsea/dge_exposure_zinbwave_deseq2/deseq2_lrt_lfc_raw_count_lengthScaledTPM_neuron_B856_B156_B067_included_solo_recip_pc_genes_min_cnt_excl5_min_freq_incl0.2_scran_sf_model_formula~_line_name_+_exposure_padj0.05.rds"

# Path to MAST results file (rds)
mast_res_path <- "/scratch/ben/rnaseq/data/ben/coculture_div67/output/dge_gsea/dge_exposure_mast/mast_glmer_lfc_log2_tpm_lengthScaledTPM1p_neuron_B856_B156_B067_included_solo_recip_pc_genes_min_cnt_excl0_min_freq_incl0.2_model_formula~_n_gene_on_+_exposure_+_(1_|_line_name)_padj0.05.rds"

# Path to log2 COUNT+1 matrix file (csv)
log2_tpm1p_path <- "/scratch/ben/rnaseq/data/ben/coculture_div67/output/dge_gsea/dge_exposure_mast/mast_log2_tpm_lengthScaledTPM1p_neuron_B856_B156_B067_included_solo_recip_pc_genes_min_cnt_excl0_min_freq_incl0.2.csv"

# Path to sample metadata file (rds)
sample_metadata_path <- "/scratch/ben/rnaseq/data/ben/coculture_div67/output/cell_clustering/tidy_sample_metadata_filtered_cells/sample_metadata.rds"

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
gencode_file_path <- paste0(persistent_data_dir_path, "gencode.v46.primary_assembly.annotation.gtf.gz")

gencode_annotation <- import_gencode_annotation(file = gencode_file_path)

# Make gene metadata from GENCODE annotation
gene_metadata <- make_gene_metadata(input = gencode_annotation)

use_data(gene_metadata, overwrite = T)

# Clean up
file.remove(gencode_file_path)
rm(gencode_file_path, gencode_annotation)

# Salmon-related datasets -------------------------------------------------

# Create example meta_info.json file
file.copy(from = salmon_meta_info_path,
          to = paste0(raw_data_dir_path, "meta_info.json"),
          overwrite = T)

# DGE-related datasets ----------------------------------------------------

# log2(TPM+1) matrix
log2_tpm1p <- data.table::fread(input = log2_tpm1p_path) %>%
  left_join(gene_metadata[, c("ensembl_gene_id_version", "ensembl_gene_id")],
            by = join_by(ensembl_gene_id_version)) %>%
  filter(!is.na(ensembl_gene_id)) %>%
  select(-ensembl_gene_id_version) %>%
  column_to_rownames(var = "ensembl_gene_id") %>%
  as.matrix()

keep <- names(head(x = sort(x = rowMeans(log2_tpm1p), decreasing = T), n = 1000L))

log2_tpm1p <- log2_tpm1p[keep,]

use_data(log2_tpm1p, overwrite = T)

# MAST DGE results
mast_results <- read_rds(file = mast_res_path) %>%
  select(ensembl_gene_id_version, ensembl_gene_id, gene_symbol = gene_name, `Pr(>Chisq)`:sign_log2fc_times_minus_log10pvalue)

mast_results <- mast_results[match(x = keep, table = mast_results$ensembl_gene_id),]

use_data(mast_results, overwrite = T)

# DESeq2 DGE results
deseq2_results <- read_rds(file = deseq2_res_path) %>%
  select(ensembl_gene_id_version, ensembl_gene_id, gene_symbol = gene_name, baseMean:padj, sign_log2fc_times_minus_log10pvalue)

deseq2_results <- deseq2_results[match(x = keep, table = deseq2_results$ensembl_gene_id),]

use_data(deseq2_results, overwrite = T)

# Sample metadata
sample_metadata <- read_rds(file = sample_metadata_path)

sample_metadata <- sample_metadata[match(x = colnames(log2_tpm1p), table = sample_metadata$barcode),] %>%
  select(barcode, donor_id = line_name, treatment = exposure) %>%
  mutate(donor_id = case_when(donor_id == "B856" ~ "donor1",
                              donor_id == "B156" ~ "donor2",
                              donor_id == "B067" ~ "donor3",
                              T ~ NA_character_) %>%
           fct_relevel(c("donor1", "donor2", "donor3")),
         treatment = case_when(treatment == "TRIP-exposed" ~ "treated",
                               treatment == "unexposed" ~ "untreated",
                               T ~ NA_character_) %>%
           fct_relevel(c("untreated", "treated")))

use_data(sample_metadata, overwrite = T)

# GSEA-related datasets ---------------------------------------------------

# Multi fgsea results
multi_fgsea_results <- multi_fgsea(input = get_msigdb_collections(),
                             collections = c("MSigDB_H", "MSigDB_C2_CP:REACTOME", "MSigDB_C2_CP:KEGG", "MSigDB_C5_GO:BP"),
                             stats = deseq2_results %>%
                               dplyr::arrange(dplyr::desc(sign_log2fc_times_minus_log10pvalue)) %>%
                               dplyr::pull(sign_log2fc_times_minus_log10pvalue, ensembl_gene_id),
                             gene_id = "ensembl_gene_id",
                             padj_threshold = 0.05,
                             nPermSimple = 10000)

use_data(multi_fgsea_results, overwrite = T)

# Multi fora results
multi_fora_results <- multi_fora(input = get_msigdb_collections(),
                                 collections = c("MSigDB_H", "MSigDB_C2_CP:REACTOME", "MSigDB_C2_CP:KEGG", "MSigDB_C5_GO:BP"),
                                 genes = mast_results %>%
                                   dplyr::filter(log2fc < 0,
                                                 padj < 0.05) %>%
                                   dplyr::pull(ensembl_gene_id),
                                 universe = mast_results %>%
                                   dplyr::filter(!is.na(log2fc)) %>%
                                   dplyr::pull(ensembl_gene_id),
                                 gene_id = "ensembl_gene_id",
                                 padj_threshold = 0.05)

use_data(multi_fora_results, overwrite = T)
