# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to output directory for raw data
raw_data_dir_path <- "inst/extdata/"

# Path to output directory for persistent user data
persistent_data_dir_path <- tools::R_user_dir(package = "rnaseqtools", which = "data")

# Path to sample metadata file (rds)
sample_metadata_path <- "/scratch/ben/rnaseq/data/ben/coculture_div67/output/cell_clustering/tidy_sample_metadata_filtered_cells/sample_metadata.rds"

# Path to Salmon's meta_info file (json)
salmon_meta_info_path <- "/scratch/ben/rnaseq/seq_data/ben/coculture_div67/gencode_v46/04_salmon_quant/Sample_2_S36_L003_CGTAGAACCTTGAGTT/aux_info/meta_info.json"

# Path to log2 COUNT+1 matrix file (csv)
log2_tpm1p_path <- "/scratch/ben/rnaseq/data/ben/coculture_div67/output/dge_gsea/dge_exposure_mast/mast_log2_tpm_lengthScaledTPM1p_neuron_B856_B156_B067_included_solo_recip_pc_genes_min_cnt_excl0_min_freq_incl0.2.csv"

# Path to MAST ZlmFit file (rds)
mast_ZlmFit_path <- "/scratch/ben/rnaseq/data/ben/coculture_div67/output/dge_gsea/dge_exposure_mast/mast_glmer_zlm_log2_tpm_lengthScaledTPM1p_neuron_B856_B156_B067_included_solo_recip_pc_genes_min_cnt_excl0_min_freq_incl0.2_model_formula~_n_gene_on_+_exposure_+_(1_|_line_name)_padj0.05.rds"

# Path to DESeq2 DESeqResults file (rds)
deseq2_DESeqResults_path <- "/scratch/ben/rnaseq/data/ben/coculture_div67/output/dge_gsea/dge_exposure_zinbwave_deseq2/deseq2_lrt_res_raw_count_lengthScaledTPM_neuron_B856_B156_B067_included_solo_recip_pc_genes_min_cnt_excl0_min_freq_incl0.2_scran_sf_model_formula~_line_name_+_exposure_padj0.05.rds"

# Path to directory containing Salmon quant files
salmon_quant_dir_path <- "/scratch/ben/rnaseq/seq_data/ben/coculture_div67/gencode_v46/05_salmon_splici_quant/"

# Path to transcript metadata file (csv)
transcript_metadata_path <- "/scratch/ben/rnaseq/ref_data/2024.08_reprocess/feature_metadata/transcript_metadata.csv"

# Path to splici metadata file (tsv)
splici_metadata_path <- "/scratch/ben/rnaseq/ref_data/2024.08_reprocess/gencode.v46.tx2gene_expanded.tsv"

# Path to spliced - unspliced gene ID correspondence file (tsv)
split_df_path <- "/scratch/ben/rnaseq/ref_data/2024.08_reprocess/gencode.v46.features_expanded.tsv"

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

# Set n cores for parallel processing
options(mc.cores = 48)

# GENCODE -----------------------------------------------------------------

# Download GENCODE annotation file
download_gencode_annotation(species = "human",
                            release = 46,
                            out_dir_path = persistent_data_dir_path)

# Import GENCODE annotation from downloaded file
gencode_file_path <- paste0(persistent_data_dir_path, "gencode.v46.primary_assembly.annotation.gtf.gz")

gencode_annotation <- import_gencode_annotation(file = gencode_file_path)

# Make gene metadata from GENCODE annotation
gene_metadata_ex <- make_gene_metadata(input = gencode_annotation)

use_data(gene_metadata_ex, overwrite = T)

# Clean up
file.remove(gencode_file_path)
rm(gencode_file_path, gencode_annotation)

# Salmon ------------------------------------------------------------------

# Create example meta_info_ex.json file
file.copy(from = salmon_meta_info_path,
          to = paste0(raw_data_dir_path, "meta_info_ex.json"),
          overwrite = T)

# Counts ------------------------------------------------------------------

# log2(TPM+1) matrix
log2_tpm1p_ex <- data.table::fread(input = log2_tpm1p_path) %>%
  left_join(gene_metadata_ex[, c("ensembl_gene_id_version", "ensembl_gene_id")],
            by = join_by(ensembl_gene_id_version)) %>%
  filter(!is.na(ensembl_gene_id)) %>%
  select(-ensembl_gene_id_version) %>%
  column_to_rownames(var = "ensembl_gene_id") %>%
  as.matrix()

keep <- names(head(x = sort(x = rowMeans(log2_tpm1p_ex), decreasing = T), n = 500L))

log2_tpm1p_ex <- log2_tpm1p_ex[keep,]

# log2_tpm1p_ex <- Matrix::Matrix(log2_tpm1p_ex, sparse = TRUE)

use_data(log2_tpm1p_ex, compress = "xz", overwrite = T)

# Samples -----------------------------------------------------------------

# Sample metadata
sample_metadata_ex <- read_rds(file = sample_metadata_path)

sample_metadata_ex <- sample_metadata_ex[match(x = colnames(log2_tpm1p_ex),
                                               table = sample_metadata_ex$barcode),] %>%
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

use_data(sample_metadata_ex, overwrite = T)

# MAST --------------------------------------------------------------------

# MAST summaryZlmFit
zlm <- read_rds(mast_ZlmFit_path)

summaryZlmFit_ex <- MAST::summary(object = zlm, logFC = T,
                                  doLRT = "exposureTRIP-exposed",
                                  level = 0.95, parallel = T)

temp <- summaryZlmFit_ex$datatable

temp$primerid <- gsub(pattern = "\\..*$",
                      replacement = "",
                      x = temp$primerid)

temp <- temp[temp$primerid %in% keep,]

temp <- temp[temp$contrast == "exposureTRIP-exposed" & temp$component %in% c("H", "logFC"),
             c("primerid", "component", "contrast", "Pr(>Chisq)", "ci.hi", "ci.lo", "coef")]

temp$contrast <- "treatmenttreated"

temp$contrast <- as_factor(temp$contrast)

summaryZlmFit_ex$datatable <- temp

rm(zlm, temp)

use_data(summaryZlmFit_ex, overwrite = T)

# MAST DGE results
mast_results_ex <- tidy_summaryZlmFit(input = summaryZlmFit_ex,
                                      coef = "treatmenttreated",
                                      level = 0.95,
                                      gene_id = "ensembl_gene_id")

mast_results_ex <- mast_results_ex %>%
  left_join(gene_metadata_ex[, c("ensembl_gene_id", "ensembl_gene_id_version", "gene_symbol")],
            by = join_by("ensembl_gene_id")) %>%
  select(ensembl_gene_id_version, ensembl_gene_id, gene_symbol, everything())

mast_results_ex <- mast_results_ex[match(x = keep, table = mast_results_ex$ensembl_gene_id),]

use_data(mast_results_ex, overwrite = T)

# DESeq2 ------------------------------------------------------------------

# DESeq2 DESeqResults
DESeqResults_ex <- read_rds(file = deseq2_DESeqResults_path)

DESeqResults_ex@rownames <- gsub(pattern = "\\..*$",
                                 replacement = "",
                                 x = DESeqResults_ex@rownames)

DESeqResults_ex <- DESeqResults_ex[DESeqResults_ex@rownames %in% keep,]

DESeqResults_ex@elementMetadata@listData <- list()

use_data(DESeqResults_ex, overwrite = T)

# DESeq2 DGE results
deseq2_results_ex <- tidy_DESeqResults(input = DESeqResults_ex,
                                       gene_id = "ensembl_gene_id")

deseq2_results_ex <- deseq2_results_ex %>%
  left_join(gene_metadata_ex[, c("ensembl_gene_id", "ensembl_gene_id_version", "gene_symbol")],
            by = join_by("ensembl_gene_id")) %>%
  select(ensembl_gene_id_version, ensembl_gene_id, gene_symbol, everything())

deseq2_results_ex <- deseq2_results_ex[match(x = keep, table = deseq2_results_ex$ensembl_gene_id),]

use_data(deseq2_results_ex, overwrite = T)

# fgsea -------------------------------------------------------------------

# Multi fgsea results
multi_fgsea_results_ex <- multi_fgsea(input = get_msigdb_collections(),
                                      collections = c("MSigDB_H", "MSigDB_C2_CP:REACTOME", "MSigDB_C2_CP:KEGG", "MSigDB_C5_GO:BP"),
                                      stats = deseq2_results_ex %>%
                                        dplyr::arrange(dplyr::desc(sign_log2fc_times_minus_log10pvalue)) %>%
                                        dplyr::pull(sign_log2fc_times_minus_log10pvalue, ensembl_gene_id),
                                      gene_id = "ensembl_gene_id",
                                      padj_threshold = 0.05,
                                      nPermSimple = 10000)

use_data(multi_fgsea_results_ex, overwrite = T)

# fora --------------------------------------------------------------------

# Multi fora results
multi_fora_results_ex <- multi_fora(input = get_msigdb_collections(),
                                    collections = c("MSigDB_H", "MSigDB_C2_CP:REACTOME", "MSigDB_C2_CP:KEGG", "MSigDB_C5_GO:BP"),
                                    genes = mast_results_ex %>%
                                      dplyr::filter(log2FoldChange < 0,
                                                    padj < 0.05) %>%
                                      dplyr::pull(ensembl_gene_id),
                                    universe = mast_results_ex %>%
                                      dplyr::filter(!is.na(log2FoldChange)) %>%
                                      dplyr::pull(ensembl_gene_id),
                                    gene_id = "ensembl_gene_id",
                                    padj_threshold = 0.05)

use_data(multi_fora_results_ex, overwrite = T)

# splici ------------------------------------------------------------------

# Import sample metadata
sample_metadata <- read_rds(file = sample_metadata_path)

sample_metadata <- sample_metadata[101:105, c("barcode", "barcode_spe_dirs")]

sample_metadata <- sample_metadata %>%
  mutate(files = setNames(object = map_chr(.x = barcode_spe_dirs,
                                           .f = ~ list.files(path = paste0(salmon_quant_dir_path, .x),
                                                             pattern = "quant.sf",
                                                             full.names = TRUE)),
                          nm = barcode)) %>%
  select(files, everything())

# Import transcript metadata
transcript_metadata <- read_csv(transcript_metadata_path)

keep <- unique(transcript_metadata$ensembl_gene_id_version)[101:110]

transcript_metadata <- transcript_metadata[transcript_metadata$ensembl_gene_id_version %in% keep,]

# Import splici metadata
splici_metadata <- read_tsv(splici_metadata_path, col_names = c("ensembl_transcript_id_version", "ensembl_gene_id_version"))

subset_splici_metadata <- function(splici_metadata, gene_ids) {

  keep <- lapply(X = gene_ids,
                 FUN = function(x) {
                   str_detect(splici_metadata$ensembl_gene_id_version, x)
                 })

  keep <- as.data.frame(keep, col.names = gene_ids)

  keep <- rowSums(keep)

  keep <- ifelse(keep == 1, T, F)

  splici_metadata[keep,]

}

splici_metadata <- subset_splici_metadata(splici_metadata = splici_metadata,
                                          gene_ids = keep)

# Import spliced - unspliced gene ID correspondence
split_df_ex <- read.delim(split_df_path, header = T, as.is = T) %>%
  dplyr::rename(unspliced = intron)

split_df_ex <- split_df_ex[split_df_ex$spliced %in% keep,]

use_data(split_df_ex, overwrite = T)

# Construct txi
txi_ex <- tximport::tximport(files = sample_metadata$files,
                             type = "salmon",
                             txOut = FALSE,
                             countsFromAbundance = "no",
                             tx2gene = splici_metadata)

use_data(txi_ex, overwrite = T)
