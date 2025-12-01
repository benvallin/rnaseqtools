# Define paths ------------------------------------------------------------

# ==> To be set by user <==

# Path to directory containing source data
example_data_dir_path <- "/scratch/ben/rnaseq/tools/r_packages/example_data/rnaseqtools/"

# Path to output directory for raw data
raw_data_dir_path <- "inst/extdata/"

# Path to output directory for persistent user data
persistent_data_dir_path <- tools::R_user_dir(package = "rnaseqtools", which = "data")

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

# Define source data paths ------------------------------------------------

# Path to Salmon's meta_info file (json)
# salmon_meta_info_path <- "/scratch/ben/rnaseq/seq_data/ben/coculture_div67/gencode_v46/04_salmon_quant/Sample_2_S36_L003_CGTAGAACCTTGAGTT/aux_info/meta_info.json"
salmon_meta_info_path <- paste0(example_data_dir_path, "meta_info.json")

# Path to transcript metadata file (csv)
# transcript_metadata_path <- "/scratch/ben/rnaseq/ref_data/2024.08_reprocess/feature_metadata/transcript_metadata.csv"
transcript_metadata_path <- paste0(example_data_dir_path, "transcript_metadata.csv")

# Path to bulk RNA-seq sample metadata file (rds)
# bulk_sample_metadata_path <- "/scratch/ben/rnaseq/data/ben/pff_bulk/gencode_v46/output/make_sample_metadata/sample_metadata.rds"
bulk_sample_metadata_path <- paste0(example_data_dir_path, "bulk_sample_metadata.rds")

# Path to scRNA-seq sample metadata file (rds)
# sc_sample_metadata_path <- "/scratch/ben/rnaseq/data/ben/coculture_div67/output/cell_clustering/tidy_sample_metadata_filtered_cells/sample_metadata.rds"
sc_sample_metadata_path <- paste0(example_data_dir_path, "sc_sample_metadata.rds")

# Path to scRNA-seq log2(TPM+1) matrix file (csv)
# sc_log2_tpm1p_path <- "/scratch/ben/rnaseq/data/ben/coculture_div67/output/dge_gsea/dge_exposure_mast/mast_log2_tpm_lengthScaledTPM1p_neuron_B856_B156_B067_included_solo_recip_pc_genes_min_cnt_excl0_min_freq_incl0.2.csv"
sc_log2_tpm1p_path <- paste0(example_data_dir_path, "sc_log2_tpm1p.csv")

# Path to scRNA-seq raw count matrix file (csv)
# sc_cnt_path <- "/scratch/ben/rnaseq/data/ben/coculture_div67/output/count_processing/summarise_to_gene_level/raw_count.csv"
sc_cnt_path <- paste0(example_data_dir_path, "sc_cnt.csv")

# Path to MAST ZlmFit file (rds)
# mast_ZlmFit_path <- "/scratch/ben/rnaseq/data/ben/coculture_div67/output/dge_gsea/dge_exposure_mast/mast_glmer_zlm_log2_tpm_lengthScaledTPM1p_neuron_B856_B156_B067_included_solo_recip_pc_genes_min_cnt_excl0_min_freq_incl0.2_model_formula~_n_gene_on_+_exposure_+_(1_|_line_name)_padj0.05.rds"
mast_ZlmFit_path <- paste0(example_data_dir_path, "mast_ZlmFit.rds")

# Path to DESeq2 DESeqResults file (rds)
# deseq2_DESeqResults_path <- "/scratch/ben/rnaseq/data/ben/coculture_div67/output/dge_gsea/dge_exposure_zinbwave_deseq2/deseq2_lrt_res_raw_count_lengthScaledTPM_neuron_B856_B156_B067_included_solo_recip_pc_genes_min_cnt_excl0_min_freq_incl0.2_scran_sf_model_formula~_line_name_+_exposure_padj0.05.rds"
deseq2_DESeqResults_path <- paste0(example_data_dir_path, "deseq2_DESeqResults.rds")

# Path to directory containing Salmon quant files
salmon_quant_dir_path <- "/scratch/ben/rnaseq/seq_data/ben/coculture_div67/gencode_v46/05_salmon_splici_quant/"

# Path to splici metadata file (tsv)
# splici_metadata_path <- "/scratch/ben/rnaseq/ref_data/2024.08_reprocess/gencode.v46.tx2gene_expanded.tsv"
splici_metadata_path <- paste0(example_data_dir_path, "splici_metadata.tsv")

# Path to spliced - unspliced gene ID correspondence file (tsv)
# split_df_path <- "/scratch/ben/rnaseq/ref_data/2024.08_reprocess/gencode.v46.features_expanded.tsv"
split_df_path <- paste0(example_data_dir_path, "split_df.tsv")

# Path to Mari's list of target genes (csv)
# mcc_tars_path <- "/scratch/ben/rnaseq/data/mari/input/mari_target_genes2.csv"
mcc_tars_path <- paste0(example_data_dir_path, "mcc_tars.csv")

# Path to SNCA metadata (csv)
# snca_metadata_path <- "/scratch/ben/rnaseq/data/ben/pff_bulk/gencode_v46/output/make_feature_metadata/snca_metadata.csv"
snca_metadata_path <- paste0(example_data_dir_path, "snca_metadata.csv")

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

# Bulk RNA-seq data -------------------------------------------------------

# Sample metadata
bulk_sample_metadata <- read_rds(file = bulk_sample_metadata_path)

bulk_sample_metadata <- bulk_sample_metadata %>%
  filter(treatment == "NT", genotype %in% c("CTRL", "TRIP")) %>%
  mutate(donor_id = paste0("donor", seq(6)),
         disease_status = ifelse(genotype == "CTRL", "healthy", "diseased"),
         files = setNames(paste0("/scratch/ben/rnaseq/", files), donor_id))

bulk_sample_metadata_ex <- bulk_sample_metadata %>%
  dplyr::select(donor_id, disease_status)

use_data(bulk_sample_metadata_ex, overwrite = T)

# Import transcript metadata
transcript_metadata <- read_csv(transcript_metadata_path)

# Import counts from salmon output
txi <- tximport::tximport(files = bulk_sample_metadata$files,
                          type = "salmon",
                          txOut = FALSE,
                          countsFromAbundance = "no",
                          tx2gene = transcript_metadata)

# TPM matrix
bulk_tpm_ex <- txi$abundance

# Filter out non protein-coding genes
bulk_tpm_ex <- bulk_tpm_ex[rownames(bulk_tpm_ex) %in% gene_metadata_ex[gene_metadata_ex$gene_type == "protein_coding",
                                                                       "ensembl_gene_id_version"][[1]],]

# Filter out low TPM genes
bulk_tpm_ex <- bulk_tpm_ex[rowSums(bulk_tpm_ex >= 0) >= 3,]

# Rename samples
rownames(bulk_tpm_ex) <- str_replace(rownames(bulk_tpm_ex), "\\.\\d+$", "")

use_data(bulk_tpm_ex, compress = "xz", overwrite = T)

# scRNA-seq data ----------------------------------------------------------

# log2(TPM+1) matrix
sc_log2_tpm1p_ex <- data.table::fread(input = sc_log2_tpm1p_path) %>%
  left_join(gene_metadata_ex[, c("ensembl_gene_id_version", "ensembl_gene_id")],
            by = join_by(ensembl_gene_id_version)) %>%
  filter(!is.na(ensembl_gene_id)) %>%
  select(-ensembl_gene_id_version) %>%
  column_to_rownames(var = "ensembl_gene_id") %>%
  as.matrix()

keep <- names(head(x = sort(x = rowMeans(sc_log2_tpm1p_ex), decreasing = T), n = 500L))

sc_log2_tpm1p_ex <- sc_log2_tpm1p_ex[keep,]

use_data(sc_log2_tpm1p_ex, compress = "xz", overwrite = T)

# Raw count matrix
sc_cnt_ex <- data.table::fread(input = sc_cnt_path) %>%
  left_join(gene_metadata_ex[, c("ensembl_gene_id_version", "ensembl_gene_id")],
            by = join_by(ensembl_gene_id_version)) %>%
  filter(!is.na(ensembl_gene_id)) %>%
  select(-ensembl_gene_id_version) %>%
  column_to_rownames(var = "ensembl_gene_id") %>%
  as.matrix()

sc_cnt_ex <- sc_cnt_ex[rownames(sc_log2_tpm1p_ex), colnames(sc_log2_tpm1p_ex)]

use_data(sc_cnt_ex, compress = "xz", overwrite = T)

# Sample metadata
sc_sample_metadata_ex <- read_rds(file = sc_sample_metadata_path)

sc_sample_metadata_ex <- sc_sample_metadata_ex[match(x = colnames(sc_log2_tpm1p_ex),
                                                     table = sc_sample_metadata_ex$barcode),] %>%
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

use_data(sc_sample_metadata_ex, overwrite = T)

# # SingleCellAssay object
# sca_ex <- MAST::FromMatrix(exprsArray = sc_log2_tpm1p_ex,
#                            cData = as.data.frame(sc_sample_metadata_ex))
#
# use_data(sca_ex, compress = "xz", overwrite = T)

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
                                      collections = c("MSigDB_H", "MSigDB_C2_CP:REACTOME", "MSigDB_C2_CP:KEGG_LEGACY", "MSigDB_C5_GO:BP"),
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
                                    collections = c("MSigDB_H", "MSigDB_C2_CP:REACTOME", "MSigDB_C2_CP:KEGG_LEGACY", "MSigDB_C5_GO:BP"),
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
sc_sample_metadata <- read_rds(file = sc_sample_metadata_path)

sc_sample_metadata <- sc_sample_metadata[101:105, c("barcode", "barcode_spe_dirs")]

sc_sample_metadata <- sc_sample_metadata %>%
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
txi_ex <- tximport::tximport(files = sc_sample_metadata$files,
                             type = "salmon",
                             txOut = FALSE,
                             countsFromAbundance = "no",
                             tx2gene = splici_metadata)

use_data(txi_ex, overwrite = T)

# Calcium genes metadata --------------------------------------------------

calcium_genes_ex <- read_csv(mcc_tars_path)

calcium_genes_ex <- calcium_genes_ex %>%
  mutate(ensembl_gene_id = str_remove(ensembl_gene_id_version, "\\.\\d+$")) %>%
  select(ensembl_gene_id, gene_symbol = gene_name, protein_complex, functional_pathway)

use_data(calcium_genes_ex, overwrite = T)

# msa ---------------------------------------------------------------------

snca_metadata <- read_csv(snca_metadata_path)

snca_metadata <- snca_metadata %>%
  dplyr::filter(peptide_seq != "Sequence unavailable") %>%
  mutate(peptide_seq = str_remove(peptide_seq, "\\*")) %>%
  arrange(desc(width_peptide), desc(width_transcript))

asyn_metadata <- snca_metadata %>%
  mutate(isoform_id = paste0("alpha-synuclein - ", width_peptide, "aa")) %>%
  select(isoform_id, peptide_seq) %>%
  unique()

prot_seq <- lapply(X = asyn_metadata$peptide_seq,
                   FUN = function(x)  Biostrings::AAString(x = x))

prot_seq <- setNames(object = prot_seq,
                     nm = asyn_metadata$isoform_id)

prot_seq <- Biostrings::AAStringSet(x = prot_seq, use.names = TRUE)

msa_ex <- msa::msa(inputSeqs = prot_seq,
                   method = "ClustalW",
                   cluster = "upgma",
                   gapOpening = 7,
                   gapExtension = 0,
                   type = "protein",
                   order = "input")

use_data(msa_ex, overwrite = T)
