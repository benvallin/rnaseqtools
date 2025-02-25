import_salmon_meta_info <- function(file) {

  meta_info <- rjson::fromJSON(file = file)

  meta_info <- meta_info[c("library_types",
                           "frag_dist_length", "frag_length_mean", "frag_length_sd",
                           "seq_bias_correct", "gc_bias_correct", "keep_duplicates",
                           "num_valid_targets", "num_decoy_targets",
                           "num_processed", "num_mapped",
                           "num_decoy_fragments", "num_dovetail_fragments",
                           "num_fragments_filtered_vm", "num_alignments_below_threshold_for_mapped_fragments_vm",
                           "percent_mapped")]

  meta_info <- as_tibble(meta_info) %>%
    mutate(file = file) %>%
    select(file, everything())

}
