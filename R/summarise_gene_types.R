#' Summarise gene types
#'
#' @param input data.frame or tibble with character column <gene_id> and double columns representing sample-specific gene counts.
#' @param gene_id character vector of length 1 or name representing gene ID. Must be a column name of input.
#' @param sample_id character vector of length 1 or name representing sample ID.
#' @param gene_metadata data.frame or tibble with character column <gene_id>, gene_symbol and gene_type. Typically produced by rnaseqtools::make_gene_metadata.
#'
#' @return a data.frame summarising the number and percentage of protein coding, mitochondrial and ribosomal genes detected in each sample.
#' @export
#'
#' @examples
#' # Convert TPM matrix to tibble with gene IDs as column
#' tpm_df <- tibble::as_tibble(x = bulk_tpm_ex, rownames = "ensembl_gene_id")
#'
#' # Compute mean TPMs per disease status
#' gene_type_summary <- summarise_gene_types(input = tpm_df,
#'                                           gene_id = ensembl_gene_id,
#'                                           sample_id = donor_id,
#'                                           gene_metadata = gene_metadata_ex)
#'
summarise_gene_types <- function(input, gene_id, sample_id, gene_metadata) {

  gene_id <- substitute(gene_id)
  gene_id <- ifelse(is.symbol(gene_id), deparse(gene_id), eval(gene_id))

  sample_id <- substitute(sample_id)
  sample_id <- ifelse(is.symbol(sample_id), deparse(sample_id), eval(sample_id))

  if(!is.character(gene_id) ||
     length(gene_id) != 1L) {
    stop("Invalid gene_id argument.",
         call. = F)
  }

  if(!is.character(sample_id) ||
     length(sample_id) != 1L) {
    stop("Invalid sample_id argument.",
         call. = F)
  }

  if(!is.data.frame(input) ||
     !gene_id %in% colnames(input)) {
    stop("Input must be a data.frame or tibble with character column <gene_id> and numeric columns representing sample-specific gene counts.",
         call. = F)
  }

  if(!is.character(input[[gene_id]]) ||
     !all(vapply(X = setdiff(colnames(input), gene_id),
                 FUN = function(x) { is.numeric(input[[x]]) },
                 FUN.VALUE = NA))) {
    stop("Invalid column type in input.",
         call. = F)
  }

  if(!is.data.frame(gene_metadata) ||
     !all(vapply(X = c(gene_id, "gene_symbol", "gene_type"),
                 FUN = function(x) { x %in% colnames(gene_metadata) },
                 FUN.VALUE = NA))) {
    stop("gene_metadata must be a data.frame or tibble with character columns <gene_id>, gene_symbol and gene_type.",
         call. = F)
  }

  if(!all(vapply(X = c(gene_id, "gene_symbol", "gene_type"),
                 FUN = function(x) { is.character(gene_metadata[[x]]) },
                 FUN.VALUE = NA))) {
    stop("Invalid column type in gene_metadata.",
         call. = F)
  }

  input_ids <- input[[gene_id]]
  meta_ids <- gene_metadata[[gene_id]]

  if(length(input_ids) > length(stats::na.omit(meta_ids[match(x = input_ids, table = meta_ids)]))) {
    stop("Input contains gene_id value(s) not present in gene_metadata.",
         call. = F)
  }

  cnt <- as.data.frame(input)

  rownames(cnt) <- cnt[[gene_id]]

  cnt <- cnt[, -match(x = gene_id, table = colnames(cnt))]

  libsize <- colSums(cnt)

  df1 <- data.frame(libsize = libsize,
                    n_genes = colSums(cnt > 0))

  df1[[sample_id]] <- rownames(df1)

  pc_genes <- with(gene_metadata,
                   gene_metadata[gene_type == "protein_coding", gene_id][[1]])

  mt_genes <- with(gene_metadata,
                   gene_metadata[grepl(pattern = "^MT-", x = gene_symbol), gene_id][[1]])

  rb_genes <- with(gene_metadata,
                   gene_metadata[(grepl(pattern = "^RPS", x = gene_symbol) | grepl(pattern = "^RPL", x = gene_symbol)), gene_id][[1]])

  subcnts <- lapply(X = list(pc_genes, mt_genes, rb_genes),
                    FUN = function(x) {

                      cnt[rownames(cnt) %in% x,]

                    })

  names(subcnts) <- c("pc", "mt", "rb")


  df2 <- lapply(X = subcnts,
                FUN = function(x) { colSums(x > 0) })

  names(df2) <- paste0("n_", names(df2), "_genes")

  df2 <- as.data.frame(df2)

  df2[[sample_id]] <- rownames(df2)

  df3 <- lapply(X = subcnts,
                FUN = function(x) { colSums(x) / libsize * 100 })

  names(df3) <- paste0("pct_", names(df3), "_genes")

  df3 <- as.data.frame(df3)

  df3[[sample_id]] <- rownames(df3)

  Reduce(f = function(...) { merge(..., all = T) },
         x = list(df1, df2, df3))

}
