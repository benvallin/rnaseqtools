#' Compute pseudobulk profiles
#'
#' @param input data.frame or tibble with character / factor column <cell_var> and numeric columns representing feature-specific counts.
#' @param cell_metadata data.frame or tibble with character / factor columns <cell_var> and <group_vars>.
#' @param cell_var character vector of length 1 or name representing cell ID. Must be a column name of input and cell_metadata.
#' @param group_vars character vector representing grouping variables for pseudobulking. Must be column name(s) of cell_metadata.
#' @param discard_not_expressed logical vector of length 1 indicating if features with 0 counts for all the pseudobulks should be discarded.
#' @param parallel logical vector of length 1 indicating if parallelisation should be enabled. If TRUE, use n detected cores - 1. Not implemented for windows.
#' @param skip_table_checks logical vector of length 1 indicating if format checking of input / cell_metadata should be skipped. If TRUE, malformed tables won't be detected but can greatly improve speed.
#'
#' @return a tibble with pseudobulk profiles defined by <group_vars>. Pseudobulks are produced by summing the feature counts for all cells in a group.
#' @export
#'
#' @examples
#' # Convert raw count matrix to tibble with gene IDs as column
#' sc_cnt_df <- tibble::as_tibble(x = sc_cnt_ex,
#'                                rownames = "ensembl_gene_id")
#'
#' # Transpose raw count tibble with barcode as column
#' sc_cnt_df <- transpose_cnt_df(input = sc_cnt_df,
#'                               gene_id = "ensembl_gene_id",
#'                               sample_id = "barcode")
#'
#' # Compute pseudobulk profiles for donor IDs
#' pseudobulks <- compute_pseudobulk(input = sc_cnt_df,
#'                                   cell_metadata = sc_sample_metadata_ex,
#'                                   cell_var = barcode,
#'                                   group_vars = "donor_id",
#'                                   parallel = FALSE)
#'
compute_pseudobulk <- function(input,
                               cell_metadata,
                               cell_var,
                               group_vars,
                               discard_not_expressed = TRUE,
                               parallel = TRUE,
                               skip_table_checks = FALSE) {

  cell_var <- substitute(cell_var)
  cell_var <- ifelse(is.symbol(cell_var), deparse(cell_var), eval(cell_var))

  if(!is.character(cell_var) ||
     length(cell_var) != 1L) {

    stop("Invalid cell_var argument.",
         call. = F)

  }

  if(!is.character(group_vars) ||
     (length(group_vars) == 1L && is.na(group_vars)) ||
     is.null(group_vars)) {

    stop("Invalid group_vars argument.",
         call. = F)

  }

  if(!is.logical(discard_not_expressed) ||
     length(discard_not_expressed) != 1L) {

    stop("Invalid discard_not_expressed argument.",
         call. = F)

  }

  if(!is.logical(parallel) ||
     length(parallel) != 1L) {

    stop("Invalid parallel argument.",
         call. = F)

  }

  if(!is.logical(skip_table_checks) ||
     length(skip_table_checks) != 1L) {

    stop("Invalid skip_table_checks argument.",
         call. = F)

  }

  if(!skip_table_checks) {

    if(!is.data.frame(input) ||
       !cell_var %in% colnames(input) ||
       (!is.character(input[[cell_var]]) && !is.factor(input[[cell_var]])) ||
       !all(vapply(X = setdiff(colnames(input), cell_var),
                   FUN = function(x) is.numeric(input[[x]]),
                   FUN.VALUE = NA))) {

      stop("Input must be a data.frame or tibble with only character / factor column <cell_var> and numeric columns representing feature-specific counts.",
           call. = F)

    }

    if(!is.data.frame(cell_metadata) ||
       !all(c(group_vars, cell_var) %in% colnames(cell_metadata)) ||
       (!is.character(cell_metadata[[cell_var]]) && !is.factor(cell_metadata[[cell_var]])) ||
       any(vapply(X = group_vars,
                  FUN = function(x) (!is.character(cell_metadata[[x]]) && !is.factor(cell_metadata[[x]])),
                  FUN.VALUE = NA))) {

      stop("cell_metadata must be a data.frame or tibble with character / factor columns <cell_var> and <group_vars>.",
           call. = F)

    }

    if(!all(unique(as.character(cell_metadata[[cell_var]])) %in% unique(as.character(input[[cell_var]])))) {

      stop("The <cell_var> values in cell_metadata must all be present in column <cell_var> of input.",
           call. = F)

    }

  }

  cnt <- as.data.frame(x = input)

  rownames(cnt) <- cnt[[cell_var]]

  cnt[[cell_var]] <- NULL

  pb_metadata <- unique(cell_metadata[, c(group_vars, cell_var)])

  pb_var_nm <- paste0(group_vars, collapse = "_")

  pb_var_val <- do.call(what = paste,
                        args = c(pb_metadata[, group_vars],
                                 sep = "_"))

  pb_metadata[[pb_var_nm]] <- pb_var_val

  cell_grps <- lapply(X = unique(pb_var_val),
                      FUN = function(x) {

                        pb_metadata[pb_metadata[[pb_var_nm]] == x,][[cell_var]]

                      })

  cell_grps <- stats::setNames(object = cell_grps,
                               nm = unique(pb_var_val))

  list(cell_grps = cell_grps,
       cnt = cnt)

  if(parallel) {

    windows <- .Platform$OS.type == "windows"

    if(windows) {

      message("Parallelisation not implemented for windows... Setting parallel = FALSE.")

      parallel <- FALSE

    }

  }

  if(parallel) {

    output <- parallel::mclapply(X = cell_grps,
                                 FUN = function(x) colSums(x = cnt[x, ]),
                                 mc.cores = parallel::detectCores() - 1)

  } else {

    output <- lapply(X = cell_grps,
                     FUN = function(x) colSums(x = cnt[x, ]))

  }

  output <- do.call(what = rbind,
                    args = output)

  if(discard_not_expressed) {

    output <- output[, colSums(output) > 0L]

  }

  tibble::as_tibble(x = output,
                    rownames = pb_var_nm)

}
