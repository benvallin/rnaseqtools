#' Set factor levels color
#'
#' @param input a data.frame or tibble.
#' @param variable character vector of length 1 representing the name of the variable to set as factor.
#' @param levels character vector with factor levels in desired order.
#' @param colors character vector with factor level colors in desired order.
#'
#' @return a tibble with variable set as factor and augmented with a variable_color column containing level colors.
#' @export
#'
#' @examples
#' # Set donor_id as factor and assign color levels
#' sample_metadata <- set_level_color(input = sc_sample_metadata_ex,
#'                                    variable = "donor_id",
#'                                    levels = c("donor3", "donor2", "donor1"),
#'                                    colors = c("#0077b6", "#a4133c", "#29a655"))
#'
set_level_color <- function(input, variable, levels, colors) {

  if(!is.data.frame(input) ||
     !variable %in% colnames(input)) {
    stop("Input must be a data.frame or tibble.",
         call. = F)
  }

  if(!is.character(variable) ||
     length(variable) != 1L ||
     !variable %in% colnames(input)) {
    stop("Invalid variable argument.",
         call. = F)
  }

  if(!is.character(levels) ||
     length(levels) != length(unique(levels)) ||
     length(levels) != length(unique(input[[variable]]))) {
    stop("Invalid levels argument.",
         call. = F)
  }

  if(!is.character(colors) ||
     length(colors) != length(levels)) {
    stop("Invalid colors argument.",
         call. = F)
  }

  if(paste0(variable, "_color") %in% colnames(input)) {
    stop("Column ", paste0(variable, "_color"), " already present in input.",
         call. = F)
  }

  output <- input

  fct_df <- tibble::tibble(factor = forcats::as_factor(x = levels),
                           color = forcats::as_factor(x = colors))

  fct_df <- fct_df[fct_df$factor %in% unique(output[[variable]]),]

  fct_df <- dplyr::mutate(.data = fct_df,
                          dplyr::across(.cols = tidyselect::everything(),
                                        .fns = forcats::fct_drop))

  colnames(fct_df) <- c(variable, paste0(variable, "_color"))

  output[[variable]] <- forcats::fct_relevel(output[[variable]],
                                             levels(fct_df[[variable]]))

  dplyr::left_join(x = output,
                   y = fct_df,
                   by = variable)

}
