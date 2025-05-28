#' Extract Bioanalyzer results
#'
#' @param file a Results.csv file produced by Agilent 2100 Expert software.
#'
#' @return a tibble with summary of bioanalyzer results.
#' @export
#'
#' @examples
#' # Extract Bioanalyzer results
#' bio_res <- extract_bioanalyzer_results(file = rnaseqtools_example(file = "bio_res_ex.csv"))
#'
extract_bioanalyzer_results <- function(file) {

  output <- utils::read.csv(file = file, header = F)

  output <- tibble::tibble(sample_id = output[output$V1 == "Sample Name", "V2"][1:length(output[output$V1 == "Sample Name", "V2"])-1],
                           concentration_pg.ul = output[output$V1 == "RNA Concentration:", "V2"][1:length(output[output$V1 == "RNA Concentration:", "V2"])-1],
                           rin = output[output$V1 == "RNA Integrity Number (RIN):", "V2"],
                           rrna_28s_18s_ratio = output[output$V1 == "rRNA Ratio [28s / 18s]:", "V2"],
                           anomaly_thres_adapted = output[which(output$V1 == "RNA Integrity Number (RIN):")+1, "V1"])

  output$concentration_pg.ul <- as.double(output$concentration_pg.ul)
  output$rin <- as.double(gsub(pattern = "\\s.*$", replacement = "", x = output$rin))
  output$rrna_28s_18s_ratio <- as.double(output$rrna_28s_18s_ratio)
  output$anomaly_thres_adapted <- ifelse(output$anomaly_thres_adapted == " Anomaly Threshold(s) manually adapted)", T, F)

  output

}
