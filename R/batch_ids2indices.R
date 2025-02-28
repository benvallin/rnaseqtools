# batch_ids2indices <- function(input,
#                               collections,
#                               identifiers) {
#
#   output <- lapply(X = collections,
#                    FUN = function(x) {
#
#                      temp <- tibble::as_tibble(input) %>%
#                        dplyr::select(ensembl_gene_id_version, tidyselect::all_of(x)) %>%
#                        unnest(all_of(x), keep_empty = F) %>%
#                        nest(data = -set_name)
#
#                      temp <- setNames(object = map(.x = temp$data,
#                                                    .f = ~ .x$ensembl_gene_id_version),
#                                       nm = temp$set_name)
#
#                      temp <- limma::ids2indices(gene.sets = temp,
#                                                 identifiers = identifiers,
#                                                 remove.empty = T)
#
#                    })
#
#   setNames(object = output, nm = collections)
#
# }
