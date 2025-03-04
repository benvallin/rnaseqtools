#' Download GENCODE annotation
#'
#' @param species character vector of length 1 representing species. Must be "human" or "mouse".
#' @param release "latest" or numeric vector of length 1 representing release number.
#' @param out_dir_path output directory.
#'
#' @return download requested GENCODE annotation to output directory.
#' @export
#'
#' @examples
#' # Download GENCODE annotation file
#' persistent_data_dir_path <- tools::R_user_dir(package = "rnaseqtools", which = "data")
#'
#' download_gencode_annotation(species = "human",
#'                             release = 46,
#'                             out_dir_path = persistent_data_dir_path)
#'
#' # Clean up
#' gencode_file_path <- paste0(persistent_data_dir_path,
#'                             "/",
#'                             "gencode.v46.primary_assembly.annotation.gtf.gz")
#'
#' file.remove(gencode_file_path)
#'
#' rm(persistent_data_dir_path, gencode_file_path)
#'
download_gencode_annotation <- function(species = "human",
                                        release = "latest",
                                        out_dir_path) {

  if(!requireNamespace("RCurl", quietly = TRUE)) {

    stop("Package \"RCurl\" must be installed to use this function.",
         call. = F)

  }

  if(!species %in% c("human", "mouse")) {

    stop("Invalid species argument.",
         call. = F)

  }

  if(release != "latest") {

    if(!is.numeric(release) ||
       length(release) != 1L ||
       release < 1L) {

      stop("Invalid release argument.",
           call. = F)

    }

  }

  if(release == "latest") {

    url <- paste0("https://www.gencodegenes.org/", species, "/")

    release <- RCurl::getURL(url = url,
                             ftp.use.epsv = T,
                             dirlistonly = T,
                             .encoding = "UTF-8")

    release <- regmatches(x = release,
                          m = regexpr(pattern = ifelse(species == "human",
                                                       "Release \\d+",
                                                       "Release M\\d+"),
                                      text = release))

    release <- gsub(pattern = "Release ",
                    replacement = "",
                    x = release)

    message("GENCODE latest release for ", species, " is ", release, ".")

  }

  if(species == "mouse" &&
     !grepl(pattern = "^M\\d+$", x = release)) {

    release <- paste0("M", release)

  }

  message("Looking for:\nSpecies: ", species, "\nRelease: ", release)

  url_prefix <- paste0("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_", species,
                       "/release_", release, "/gencode.v", release, ".")

  url_suffix <- ".gtf.gz"

  annotations <- c("primary_assembly.annotation", "annotation")

  if(RCurl::url.exists(paste0(url_prefix, annotations[[1]], url_suffix))) {

    message("\"primary_assembly.annotation\" file is available.")

    annotation <- annotations[[1]]

  } else if(RCurl::url.exists(paste0(url_prefix, annotations[[2]], url_suffix))) {

    message("\"primary_assembly.annotation\" file is not available.\nDownloading \"annotation\" file instead.")

    annotation <- annotations[[2]]

  } else {

    stop("GENCODE release ", release, " is not available for ", species, ".",
         call. = F)

  }

  url <- paste0(url_prefix, annotation, url_suffix)

  if(!grepl(pattern = "^.*/$", x = out_dir_path)) {

    out_dir_path <- paste0(out_dir_path, "/")

  }

  if(!dir.exists(out_dir_path)) {

    dir.create(out_dir_path, recursive = T)

  }

  destfile <- paste0(out_dir_path, basename(url))

  utils::download.file(url = url,
                       destfile = destfile,
                       quiet = T)

  message("GENCODE annotation file ", basename(url), " added to ", out_dir_path, ".")

}

