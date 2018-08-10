#' @include CompDb.R

#' @description `dbconn` returns the connection (`DBIConnection`) to the
#'     database.
#'
#' @importMethodsFrom BiocGenerics dbconn
#'
#' @export
#'
#' @md
#' 
#' @rdname CompDb
setMethod("dbconn", "CompDb", function(x) {
    .dbconn(x)
})

#' @importFrom methods show
#'
#' @md
#' 
#' @export
setMethod("show", "CompDb", function(object) {
    cat("class:", class(object), "\n")
    con <- .dbconn(object)
    if (!is.null(con)) {
        cat(" data source:", .metadata_value(con, "source"), "\n")
        cat(" version:", .metadata_value(con, "source_version"), "\n")
        cat(" organism:", .metadata_value(con, "organism"), "\n")
        cmp_nr <- dbGetQuery(con, paste0("select count(distinct compound_id) ",
                                         "from compound"))
        cat(" compound count:", cmp_nr[1, 1], "\n")
        if (hasSpectra(object)) {
            ## spctra <- dbGetQuery(con, paste0("select count(distinct spectrum_",
            ##                                  "id) from msms_spectrum_metadata"))
            spctra <- dbGetQuery(con, paste0("select count(distinct spectrum_",
                                             "id) from msms_spectrum"))
            cat(" MS/MS spectra count:", spctra[1, 1], "\n")
        }
    } else cat(" no database connection available\n")
})

# organism


#' @importFrom ProtGenerics spectra
#'
#' @export
#'
#' @md
#'
#' @rdname CompDb
setMethod("spectra", "CompDb", function(object, columns, filter,
                                        return.type = c("Spectrum2List",
                                                        "data.frame",
                                                        "tibble")) {
    if (!hasSpectra(object))
        stop("No spectrum data available in the provided database",
             call. = FALSE)
    return.type <- match.arg(return.type)
    if (missing(columns))
        columns <- .tables(object, "msms_spectrum")[[1]]
    ## columns <- unique(unlist(
    ##     .tables(object, c("msms_spectrum_peak", "msms_spectrum_metadata"))))
    ## ordr <- "msms_spectrum_peak.spectrum_id"
    ordr <- "msms_spectrum.spectrum_id"
    if (return.type == "Spectrum2List") {
        columns <- unique(c("mz", "intensity", "polarity", "collision_energy",
                            columns))
        ## ordr <- paste0(ordr, ", msms_spectrum_peak.mz")
    }
    if (!any(columns == "spectrum_id"))
        columns <- c("spectrum_id", columns)
    res <- .fetch_data(object, columns = columns, filter = filter,
                       start_from = "msms_spectrum", order = ordr)
    ## start_from = "msms_spectrum_metadata", order = ordr)
    if (return.type == "tibble")
        res <- as_tibble(res)
    if (return.type == "Spectrum2List")
        res <- Spectrum2List(res)
    res
})
