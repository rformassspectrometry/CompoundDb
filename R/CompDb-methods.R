#' @include CompDb.R

#' @description `dbconn` returns the connection (`DBIConnection`) to the
#'     database.
#'
#' @importMethodsFrom BiocGenerics dbconn
#'
#' @export
#'
#' @rdname CompDb
setMethod("dbconn", "CompDb", function(x) {
    .dbconn(x)
})

#' @importFrom methods show
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
        if (.has_msms_spectra(object)) {
            ## spctra <- dbGetQuery(con, paste0("select count(distinct spectrum_",
            ##                                  "id) from msms_spectrum_metadata"))
            spctra <- dbGetQuery(con, paste0("select count(distinct spectrum_",
                                             "id) from msms_spectrum"))
            cat(" MS/MS spectra count:", spctra[1, 1], "\n")
        }
    } else cat(" no database connection available\n")
})

# organism


#' @importMethodsFrom Spectra Spectra
#'
#'
#' @importClassesFrom Spectra Spectra
#'
#' @export
#'
#' @rdname CompDb
setMethod("Spectra", "CompDb", function(object, columns, filter, ...) {
    if (!.has_msms_spectra(object)) {
        warning("No spectrum data available in the provided database",
                call. = FALSE)
        return(Spectra())
    }
    sps <- new("Spectra")
    sps@backend <- backendInitialize(MsBackendCompDb(), x = object,
                                     columns = columns, filter = filter, ...)
    sps
})

#' @importMethodsFrom AnnotationFilter supportedFilters
#'
#' @export
#'
#' @rdname CompDb
setMethod("supportedFilters", "CompDb", function(object) {
    .supported_filters(object)
})
