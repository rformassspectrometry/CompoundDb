#' @include CompDb.R

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
                                         "from ms_compound"))
        cat(" compound count:", cmp_nr[1, 1], "\n")
        if (.has_msms_spectra(object)) {
            spctra <- dbGetQuery(con, paste0("select count(distinct spectrum_",
                                             "id) from msms_spectrum"))
            cat(" MS/MS spectra count:", spctra[1, 1], "\n")
        }
    } else cat(" no database connection available\n")
})

#' @importMethodsFrom Spectra Spectra
#'
#' @importClassesFrom Spectra Spectra
#'
#' @export
#'
#' @rdname CompDb
setMethod("Spectra", "CompDb", function(object,
                                        columns = spectraVariables(object),
                                        filter, ...) {
    if (!.has_msms_spectra(object)) {
        warning("No spectrum data available in the provided database",
                call. = FALSE)
        return(Spectra())
    }
    if (!requireNamespace("Spectra", quietly = TRUE))
        stop("The use of 'Spectra' requires package 'Spectra'. Please install ",
             "with 'Biobase::install(\"RforMassSpectrometry/Spectra\")'")
    sps <- new("Spectra")
    columns <- columns[!columns %in% c("mz", "intensity")]
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

#' @importMethodsFrom S4Vectors metadata
#'
#' @export
#'
#' @rdname CompDb
setMethod("metadata", "CompDb", function(x, ...) {
    .metadata(x)
})

#' @importMethodsFrom ProtGenerics spectraVariables
#'
#' @export
#'
#' @rdname CompDb
setMethod("spectraVariables", "CompDb", function(object, ...) {
    if (hasMsMsSpectra(object))
        .tables(object)$msms_spectrum
    else character()
})

#' @export
#'
#' @rdname CompDb
setMethod("compoundVariables", "CompDb", function(object,
                                                  includeId = FALSE, ...) {
    if (length(.tables(object))) {
        cn <- .tables(object)$ms_compound
        if (includeId)
            cn
        else cn[!cn %in% "compound_id"]
    } else character()
})

#' @importFrom tibble as_tibble
#'
#' @importMethodsFrom ProtGenerics compounds
#'
#' @export
#'
#' @rdname CompDb
setMethod("compounds", "CompDb", function(object,
                                          columns = compoundVariables(object),
                                          filter,
                                          return.type = c("data.frame",
                                                          "tibble"), ...) {
    return.type <- match.arg(return.type)
    if (length(columns))
        res <- .fetch_data(object, columns = columns, filter = filter,
                           start_from = "ms_compound")
    else res <- data.frame()
    if (return.type == "tibble")
        as_tibble(res)
    else res
})