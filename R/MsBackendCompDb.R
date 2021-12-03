#' @title CompDb-based MS spectrum backend
#'
#' @aliases MsBackendCompDb-class
#'
#' @description
#'
#' The `MsBackendCompDb` allows to retrieve MS2 spectra from an [CompDb()]
#' object/database. The object keeps only a limited amount of data in memory
#' and retrieves the m/z and intensity values from the database *on-demand*. By
#' extending the [MsBackendCached()] class directly, `MsBackendCompDb` supports
#' adding/replacing spectra variables. These values are however only cached
#' within the object and not propagated (written) to the database.
#'
#' It is not intended that users create or use instances of this class directly,
#' the [Spectra()] call on [CompDb()] will return a `Spectra` object that uses
#' this backend.
#'
#' @param columns for `spectraData`: `character` with names of columns/spectra
#'     variables that should be returned. Defaults to
#'     `spectraVariables(object)`. Database columns `"ms_level"`,
#'     `"precursor_mz"`, `"precursor_intensity"`, `"precursor_charge"` are
#'     mapped to the core `Spectra` variables `msLevel`, `precursorMz`,
#'     `precursorIntensity` and `precursorCharge`, respectively.
#'
#' @param drop For `[`: not considered.
#'
#' @param filter for `backendInitialize`: optional filter expression to specify
#'     which elements to retrieve from the database.
#'
#' @param i For `[`: `integer`, `logical` or `character` to subset the object.
#'
#' @param j For `[`: not supported.
#'
#' @param name for `$<-`: the name of the spectra variable to replace.
#'
#' @param object an `MsBackendCompDb` instance.
#'
#' @param value for `$<-`: the replacement values.
#'
#' @param x an `MsBackendCompDb` instance.
#'
#' @param ... ignored.
#'
#' @section Methods implemented for `MsBackendCompDb`:
#'
#' The methods listed here are implemented for the `MsBackendCompDb`. All other
#' methods are inherited directly from the parent [MsBackendCached()] class.
#' See the help of [MsBackend()] in the `Spectra` package for a
#' complete listing of methods.
#'
#' - `peaksData`: gets the full list of peak matrices. Returns a [list()],
#'   length equal to the number of spectra and each element being a `matrix`
#'   with columns `"mz"` and `"intensity"` with the spectra's m/z and intensity
#'   values.
#'
#' - `intensity<-`: not supported.
#'
#' - `mz<-`: not supported.
#'
#' - `spectraData`: returns the complete spectrum data including m/z and
#'   intensity values as a [DataFrame()].
#'
#' - `$<-`: replace or add a spectrum variable. Note that `mz`, `intensity` and
#'   `spectrum_id` variables can not be replaced.
#'
#' @name MsBackendCompDb
#'
#' @author Johannes Rainer
#'
#' @exportClass MsBackendCompDb
#'
#' @importClassesFrom Spectra MsBackendCached
#'
#' @importClassesFrom DBI DBIConnection
#'
#' @importMethodsFrom Spectra peaksData
NULL

setClassUnion("DBIConnectionOrNULL", c("DBIConnection", "NULL"))

setClass("MsBackendCompDb",
         contains = "MsBackendCached",
         slots = c(dbcon = "DBIConnection",
                   spectraIds = "character",
                   .properties = "list"),
         prototype = prototype(
             dbcon = NULL,
             spectraIds = character(),
             .properties = list(),
             version = "0.1",
             readonly = TRUE))

setValidity("MsBackendCompDb", function(object) {
    msg <- .valid_dbcon(object@dbcon)
    if (length(msg)) msg
    else TRUE
})

#' @rdname MsBackendCompDb
#'
#' @importFrom methods callNextMethod
#'
#' @importFrom S4Vectors DataFrame
#'
#' @importMethodsFrom Spectra backendInitialize
#'
#' @exportMethod backendInitialize
setMethod("backendInitialize", "MsBackendCompDb", function(object,
                                                           x,
                                                           filter, ...) {
    if (missing(x))
        stop("Parameter 'x' is mandatory for 'MsBackendCompDb'")
    if (!is(x, "CompDb"))
        stop("Parameter 'x' has to be a 'CompDb' object")
    msg <- .valid_dbcon(.dbconn(x))
    if (length(msg))
        stop(msg)
    object@dbcon <- .dbconn(x)
    ## Get spectrum ID, precursor_mz and compound_id from db (msms_spectrum).
    ## Put that into localData.
    local_data <- .fetch_data(
        x, columns = c("spectrum_id", "precursor_mz", "compound_id"),
        filter = filter, start_from = "msms_spectrum")
    rownames(local_data) <- NULL
    colnames(local_data)[colnames(local_data) == "precursor_mz"] <-
        "precursorMz"

    object@spectraIds <- as.character(local_data$spectrum_id)
    ## Get info on tables and column names. Put them into spectraVariables.
    object@.properties$tables <- .tables(x)
    spectra_variables <- unique(unlist(.tables(object)))
    spectra_variables <- spectra_variables[!spectra_variables %in% c("peak_id")]
    object <- callNextMethod(
        object,
        data = local_data[, !colnames(local_data) == "spectrum_id"],
        nspectra = nrow(local_data),
        spectraVariables = .map_sql_to_spectraVariables(spectra_variables))
    validObject(object)
    object
})

#' @rdname MsBackendCompDb
#'
#' @exportMethod show
setMethod("show", "MsBackendCompDb", function(object) {
    callNextMethod()
    if (length(object)) {
        cat(" data source:", .metadata_value(object@dbcon, "source"), "\n")
        cat(" version:", .metadata_value(object@dbcon, "source_version"), "\n")
        cat(" organism:", .metadata_value(object@dbcon, "organism"), "\n")
    }
})

#' @importMethodsFrom Spectra peaksData
#'
#' @exportMethod peaksData
#'
#' @rdname MsBackendCompDb
setMethod("peaksData", "MsBackendCompDb", function(object) {
    .peaks_data(object)
})

#' @importMethodsFrom ProtGenerics dataStorage
#'
#' @exportMethod dataStorage
#'
#' @rdname MsBackendCompDb
setMethod("dataStorage", "MsBackendCompDb", function(object) {
    rep("<db>", length(object))
})

#' @exportMethod intensity<-
#'
#' @importMethodsFrom ProtGenerics intensity<-
#'
#' @rdname MsBackendCompDb
setReplaceMethod("intensity", "MsBackendCompDb", function(object, value) {
    stop("Can not replace original intensity values in the database.")
})

#' @exportMethod mz<-
#'
#' @importMethodsFrom ProtGenerics mz<-
#'
#' @rdname MsBackendCompDb
setReplaceMethod("mz", "MsBackendCompDb", function(object, value) {
    stop("Can not replace original data in the database.")
})

#' @importMethodsFrom ProtGenerics spectraData
#'
#' @exportMethod spectraData
#'
#' @rdname MsBackendCompDb
setMethod(
    "spectraData", "MsBackendCompDb",
    function(object, columns = spectraVariables(object)) {
        .spectra_data(object, columns = columns)
    })

#' @exportMethod spectraNames
#'
#' @importMethodsFrom ProtGenerics spectraNames
#'
#' @rdname MsBackendCompDb
setMethod("spectraNames", "MsBackendCompDb", function(object) {
    object@spectraIds
})

#' @exportMethod spectraNames<-
#'
#' @importMethodsFrom ProtGenerics spectraNames<-
#'
#' @rdname MsBackendCompDb
setReplaceMethod("spectraNames", "MsBackendCompDb",
                 function(object, value) {
                     stop(class(object)[1],
                          " does not support replacing spectra names (IDs).")
                 })

#' @exportMethod [
#'
#' @importFrom MsCoreUtils i2index
#'
#' @importFrom methods slot<-
#'
#' @importFrom S4Vectors extractROWS
#'
#' @rdname MsBackendCompDb
setMethod("[", "MsBackendCompDb", function(x, i, j, ..., drop = FALSE) {
    if (missing(i))
        return(x)
    i <- i2index(i, length(x), x@spectraIds)
    slot(x, "spectraIds", check = FALSE) <- x@spectraIds[i]
    x <- callNextMethod(x, i = i)
    x
})

#' @rdname MsBackendCompDb
#'
#' @export
setReplaceMethod("$", "MsBackendCompDb", function(x, name, value) {
    if (name %in% c("spectrum_id"))
        stop("Spectra IDs can not be changes.", call. = FALSE)
    callNextMethod()
})
