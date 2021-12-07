#' @title CompDb-based MS spectrum backend
#'
#' @description
#'
#' The `MsBackendCompDb` allows to retrieve MS2 spectra from an [CompDb()]
#' object/database. The object keeps only a limited amount of data in memory
#' and retrieves the m/z and intensity values from the database *on-demand*.
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
#' @param filter for `backendInitialize`: optional filter expression to specify
#'     which elements to retrieve from the database.
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
#' @return See the description of the respective function.
#'
#' @note
#'
#' For higher performance it is suggested to change the backend of the
#' [Spectra()] object to an [MsBackendDataFrame()] backend with the
#' [setBackend()] method of `Spectra` objects.
#'
#' @section Methods implemented for `MsBackendCompDb`:
#'
#' The methods listed here are implemented for the `MsBackendCompDb`. All other
#' methods are inherited directly from the parent [MsBackendDataFrame()] class.
#' See the help of `MsBackendDataFrame` in the `Spectra` package for a
#' complete listing of methods.
#'
#' - `peaksData`: gets the full list of peak matrices. Returns a [list()],
#'   length equal to the number of spectra and each element being a `matrix`
#'   with columns `"mz"` and `"intensity"` with the spectra's m/z and intensity
#'   values.
#'
#' - `intensity`: retrieves the intensity values for all spectra. Returns a
#'   [NumericList()], each element being the intensity values of one spectrum.
#'   The actual intensity values are retrieved from the database.
#'
#' - `intensity<-`: not supported.
#'
#' - `mz`: retrieves the m/z values for all spectra. Returns a [NumericList()],
#'   each element being a `numeric` with the m/z values of one spectrum. These
#'   values are retrieved from the database.
#'
#' - `mz<-`: not supported.
#'
#' - `spectraData`: returns the complete spectrum data including m/z and
#'   intensity values as a [DataFrame()].
#'
#' - `spectraData<-`: replace the spectrum metadata. Note that columns `"mz"`
#'   and `"intensity"` are ignored.
#'
#' - `$<-`: replace or add a spectrum variable. Note that `mz`, `intensity` and
#'   `spectrum_id` variables can not be replaced.
#'
#' @rdname MsBackendCompDb
#'
#' @author Johannes Rainer
#'
#' @exportClass MsBackendCompDb
#'
#' @importMethodsFrom Spectra peaksData
setClass("MsBackendCompDb",
         contains = "MsBackendDataFrame",
         slots = c("dbcon", "DBIConnection"),
         prototype = prototype(version = "0.1", readonly = TRUE))

setValidity("MsBackendCompDb", function(object) {
    msg <- .valid_spectra_data_required_columns(object@spectraData,
                                                c("spectrum_id"))
    msg <- c(msg, .valid_ms_backend_dbcon(object@dbcon))
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
setMethod("backendInitialize", "MsBackendCompDb", function(object,
                                                           x, columns,
                                                           filter, ...) {
    if (missing(x))
        stop("Parameter 'x' is mandatory for 'MsBackendCompDb'")
    if (!is(x, "CompDb"))
        stop("Parameter 'x' has to be a 'CompDb' object")
    if (missing(columns))
        columns <- .tables(x, "msms_spectrum")[[1]]
    columns <- columns[!(columns %in% c("mz", "intensity"))]
    msg <- .valid_ms_backend_dbcon(.dbconn(x))
    if (length(msg))
        stop(msg)
    ordr <- "msms_spectrum.spectrum_id"
    if (!any(columns == "spectrum_id"))
        columns <- c("spectrum_id", columns)
    spectraData <- .fetch_data(x, columns = columns, filter = filter,
                               start_from = "msms_spectrum", order = ordr)
    if (!nrow(spectraData)) {
        object@spectraData <- DataFrame(spectraData)
        return(object)
    }
    spectraData$dataStorage <- "<database>"
    spectraData$dataOrigin <- .metadata_value(x, "source")
    colnames(spectraData) <- sub("ms_level", "msLevel", colnames(spectraData),
                                 fixed = TRUE)
    colnames(spectraData) <- sub("precursor_mz", "precursorMz",
                                 colnames(spectraData), fixed = TRUE)
    colnames(spectraData) <- sub("precursor_intensity", "precursorIntensity",
                                 colnames(spectraData), fixed = TRUE)
    colnames(spectraData) <- sub("precursor_charge", "precursorCharge",
                                 colnames(spectraData), fixed = TRUE)
    if (any(colnames(spectraData) == "collision_energy") &&
        is.numeric(spectraData$collision_energy))
        colnames(spectraData) <- sub("collision_energy", "collisionEnergy",
                                     colnames(spectraData), fixed = TRUE)
    rownames(spectraData) <- spectraData$spectrum_id
    object@spectraData <- DataFrame(spectraData)
    object@dbcon <- .dbconn(x)
    validObject(object)
    object
})

#' @rdname MsBackendCompDb
#'
#' @export
setMethod("show", "MsBackendCompDb", function(object) {
    callNextMethod()
    if (length(object)) {
        cat(" data source:", .metadata_value(object@dbcon, "source"), "\n")
        cat(" version:", .metadata_value(object@dbcon, "source_version"), "\n")
        cat(" organism:", .metadata_value(object@dbcon, "organism"), "\n")
    }
})

#' @importFrom S4Vectors SimpleList
#'
#' @importMethodsFrom Spectra peaksData
#'
#' @rdname MsBackendCompDb
#'
#' @export
setMethod("peaksData", "MsBackendCompDb", function(object) {
    if (!length(object))
        return(list())
    .peaks(object)
})

#' @importFrom IRanges NumericList
#'
#' @importMethodsFrom ProtGenerics intensity
#'
#' @rdname MsBackendCompDb
#'
#' @export
setMethod("intensity", "MsBackendCompDb", function(object) {
    if (!length(object))
        return(NumericList())
    NumericList(.peaks(object, column = "intensity"), compress = FALSE)
})

#' @importMethodsFrom ProtGenerics intensity<-
#'
#' @rdname MsBackendCompDb
#'
#' @export
setReplaceMethod("intensity", "MsBackendCompDb", function(object, value) {
    stop(class(object), " does not support replacing intensity values")
})

#' @importMethodsFrom ProtGenerics mz
#'
#' @rdname MsBackendCompDb
#'
#' @export
setMethod("mz", "MsBackendCompDb", function(object) {
    if (!length(object))
        return(NumericList())
    NumericList(.peaks(object, column = "mz"), compress = FALSE)
})

#' @importMethodsFrom ProtGenerics mz<-
#'
#' @rdname MsBackendCompDb
#'
#' @export
setReplaceMethod("mz", "MsBackendCompDb", function(object, value) {
    stop(class(object), " does not support replacing m/z values")
})

#' @importMethodsFrom Spectra spectraData spectraVariables
#'
#' @rdname MsBackendCompDb
#'
#' @importFrom methods as
#'
#' @export
setMethod("spectraData", "MsBackendCompDb",
          function(object, columns = spectraVariables(object)) {
              have_cols <- intersect(columns, colnames(object@spectraData))
              res <- object@spectraData[, have_cols, drop = FALSE]
              miss_cols <- setdiff(columns, colnames(object@spectraData))
              if (any(miss_cols == "mz"))
                  res$mz <- mz(object)
              if (any(miss_cols == "intensity"))
                  res$intensity <- intensity(object)
              miss_cols <- miss_cols[!(miss_cols %in% c("mz", "intensity"))]
              if (length(miss_cols)) {
                  miss_res <- lapply(miss_cols, Spectra:::.get_column,
                                     x = object@spectraData)
                  names(miss_res) <- miss_cols
                  res <- cbind(res, as(miss_res, "DataFrame"))
              }
              res[, columns, drop = FALSE]
          })

#' @rdname MsBackendCompDb
#'
#' @importMethodsFrom Spectra spectraData<-
#'
#' @export
setReplaceMethod("spectraData", "MsBackendCompDb", function(object, value) {
    if (inherits(value, "DataFrame") &&
        any(colnames(value) %in% c("mz", "intensity"))) {
        warning("Ignoring columns \"mz\" and \"intensity\" ",
                "since 'MsBackendCompDb' does not support replacing them.")
        value <- value[, !(colnames(value) %in% c("mz", "intensity")),
                       drop = FALSE]
    }
    object@spectraData <- value
    validObject(object)
    object
})

#' @rdname MsBackendCompDb
#'
#' @export
setReplaceMethod("$", "MsBackendCompDb", function(x, name, value) {
    if (name == "mz" || name == "intensity" || name == "spectrum_id")
        stop("'MsBackendCompDb' does not support replacing mz, ",
             "intensity or spectrum_id values")
    value_len <- length(value)
    if (value_len == 1L || value_len == length(x))
        x@spectraData[[name]] <- value
    else
        stop("Length of 'value' has to be either 1 or ", length(x))
    validObject(x)
    x
})
