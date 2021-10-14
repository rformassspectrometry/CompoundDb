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

#' @importFrom DBI dbAppendTable dbGetQuery
#'
#' @export
#'
#' @rdname CompDb
setMethod("insertSpectra", signature(object = "CompDb", spectra = "Spectra"), 
          function(object, spectra, ...) {
              new_sD <- spectraData(spectra)
              if (!"compound_id" %in% colnames(new_sD))
                  stop("'spectra' must contain the variable 'compound_id'")
              dbcon <- object@dbcon
              if (!all(new_sD$compound_id %in%
                       dbGetQuery(dbcon, 
                                  "select compound_id from ms_compound")[, 1]))
                  stop("All values of spectra variable 'compound_id' must be",
                       " in 'compound_id' column of the database 'ms_compound'",
                       " table")
              if (!is.null(dbcon)) {
                  colnames(new_sD) <- sub("msLevel", "ms_level",
                                          colnames(new_sD), fixed = TRUE)
                  colnames(new_sD) <- sub("precursorMz", "precursor_mz",
                                          colnames(new_sD), fixed = TRUE)
                  colnames(new_sD) <- sub("precursorIntensity", 
                                          "precursor_intensity",
                                          colnames(new_sD), fixed = TRUE)
                  colnames(new_sD) <- sub("precursorCharge", "precursor_charge",
                                          colnames(new_sD), fixed = TRUE)
                  colnames(new_sD) <- sub("collisionEnergy", "collision_energy",
                                          colnames(new_sD), fixed = TRUE)
                  col_msms_spectrum <- spectraVariables(object)
                  n <- nrow(new_sD)
                  new_msms_spectrum <- matrix(NA, n, length(col_msms_spectrum))
                  colnames(new_msms_spectrum) <- col_msms_spectrum
                
                  cols <- intersect(col_msms_spectrum, colnames(new_sD))
                  new_msms_spectrum[, cols] <- as.matrix(new_sD[, cols])
                  nsp <- dbGetQuery(dbcon,
                                    paste0("select count(distinct spectrum_id)",
                                           " from msms_spectrum"))[1, 1]
                  new_msms_spectrum[, "spectrum_id"] <- nsp + seq_len(n)
                  dbAppendTable(dbcon, "msms_spectrum",
                                as.data.frame(new_msms_spectrum))
                  
                  new_pD <- Spectra:::.peaksapply(spectra)
                  np <- dbGetQuery(dbcon,
                                   paste0("select count(distinct peak_id) ",
                                          "from msms_spectrum_peak"))[1, 1]
                  nrows <- sapply(new_pD, nrow)
                  new_msms_spectrum_peak <-
                      cbind(spectrum_id = nsp + rep(seq_len(length(new_pD)),
                                                    nrows),
                            do.call(rbind, new_pD),
                            peak_id = np + seq_len(sum(nrows)))
                  dbAppendTable(dbcon, "msms_spectrum_peak",
                                as.data.frame(new_msms_spectrum_peak))
                  invisible(object)
              } else stop("Database not initialized")
          })
