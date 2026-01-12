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
        if (length(.dbname(object)))
            on.exit(dbDisconnect(con))
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
                                        filter, ...) {
    if (!.has_msms_spectra(object)) {
        warning("No spectrum data available in the provided database",
                call. = FALSE)
        return(Spectra())
    }
    if (!.require_spectra())
        stop("The use of 'Spectra' requires package 'Spectra'. Please install ",
             "with 'BiocManager::install(\"Spectra\")'")
    sps <- new("Spectra")
    sps@backend <- backendInitialize(
        MsBackendCompDb(), x = object, filter = filter, ...)
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
#' @exportMethod compounds
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
          function(object, spectra, columns = spectraVariables(spectra), ...) {
              con <- .dbconn(object)
              if (is.null(con))
                  stop("Database not initialized")
              if (length(.dbname(object)))
                  on.exit(dbDisconnect(con))
              new_sD <- as.data.frame(spectraData(spectra, columns))
              if (!any(colnames(new_sD) == "compound_id"))
                  stop("Column 'compound_id' needs to be provided (as a ",
                       "spectra variable to insert to the database).")
              if (!all(new_sD$compound_id %in%
                       dbGetQuery(
                           con, "select compound_id from ms_compound")[, 1]))
                  stop("All values of spectra variable 'compound_id' must be",
                       " in 'compound_id' column of the database 'ms_compound'",
                       " table")
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
              new_sD$mz <- as.list(spectra$mz)
              new_sD$intensity <- as.list(spectra$intensity)
              if (hasMsMsSpectra(object))
                  .append_msms_spectra(con, new_sD)
              else {
                  if("spectrum_id" %in% colnames(new_sD))
                      warning("'spectrum_id' variable in 'spectra' will be
                              replaced with internal indexes")
                  new_sD$spectrum_id <- seq_len(nrow(new_sD))
                  .insert_msms_spectra(con, new_sD)
              }
              object@.properties$tables$msms_spectrum <- colnames(
                  dbGetQuery(con, "select * from msms_spectrum limit 1"))
              object@.properties$tables$msms_spectrum_peak <- colnames(
                  dbGetQuery(con, "select * from msms_spectrum_peak limit 1"))
              object
          })


#' @importFrom DBI dbGetQuery dbExecute
#'
#' @export
#'
#' @rdname CompDb
setMethod("deleteSpectra", signature(object = "CompDb"),
          function(object, ids = integer(0), ...) {
              dbcon <- .dbconn(object)
              if (is.null(dbcon))
                  stop("Database not initialized")
              if (length(.dbname(object)))
                  on.exit(dbDisconnect(dbcon))
              if (hasMsMsSpectra(object)) {
                  if(any(!ids %in%
                         dbGetQuery(dbcon, paste0("select spectrum_id ",
                                                  "from msms_spectrum"))[, 1]))
                      warning("Some IDs in 'ids' not valid and will be ignored")
                  dbExecute(dbcon, paste0("delete from msms_spectrum ",
                                          "where spectrum_id in (",
                                          toString(ids), ")"))
                  dbExecute(dbcon, paste0("delete from msms_spectrum_peak ",
                                          "where spectrum_id in (",
                                          toString(ids), ")"))
              } else {
                  stop("'object' does not contain msms spectra")
              }
              object
          })

#' @inherit MetaboCoreUtils::mass2mz
#'
#' @importFrom MetaboCoreUtils mass2mz
#'
#' @export
#'
#' @rdname CompDb
setMethod("mass2mz", signature = c("CompDb"),
          function(x, adduct = c("[M+H]+"), name = "formula") {
              cmps <- compounds(x, c(name, "exactmass"))
              res <- MetaboCoreUtils::mass2mz(cmps$exactmass, adduct)
              rownames(res) <- cmps[[name]]
              res
          })
setMethod("mass2mz", signature = c("ANY"),
          function(x, ...) {
              MetaboCoreUtils::mass2mz(x, ...)
          })

#' @importFrom MsCoreUtils rbindFill
#'
#' @rdname CompDb
#'
#' @export
setMethod("insertCompound", "CompDb", function(object, compounds = data.frame(),
                                               addColumns = FALSE) {
    if (!is.data.frame(compounds))
        stop("'compounds' is expected to be a data.frame")
    dbcon <- .dbconn(object)
    if (is.null(dbcon)) stop("Database not initialized")
    if (length(.dbname(object)))
        on.exit(dbDisconnect(dbcon))
    if (!nrow(compounds)) return(object)
    ref <- data.frame(name = character(), inchi = character(),
                      inchikey = character(), formula = character(),
                      exactmass = numeric(), synonyms = character())
    suppressWarnings(compounds <- rbindFill(compounds, ref))
    compounds$compound_id <- as.character(compounds$compound_id)
    .valid_compound(compounds, db = FALSE)
    if (is.list(compounds$synonyms) | !all(is.na(compounds$synonyms))) {
        if (is.list(compounds$synonyms)) {
            syn <- data.frame(
                compound_id = rep(compounds$compound_id,
                                  lengths(compounds$synonyms)),
                synonym = as.character(unlist(compounds$synonyms)))
        } else syn <- data.frame(compound_id = compounds$compound_id,
                                 synonym = compounds$synonyms)
        syn <- syn[!is.na(syn$synonym), ]
        if (nrow(syn)) dbAppendTable(dbcon, "synonym", syn)
    }
    dbcols <- colnames(dbGetQuery(dbcon, "select * from ms_compound limit 1"))
    compounds$synonyms <- NULL
    new_cols <- colnames(compounds)[!colnames(compounds) %in% dbcols]
    if (addColumns && length(new_cols)) {
        dtype <- dbDataType(dbcon, compounds[, new_cols, drop = FALSE])
        dtype <- paste(names(dtype), dtype)
        for (dt in dtype)
            dbExecute(dbcon, paste("alter table ms_compound add", dt))
        cols <- colnames(dbGetQuery(dbcon,
                                    "select * from ms_compound limit 1"))
        object@.properties$tables$ms_compound <- cols
        dbcols <- colnames(
            dbGetQuery(dbcon, "select * from ms_compound limit 1"))
    }
    dbAppendTable(dbcon, "ms_compound",
                  compounds[, intersect(dbcols, colnames(compounds))])
    object
})

#' @export
#'
#' @rdname CompDb
setMethod("deleteCompound", "CompDb", function(object, ids = character(),
                                               recursive = FALSE, ...) {
    dbcon <- .dbconn(object)
    if (is.null(dbcon)) stop("Database not initialized")
    if (length(.dbname(object)))
        on.exit(dbDisconnect(dbcon))
    if (!length(ids)) return(object)
    id_string <- paste0("'", ids, "'", collapse = ",")
    if (hasMsMsSpectra(object)) {
        ms <- dbGetQuery(dbcon, paste0("select compound_id, spectrum_id from ",
                                       "msms_spectrum where compound_id in (",
                                       id_string, ");"))
        if (nrow(ms)) {
            if (recursive)
                object <- deleteSpectra(object, ids = ms$spectrum_id)
            else stop("MS2 spectra for ", length(unique(ms$compound_id)),
                      " of the specified compounds present. Use parameter ",
                      "'recursive = TRUE' to delete compounds and all ",
                      "related MS2 spectra from the database.")

        }
    }
    dbExecute(dbcon, paste0("delete from synonym where compound_id in (",
                            id_string, ");"))
    dbExecute(dbcon, paste0("delete from ms_compound where compound_id in (",
                            id_string, ");"))
    object
})
