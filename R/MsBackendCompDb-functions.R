#' @importFrom DBI dbListTables
#'
#' @noRd
.valid_dbcon <- function(x) {
    if (length(x)) {
        if (!inherits(x, "DBIConnection"))
            return("'dbcon' is expected to be a connection to a database")
        tables <- dbListTables(x)
        if (!all(c("msms_spectrum", "msms_spectrum_peak") %in% tables))
            return(paste0("Database has no MS/MS spectra data available."))
    }
    NULL
}

#' @rdname MsBackendCompDb
#'
#' @export
MsBackendCompDb <- function() {
    new("MsBackendCompDb")
}

.columns_sql <- c(
    precursorIntensity = "precursor_intensity",
    precursorMz = "precursor_mz",
    msLevel = "ms_level",
    compound_id = "compound_id",
    collisionEnergy = "collision_energy"
)

.map_spectraVariables_to_sql <- function(x) {
    for (i in seq_along(.columns_sql))
        x <- sub(names(.columns_sql)[i], .columns_sql[i], x, fixed = TRUE)
    x
}

.map_sql_to_spectraVariables <- function(x) {
    for (i in seq_along(.columns_sql))
        x <- sub(.columns_sql[i], names(.columns_sql[i]), x, fixed = TRUE)
    x
}

#' Get columns from the msms_spectrum_peak database table (dropping spectrum_id)
#'
#' @param x `MsBackendCompDb`
#'
#' @noRd
.available_peaks_variables <- function(x) {
    con <- .dbconn(x)
    if (length(con)) {
        if (length(.dbname(x)))
            on.exit(dbDisconnect(con))
        res <- dbGetQuery(con, "select * from msms_spectrum_peak limit 1")
        colnames(res)[!colnames(res) %in% c("spectrum_id", "peak_id")]
    } else character()
}

#' Returns a `data.frame` with the peaks data for spectra IDs in `x`. Note that
#' re-odering of the data needs to happen later.
#'
#' @param x `MsBackendCompDb`
#'
#' @noRd
.fetch_peaks <- function(x, columns = c("mz", "intensity")) {
    con <- .dbconn(x)
    if (length(con)) {
        if (length(.dbname(x)))
            on.exit(dbDisconnect(con))
        dbGetQuery(
            con,
            paste0("select spectrum_id,", paste(columns, collapse = ","),
                   " from msms_spectrum_peak where spectrum_id in (",
                   paste0("'", unique(x@spectraIds), "'", collapse = ","), ")"))
    } else {
        data.frame(spectrum_id = character(), mz = numeric(),
                   intensity = numeric())[, c("spectrum_id", columns)]
    }
}

#' Fetches the m/z and intensity values from the database and returns a list
#' of two column matrices (m/z, intensity). The function ensures that the data
#' is returned in the same order than x@spectraIds (also allowing duplicated
#' entries).
#'
#' @param x `MsBackendCompDb`.
#'
#' @author Johannes Rainer
#'
#' @noRd
.peaks_data <- function(x, columns = c("mz", "intensity")) {
    p <- .fetch_peaks(x, columns = columns)
    p <- unname(split.data.frame(p, as.factor(p$spectrum_id))[x@spectraIds])
    emat <- matrix(ncol = length(columns), nrow = 0,
                   dimnames = list(character(), columns))
    idx <- seq(2, (length(columns) + 1L))
    if (length(idx) == 1) {
        lapply(p, function(z) {
            if (nrow(z))
                matrix(z[, idx], dimnames = list(c(), columns))
            else emat
        })
    } else {
        lapply(p, function(z) {
            if (nrow(z))
                as.matrix(z[, idx], rownames.force = FALSE)
            else emat
        })
    }
}

#' @importFrom S4Vectors make_zero_col_DFrame extractCOLS
#'
#' @importFrom methods getMethod as
#'
#' @importFrom IRanges CharacterList NumericList
#'
#' @author Johannes Rainer
#'
#' @noRd
.spectra_data <- function(x, columns = spectraVariables(x)) {
    res <- getMethod("spectraData", "MsBackendCached")(x, columns = columns)
    if (is.null(res))
        res <- make_zero_col_DFrame(length(x))
    ## Define what needs to be still retrieved.
    db_cols <- intersect(columns, x@spectraVariables)
    db_cols <- db_cols[!db_cols %in% c("mz", "intensity", colnames(res))]
    peaks_cols <- intersect(columns, c("mz", "intensity"))

    if (length(db_cols)) {
        if (have_synonym <- any(db_cols == "synonym"))
            db_cols <- db_cols[!db_cols %in% c("synonym")]
        sp_data <- .fetch_data(
            x,
            columns = union("compound_id",
                            .map_spectraVariables_to_sql(db_cols)),
            filter = SpectrumIdFilter(unique(x@spectraIds)),
            start_from = "msms_spectrum")
        idx <- match(x@spectraIds, sp_data$spectrum_id)
        sp_data <- sp_data[idx[!is.na(idx)], , drop = FALSE]
        rownames(sp_data) <- NULL
        ## ? change data types for some variables ?
        if (any(colnames(sp_data) == "collision_energy"))
            sp_data$collision_energy <- as.numeric(sp_data$collision_energy)
        colnames(sp_data) <- .map_sql_to_spectraVariables(colnames(sp_data))
        res <- cbind(res, as(sp_data, "DataFrame"))
        if (have_synonym) {
            con <- .dbconn(x)
            tmp <- dbGetQuery(
                con,
                paste0("select * from synonym where compound_id in (",
                       paste0("'", unique(res$compound_id), "'",
                              collapse = ","), ")"))
            dbDisconnect(con)
            res$synonym <- CharacterList(
                unname(split(tmp$synonym, as.factor(tmp$compound_id))[
                    as.character(res$compound_id)]), compress = FALSE)
        }
    }
    if (length(peaks_cols)) {
        pks <- .fetch_peaks(x, columns = peaks_cols)
        if (any(peaks_cols == "mz"))
            res$mz <- NumericList(
                unname(split(pks$mz, as.factor(pks$spectrum_id))[x@spectraIds]),
                compress = FALSE)
        if (any(peaks_cols == "intensity"))
            res$intensity <- NumericList(
                unname(split(pks$intensity,
                             as.factor(pks$spectrum_id))[x@spectraIds]),
                compress = FALSE)
    }
    if (!all(columns %in% colnames(res)))
        stop("Column(s) ", paste0(columns[!columns %in% names(res)],
                                  collapse = ", "), " not available.",
             call. = FALSE)
    extractCOLS(res, columns)
}
