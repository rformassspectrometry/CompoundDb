#' @description
#'
#' Check if spectraData has all required columns.
#'
#' @noRd
#'
#' @param x spectraData `DataFrame`
.valid_spectra_data_required_columns <- function(x,
                                                 columns = c("dataStorage")) {
    if (nrow(x)) {
        missing_cn <- setdiff(columns, colnames(x))
        if (length(missing_cn))
            return(paste0("Required column(s): ",
                          paste(missing_cn, collapse = ", "),
                          " is/are missing"))
    }
    NULL
}

.valid_ms_backend_dbcon <- function(x) {
    tbls <- dbListTables(x)
    if (!any(tbls == "msms_spectrum"))
        return(paste0("Database has no MS/MS data available"))
    NULL
}

#' @rdname MsBackendCompDb
#'
#' @export
MsBackendCompDb <- function() {
    new("MsBackendCompDb")
}

#' Get m/z and/or intensity values.
#'
#' @param x `MsBackendCompDb`
#'
#' @param column `character` either empty, `"mz"` or `"intensity"` to get the
#'     data as a `SimpleList` of `matrix`, `NumericList` of m/z or `NumericList`
#'     of intensity values.
#'
#' @noRd
.peaks <- function(x, column = character()) {
    ids <- x@spectraData$spectrum_id
    pks <- dbGetQuery(
        .dbconn(x), paste0("select * from msms_spectrum_peak where spectrum_id",
                           " in (", paste0(ids, collapse = ","), ") order by ",
                           "peak_id"))
    rownames(pks) <- NULL
    if (length(column)) {
        if (column == "mz")
            split(pks$mz, factor(pks$spectrum_id, levels = ids))
        else split(pks$intensity, factor(pks$spectrum_id, levels = ids))
    } else
        split.data.frame(as.matrix(pks[, c("mz", "intensity")]),
                         factor(pks$spectrum_id, levels = ids))
}
