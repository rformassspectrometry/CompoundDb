#' @description
#'
#' Function to collapse `mz` and `intensity` values from a `data.frame` per
#' spectrum (identified by column `"spectrum_id"`. The result will be a
#' `data.frame` with one row per spectrum and columns `"mz"` and `"intensity"`
#' being of type `list`.
#'
#' @param x `data.frame` with spectrum data for **one** spectrum.
#'
#' @return single row `data.frame` with columns `"mz"` and `"intensity"` being
#'     of type `list`.
#'
#' @noRd
#'
#' @md
#'
#' @author Johannes Rainer
.collapse_spectrum_df <- function(x) {
    ## Simple solution is to split by required column $spectrum_id that has to
    ## uniquely identify values belonging to one spectrum. More complicated
    ## solution would be to build a unique identified from all columns except
    ## mz and intensity.
    x_new <- base::unique(x[, !colnames(x) %in% c("mz", "intensity")])
    idf <- factor(x$spectrum_id, levels = base::unique(x$spectrum_id))
    x_new$mz <- unname(base::split(x$mz, idf))
    x_new$intensity <- unname(base::split(x$intensity, idf))
    rownames(x_new) <- NULL
    x_new[, colnames(x)]
}

#' @description
#'
#' Function to expand a collapsed spectrum `data.frame` (such as collapsed with
#' the `.collapse_spectrum_df` function) or returned by the [msms_spectra_hmdb]
#' function.
#'
#' @param x `data.frame` with columns `"mz"` and `"intensity"` being of type
#'     `list`.
#'
#' @noRd
#'
#' @md
#'
#' @author Johannes Rainer
.expand_spectrum_df <- function(x) {
    x_exp <- x[rep(seq_len(nrow(x)), lengths(x$mz)),
               !colnames(x) %in% c("mz", "intensity")]
    rownames(x_exp) <- NULL
    x_exp$mz <- unlist(x$mz)
    x_exp$intensity <- unlist(x$intensity)
    x_exp[, colnames(x)]
}

#' @title Expand m/z and intensity values in a data.frame
#'
#' @description
#'
#' `expandMzIntensity` *expands* a `data.frame` with m/z and/or intensity values
#' stored as a `list` in columns `"mz"` and `"intensity"`. The resulting
#' `data.frame` has the m/z and intensity values stored as `numeric` in columns
#' `"mz"` and `"intensity"`, one value per row, with the content of other
#' columns repeated as many times as there are m/z and intensity values.
#'
#' @param x `data.frame` with *collapsed* m/z and intensity values in columns
#'     `"mz"` and `"intensity"`, such as returned by [msms_spectra_hmdb()] with
#'     parameter `collapsed = TRUE`, or by [spectra()] or [compounds()] calls.
#'
#' @return `data.frame` with `"mz"` and `"intensity"` columns *expanded*. See
#'     description for details.
#'
#' @md
#'
#' @export
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' ## Read a data.frame with collapsed columns mz and intensity columns
#' dr <- system.file("xml/", package = "CompoundDb")
#' msms_spctra <- msms_spectra_hmdb(dr)
#'
#' msms_spctra
#'
#' ## Columns mz and intensity are "collased"
#' msms_spctra$mz
#'
#' ## Expand the data.frame to get one row per m/z and intensity value
#' spctra_exp <- expandMzIntensity(msms_spctra)
#' spctra_exp
expandMzIntensity <- function(x) {
    if (!is.data.frame(x))
        stop("'x' is expected to be a data.frame")
    .expand_spectrum_df(x)
}

#' Use a pattern search to extract the value for a field from a string
#'
#' @param x `character` string
#'
#' @param field `character` defining how the field can be identified
#'
#' @param delimiter `character` defining how the field's end is defined.
#'
#' @author Johannes Rainer
#'
#' @return `character` with the field or `NA` if the field is not present.
#'
#' @noRd
#'
#' @examples
#'
#' strng <- "some=bl df;other=some nice thing;last=the last entry"
#' .extract_field_from_string(strng, "some=", delimiter = ";")
#' .extract_field_from_string(strng, "other=", delimiter = ";")
#' .extract_field_from_string(strng, "last=", delimiter = ";")
.extract_field_from_string <- function(x, field, delimiter = " __ ") {
    res <- sub(paste0(".*?", field, "(.*?)(", delimiter, ".*|$)"), "\\1", x)
    res[res == x] <- NA_character_
    res
}
