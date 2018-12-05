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

#' @title Collapsing and expanding data frames
#'
#' `collapse_table` collapses rows of a `data.frame` resulting in a
#' `data.frame` with the number of rows equal to unique elements in columns
#' `by`. Unique elements in the remaining columns will be combined into a
#' `list` per row (hence allowing for a different number of elements).
#'
#' @details
#'
#' `collapse_table` collapses the input `data.frame` based on columns `by`.
#' Reduction of the resulting `data.frame` is equivalent to
#' `unique(x[, by])`, with all other columns than `by` containing `list`s of
#' (unique) elements of these columns. See examples for details.
#'
#' @param by Columns defining how the table should be collapsed.
#'     The result `data.frame` will contain single, unique elements in these
#'     columns. See details section for more information.
#'     Can be a `character` with column names, an `integer` with column indices
#'     or a `logical` (same length then `ncol(x)`).
#'
#' @param column Columns that should
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
#'
#' @examples
#'
#' df <- data.frame(A = c("a", "a", "a", "c", "c", "d", "d", "d"),
#'     B = c("e", "e", "f", "f", "g", "g", "h", "h"),
#'     C = c("m", "n", "o", "p", "q", "r", "s", "t"),
#'     D = c(1, 2, 1, 2, 1, 1, 2, 2),
#'     stringsAsFactors = FALSE)
#'
#' ## Collapse by a single column.
#' res <- .collapse_table(df, by = "A")
#'
#' ## Column "A" contains unique elements
#' res$A
#'
#' ## Elements of all other columns have been collapsed into a list
#' ## of unique elements.
#' res$B
#' res$C
#'
#' res
#'
#' ## Collapse by the first two columns
#' res <- .collapse_table(df, by = c("A", "B"))
#'
#' ## Result for the first two columns is same as if unique was used
#' unique(df[, c("A", "B")])
#' res[, c("A", "B")]
#'
#' res$C
#'
#' res$D
.collapse_table <- function(x, by) {
    if (missing(by) || !length(by))
        return(x)
    by <- .column_indices(x, by)
    if (length(by) == 0)
        return(x)
    x <- as.list(x)
    ids <- do.call(paste, x[by])
    ids <- factor(ids, levels = base::unique(ids))
    res <- lapply(lapply(x, base::split, f = ids), function(z) {
        z_red <- lapply(z, base::unique)
        if (all(lengths(z_red) == 1))
            z_red <- unlist(z_red, use.names = FALSE)
        unname(z_red)
    })
    as.data.frame(as.tibble(res))
}

#' Simple helper function to return column indices in `x` with `y` being either
#' the column names, a `logical` or `integer` vector.
#'
#' @return `integer` with the index of the columns `y` in `x`
#'
#' @noRd
.column_indices <- function(x, y) {
    if (is.character(y))
        y <- match(y, colnames(x))
    if (is.logical(y))
        if (length(y) == ncol(x))
            y <- which(y)
        else stop("If a 'logical' is submitted its length has to match the ",
                  "number of columns of 'x'")
    y <- as.integer(y)
    if (is.integer(y))
        if (!all(y %in% seq_len(ncol(x))))
            stop("Index out of bounds")
    y[!is.na(y)]
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

## #' Convert inchi keys to an artificial compound ID
## #'
## #' @noRd
## .inchikey2id <- function(x, prefix = "CMP") {
##     x <- factor(x, levels = unique(x))
##     sprintf(paste0(prefix, "%0", ceiling(log10(length(levels(x)) + 1)), "d"),
##             as.integer(x))
## }

## #' Make a unique mona compounds table
## #'
## #' @noRd
## .reduce_mona_compounds <- function() {
##     ## HECK! don't have an inchikey for each!
##     cmp_id <- .inchikey2id(x$inchi)
##     tbl <- split.data.frame(x, f = cmp_id)
##     tbl_red <- lapply(tbl, function(z) {
##         if (length(unique(z$formula)) > 1)
##             stop("Formula is not unique")
##         data_frame(compound_id = NA_character_,
##                    compound_name = z$compound_name[1],
##                    inchi = z$inchi[1],
##                    formula = z$formula[1],
##                    mass = .aggregate_nums(z$mass),
##                    synonyms = unique(c(z$compound_name,
##                                        unlist(z$synonyms, use.names = FALSE))))
##     })
## }

## .aggregate_nums <- function(x, tolerance = 0.000001) {
##     if (length(x) == 1)
##         return(x)
##     if (all(diff(x) < tolerance))
##         x[1]
##     else stop("Difference of values is larger than allowed")
## }
