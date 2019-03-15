## Utility functions related to mass/adduct calculations and similar.

.ppm <- function(x, ppm = 10) {
    ppm * x / 1e6
}

## function with mass as input and ppm to search for compound in mass.

#' @title Match values allowing a certain difference in ppm
#'
#' @description
#'
#' Match numeric values in `x` against values in `y` and return indices of
#' **all** matches of `x` in `y` allowing differences between the values
#' specified with parameter `ppm` (part per million).
#'
#' @param x `numeric` with input values to match in `y`. Supposed to be a
#'     shorter `numeric` than `y`.
#'
#' @param y `numeric` with values to match `x` against.
#'
#' @param ppm `numeric(1)` defining the allowed difference between values to be
#'     still considered a match. Differences smaller than +/- `ppm` are
#'     accepted.
#'
#' @return `list` with `integer` representing the index in `y` where each
#'     element in `x` matches (given the provided ppm).
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @export
#'
#' @examples
#'
#' yvals <- abs(rnorm(1000))
#' yvals[123] <- yvals[2] - yvals[2] * 10 * 1e-6
#' yvals[124] <- yvals[2] + yvals[2] * 10 * 1e-6
#' yvals[125] <- yvals[2] + yvals[2] * 12 * 1e-6
#' xvals <- yvals[c(2, 3, 3, 20, 21, 20)]
#' xvals[2] <- xvals[2] + (10 * xvals[2] / 1e6)
#' xvals[3] <- xvals[3] - (10 * xvals[3] / 1e6)
#' xvals[6] <- xvals[6] + (12 * xvals[6] / 1e6)
#'
#' ## Perfect matches:
#' matchWithPpm(xvals, yvals)
#'
#' ## Match allowing +/- 10ppm difference
#' matchWithPpm(xvals, yvals, ppm = 10)
#'
#' ## Match allowing +/- 20ppm difference
#' matchWithPpm(xvals, yvals, ppm = 20)
matchWithPpm <- function(x, y, ppm = 0) {
    lapply(x, function(z, ppm) {
        which(abs(z - y) <= (.ppm(z, ppm)) + 1e-9)
    }, ppm = force(ppm))
}

#' @title Conversion functions between mass and m/z
#'
#' @description
#'
#' `mass2mz` and `mz2mass` allow to convert between (monoisotopic) mass and m/z
#' of specified ion adducts an *vice versa*. The functions return a `list` with
#' length equal to `x`, each element being a `numeric` with the m/z or mass for
#' the specified adducts. See below for examples.
#'
#' `adducts` retrieves the definitions for the supported ion adducts. The
#' function returns a `data.frame`.
#'
#' @param x `numeric` with the masses (or m/z values) to be converted.
#'
#' @param adduct either a `data.frame` with required columns `"name"`, `"nmol"`
#'     (number of molecules), `"charge"` (the total charge of the molecule)
#'     and `"massdiff"` (total mass difference), such as returned by the
#'     `adducts` function, or a `character` with the names of the adducts
#'     (see `adducts` for supported adducts).
#'
#' @rdname mass2mz
#'
#' @md
#'
#' @author Johannes Rainer, Jan Stanstrup
#'
#' @examples
#'
#' masses <- c(75.032028409, 105.042595, 162.115698458, 180.063385)
#' names(masses) <- c("a", "b", "c", "d")
#'
#' ## Calculate mz for adducts [M+H]+ and [M+Na]+
#' mzs <- mass2mz(masses, adduct = c("[M+H]+", "[M+Na]+"))
mass2mz <- function(x, adduct = adducts()) {
    if (is.character(adduct))
        adduct <- adducts()[adduct, ]
    adduct$charge <- abs(adduct$charge)
    lapply(x, function(z) {
        res <- (adduct$nmol * z + adduct$massdiff) / adduct$charge
        names(res) <- adduct$name
        res
    })
}

##mz2mass


#' @param pattern For `adducts`: optional `character(1)` specifying a pattern to
#'     be used to retrieve selected adducts, e.g. containing a hydrogen.
#'
#' @param polarity For `adducts`: optional `numeric(1)` to retrieve only adducts
#'     with positive (`polarity = 1`) or negative (`polarity = -1`) polarity.
#'
#' @param set For `adducts`: `character(1)` defining sets of adducts.
#'
#' @param ... For `adducts`: additional parameters to be passed to the [grep()]
#'     function.
#' @md
#'
#' @rdname mass2mz
adducts <- function(pattern, polarity, set, ...) {
    adds <- CompoundDb:::ADDUCTS
    if (!missing(polarity)) {
        if (polarity < 0)
            adds <- adds[adds$charge < 0, ]
        else adds <- adds[adds$charge > 0, ]
    }
    if (!missing(pattern))
        adds <- adds[grep(pattern, adds$name, ...), ]
    adds[order(adds$name), ]
}
