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

## Function to calculate masses from m/z given adduct definition(s)

## Function to calculate m/z from mass given adduct definition(s)

## Adduct definitions.
