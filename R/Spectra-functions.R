## Utility functions to work with Spectra objects.

#' Determine whether two numeric vectors are overlapping.
#'
#' @param x `numeric` with m/z values.
#'
#' @param mz `numeric` with m/z values to look for in `x`.
#'
#' @param ppm `numeric(1)` with the accepted difference between query and
#'     target m/z. A ppm of +/- `ppm` is applied to `mz` values assuming a
#'     measurement error of `ppm` on both `mz` and `x` values.
#'
#' @param tolerance `numeric(1)` with an acceptable absolute difference (i.e.
#'     differences between `x` and `mz` of +/- `tolerance are acceptable.
#'
#' @param which `character(1)` defining whether all input `mz` have to be
#'     present or only one.
#'
#' @return `logical(1)`.
#'
#' @author Johannes Rainer
#'
#' @noRd
#'
#' @importFrom MsCoreUtils ppm common
.has_mz <- function(x, mz, tolerance = 0, ppm = 10, which = c("any", "all")) {
    which <- match.arg(which)
    do.call(which, list(common(mz, x, tolerance = tolerance, ppm = ppm)))
}
