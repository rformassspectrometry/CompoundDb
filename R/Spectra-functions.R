## Utility functions to work with Spectra objects.

#' Determine whether a spectrum contains peaks with a certain m/z
#'
#' @param x `Spectrum` object
#'
#' @param mz `numeric` with m/z values to look for.
#'
#' @param ppm `numeric(1)` with the accepted difference between query and
#'     target m/z.
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
#' @md
.spectrum_has_mz <- function(x, mz, ppm = 10, which = c("any", "all")) {
    which <- match.arg(which)
    tiddle <- mz * ppm/2e6
    have_peaks <- outer(mz - tiddle, x@mz, "<=") &
        outer(mz + tiddle, x@mz, ">=")
    if (which == "any")
        any(have_peaks)
    else
        all(rowSums(have_peaks) > 0)
}
