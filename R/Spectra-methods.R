#' @title Determine whether a spectrum contains peaks with certain m/z
#'
#' @rdname hasMz
#'
#' @aliases hasMz,Spectra-method hasMz
#'
#' @description
#'
#' `hasMz` allows to determine whether a spectrum contains peaks with their mz
#' matching any or all submitted `mz` values, accepting a certain deviation.
#'
#' @param object [Spectrum-class] or [Spectra()] object.
#'
#' @param mz `numeric` with the query m/z value(s) to be looked for in `object`.
#'
#' @param ppm `numeric(1)` defining the allowed deviation of the peak's m/z
#'     from `mz`. By default (`ppm = 10`), a difference of +/-
#'     `mz * 10/2e6` is allowed.
#'
#' @param which `character(1)` defining whether peaks with m/z matching any
#'     (`which = "any"`, default) of the provided `mz` are present, or all
#'     (`which = "all"`). See examples for details.
#'
#' @return `logical` of length equal `x` whether a peak with the given m/z is
#'     present in the spectrum.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @export
#'
#' @examples
#'
#' sp1 <- new("Spectrum1", mz = c(23.231, 123.43, 255.231, 432.0952),
#'     intensity = c(123, 3432, 45432, 423))
#' sp2 <- new("Spectrum1", mz = c(123.099, 344.531, 453.2313),
#'     intensity = c(231, 431, 413))
#' sp3 <- new("Spectrum1", mz = c(123.1001, 343.4321, 432.0921),
#'     intensity = c(542, 4524, 32))
#'
#' ## Test for a single Spectrum if it contains a peak with an m/z matching
#' ## any of the provided m/z.
#' mzs <- c(123.1, 432.0931)
#'
#' ## Does it contain any peak?
#' hasMz(sp1, mzs)
#'
#' ## Does it contain peaks matching all m/z?
#' hasMz(sp1, mzs, which = "all")
#'
#' ## Same for a list of Spectrum objects
#' spl <- Spectra(sp1, sp2, sp3)
#'
#' ## Do spectra contain any of the
#' hasMz(spl, mzs)
setMethod("hasMz", "Spectrum", function(object, mz, ppm = 10,
                                              which = c("any", "all")) {
    .spectrum_has_mz(object, mz = mz, ppm = ppm, which = which)
})
setMethod("hasMz", "Spectra", function(object, mz, ppm = 10,
                                             which = c("any", "all")) {
    vapply(object, FUN = .spectrum_has_mz, FUN.VALUE = logical(1),
           mz = mz, ppm = ppm, which = which)
})
