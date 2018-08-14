#' @include Spectrum2List.R

#' @importFrom methods show
#'
#' @export
setMethod("show", "Spectrum2List", function(object) {
    .show_Spectrum2List(object, margin = "  ", print.classinfo = TRUE)
})


#' @rdname Spectrum2List
#'
#' @section Accessing metadata columns:
#'
#' Metadata columns can be accessed with `$`, that allows also to replace
#' columns or add new columns. The [DataFrame] with all metadata columns can
#' be accessed with the `mcols` function.
#' 
#' @export
#'
#' @param x `Spectrum2List`.
#'
#' @param object `Spectrum2List`.
#'
#' @param name for `$` and `$<-`: `character(1)` with the name of the
#'     metadata column to access.
#'
#' @param value for `$<-`: vector to replace/set the specified column in the
#'     object's metadata.
#' 
#' @examples
#'
#' ## Load Spectrum information from xml files
#' dr <- system.file("xml", package = "CompoundDb")
#'
#' spl <- Spectrum2List(msms_spectra_hmdb(dr))
#'
#' spl
#'
#' ## Access individual metadata columns using $
#' spl$spectrum_id
#'
#' ## Add additional metadata columns
#' spl$new_col <- 1:2
#'
#' spl
setMethod("$", "Spectrum2List", function(x, name) {
    eval(substitute(mcols(x)$NAME_ARG, list(NAME_ARG = name)))
})
#' @rdname Spectrum2List
#' 
#' @export
setReplaceMethod("$", "Spectrum2List", function(x, name, value) {
    mcols(x)[[name]] <- value
    x
})

#' @rdname Spectrum2List
#'
#' @section Accessing spectrum attributes:
#'
#' These methods allow to access the attributes and values of the individual
#' `Spectrum2` objects. For more details please refer to the corresponding
#' methods in the `MSnbase` package
#'
#' - `mz` return the m/z values of each spectrum as a `list` of `numeric`
#'   vectors.
#'
#' - `intensity` return the intensity values of each spectrum as a `list` of
#'   `numeric` vectors.
#'
#' - `rtime` return the retention time of each spectrum as a `numeric` vector
#'   with length equal to the length of `object`.
#'
#' - `precursorMz`, `precursorCharge`, `precursorIntensity`, `precScanNum`
#'   return precursor m/z values, charge, intensity and scan number for each
#'   spectrum as a `numeric` (or `integer`) vector with length equal to the
#'   length of `object`.
#'
#' - `acquisitionNum` and `scanIndex` return the acquisition number of each
#'   spectrum and its scan index as an `integer` vector with the same length
#'   than `object`.
#'
#' - `ionCount` and `tic` return the ion count and total ion current of each
#'   spectrum.
#'
#' - `peaksCount` returns the number of peaks for each spectrum as a `integer`
#'   vector.
#'
#' - `msLevel` returns the MS level of each spectrum.
#'
#' - `collisionEnergy` returns the collision energy for each spectrum.
#'
#' - `polarity` returns the spectra's polarity.
#'
#' - `fromFile` returns the index from the (e.g. mzML) file the spectra where
#'   from. This applies only for spectra read using e.g. `MSnbase`'s
#'   `readMSData` function.
#'
#' - `smoothed` whether spectra have been smoothed (i.e. processed with the
#'   [MSnbase::smooth] method. Returns a `logical` of length equal to the
#'   number of spectra.
#'
#' - `isEmpty` returns `TRUE` for spectra without peak data.
#'
#' - `centroided`, `isCentroided` returns for each spectrum whether it contains
#'   *centroided* data. While `centroided` returns the internal attribute of
#'   each spectrum, `isCentroided` tries to guess whether spectra are
#'   centroided from the actual peak data.
#' 
#' @importMethodsFrom MSnbase mz
#' 
#' @export
#'
#' @examples
#'
#' ## Extract the mz values for the individual spectra
#' mz(spl)
setMethod("mz", "Spectrum2List", function(object) {
    lapply(object, function(z) z@mz)
})

#' @rdname Spectrum2List
#'
#' @importMethodsFrom MSnbase intensity
#' 
#' @export
#'
#' @examples
#'
#' ## Extract the intensity values for the individual spectra
#' intensity(spl)
setMethod("intensity", "Spectrum2List", function(object) {
    lapply(object, function(z) z@intensity)
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase rtime
#' 
#' @export
#'
#' @examples
#'
#' ## Extract the retention time values for the individual spectra
#' rtime(spl)
setMethod("rtime", "Spectrum2List", function(object) {
    vapply(object, function(z) if(length(z@rt)) z@rt else NA_real_, numeric(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase precursorMz
#' 
#' @export
#'
#' @examples
#'
#' ## Extract the precursor m/z of each spectrum.
#' precursorMz(spl)
setMethod("precursorMz", "Spectrum2List", function(object) {
    vapply(object, function(z)
        if(length(z@precursorMz)) z@precursorMz else NA_real_, numeric(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase precursorCharge
#' 
#' @export
#'
#' @examples
#'
#' ## Extract the precursor charge of each spectrum.
#' precursorCharge(spl)
setMethod("precursorCharge", "Spectrum2List", function(object) {
    vapply(object, function(z)
        if(length(z@precursorCharge)) z@precursorCharge else NA_integer_,
        integer(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase precScanNum
#' 
#' @export
#'
#' @examples
#'
#' ## Extract the precursor scan number for each spectrum.
#' precScanNum(spl)
setMethod("precScanNum", "Spectrum2List", function(object) {
    vapply(object, function(z)
        if(length(z@precScanNum)) z@precScanNum else NA_integer_,
        integer(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase precursorIntensity
#' 
#' @export
#'
#' @examples
#'
#' ## Extract the precursor intensity of each spectrum.
#' precursorIntensity(spl)
setMethod("precursorIntensity", "Spectrum2List", function(object) {
    vapply(object, function(z)
        if(length(z@precursorIntensity)) z@precursorIntensity else NA_real_,
        numeric(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase acquisitionNum
#' 
#' @export
#'
#' @examples
#'
#' ## Extract the acquisition number of each spectrum.
#' acquisitionNum(spl)
setMethod("acquisitionNum", "Spectrum2List", function(object) {
    vapply(object, function(z)
        if(length(z@acquisitionNum)) z@acquisitionNum else NA_integer_,
        integer(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase scanIndex
#' 
#' @export
#'
#' @examples
#'
#' ## Extract the scan index of each spectrum.
#' scanIndex(spl)
setMethod("scanIndex", "Spectrum2List", function(object) {
    vapply(object, function(z)
        if(length(z@scanIndex)) z@scanIndex else NA_integer_,
        integer(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase peaksCount
#' 
#' @export
#'
#' @examples
#'
#' ## Get the number of peaks per spectrum.
#' peaksCount(spl)
setMethod("peaksCount", "Spectrum2List", function(object) {
    vapply(object, peaksCount, integer(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase msLevel
#' 
#' @export
#'
#' @examples
#'
#' ## Get the MS level of each spectrum.
#' msLevel(spl)
setMethod("msLevel", "Spectrum2List", function(object) {
    vapply(object, msLevel, integer(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase tic
#' 
#' @export
#'
#' @examples
#'
#' ## Get the total ion current for each spectrum.
#' tic(spl)
setMethod("tic", "Spectrum2List", function(object) {
    vapply(object, tic, numeric(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase ionCount
#' 
#' @export
#'
#' @examples
#'
#' ## Get the total ion current for each spectrum.
#' ionCount(spl)
setMethod("ionCount", "Spectrum2List", function(object) {
    vapply(object, ionCount, numeric(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase collisionEnergy
#' 
#' @export
#'
#' @examples
#'
#' ## Extract the collision energy for spectrum.
#' collisionEnergy(spl)
setMethod("collisionEnergy", "Spectrum2List", function(object) {
    vapply(object, function(z)
        if(length(z@collisionEnergy)) z@collisionEnergy else NA_real_,
        numeric(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase fromFile
#' 
#' @export
#'
#' @examples
#'
#' ## Extract the file index for the spectrum.
#' fromFile(spl)
setMethod("fromFile", "Spectrum2List", function(object) {
    vapply(object, function(z)
        if(length(z@fromFile)) z@fromFile else NA_integer_,
        integer(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase polarity
#' 
#' @export
#'
#' @examples
#'
#' ## Get the polarity for each spectrum.
#' polarity(spl)
setMethod("polarity", "Spectrum2List", function(object) {
    vapply(object, polarity, integer(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase smoothed
#' 
#' @export
#'
#' @examples
#'
#' ## Whether spectra are smoothed (i.e. processed with the `MSnbase::smooth`
#' ## function.
#' smoothed(spl)
setMethod("smoothed", "Spectrum2List", function(object) {
    vapply(object, smoothed, logical(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase isEmpty
#' 
#' @export
#'
#' @examples
#'
#' ## Are spectra empty (i.e. contain no peak data)?
#' isEmpty(spl)
setMethod("isEmpty", "Spectrum2List", function(x) {
    vapply(x, isEmpty, logical(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase centroided
#' 
#' @export
#'
#' @examples
#'
#' ## Do the spectra contain centroided data?
#' centroided(spl)
setMethod("centroided", "Spectrum2List", function(object) {
    vapply(object, centroided, logical(1))
})

#' @rdname Spectrum2List
#' 
#' @importMethodsFrom MSnbase isCentroided
#' 
#' @export
#'
#' @examples
#'
#' ## Do the spectra contain centroided data? Whether spectra are centroided
#' ## is estimated from the peak data.
#' isCentroided(spl)
setMethod("isCentroided", "Spectrum2List", function(object) {
    vapply(object, isCentroided, logical(1))
})

## clean, all = FALSE
## removePeaks, t
## trimMz, mzlim
## filterMz, mz
## quantify, method = c("trapezoidation", "max", "sum", reporters, strict = FALSE)
## normalize, method = c("max", "sum", "precursor", precursorIntensity)
## bin, binSize = 1L, breaks....
## pickPeaks...
## smooth...


## Fun things:
## compareSpectra
