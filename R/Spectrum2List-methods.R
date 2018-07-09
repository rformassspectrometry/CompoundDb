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
#' `Spectrum2` objects. For more details please refer to the corresoonding
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

## precursorMz
## precursorCharge
## precursorIntensity
## acquisitionNum
## scanIndex
## precScanNum
## peaksCount
## mzLevel
## collisionEnergy
## tic
## ionCount
## fromFile
## polarity
## centroided
## smoothed
## isEmpty
## isCentroided

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
