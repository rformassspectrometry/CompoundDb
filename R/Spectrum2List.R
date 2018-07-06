
#' @name Spectrum2List
#'
#' @aliases Spectrum2List-class show,Spectrum2List-method
#'
#' @title List of Spectrum 2 objects along with annotations
#'
#' @description
#'
#' `Spectrum2List` objects allow to collect one or more [Spectrum2] object(s)
#' in a `list`-like structure with the additional possibility to add arbitrary
#' annotations to each individual [Spectrum2] object. These can be accessed/set
#' with the [mcols] method.
#'
#' @md
#'
#' @seealso
#'
#' [Spectrum2] in the `MSnbase` package.
#'
#' [msms_spectra_hmdb()] for a function to import a `data.frame` from files
#' containing spectrum MS/MS information from HMDB (http://www.hmdb.ca) in
#' xml format.
#'
#' @rdname Spectrum2List
NULL

#' @importClassesFrom MSnbase Spectrum2
#' 
#' @importClassesFrom S4Vectors SimpleList
#'
#' @importFrom S4Vectors SimpleList
#'
#' @exportClass Spectrum2List
.Spectrum2List <- setClass("Spectrum2List",
                           contains = "SimpleList",
                           prototype = prototype(elementType = "Spectrum2")
                           )

setValidity("Spectrum2List", function(object) {
    ## All elements in the list have to be Spectrum2 objects.
    msg <- character()
    if (any(vapply(object, function(z) !is(z, "Spectrum2"), logical(1))))
        msg <- c(msg, "All elements have to be Spetrum2 objects")
    if (length(msg)) msg else TRUE
})

#' @importMethodsFrom S4Vectors mcols
.show_Spectrum2List <- function(x, margin = "", print.classinfo = FALSE) {
    cat("Spectrum2List with", length(x), "spectra and", ncol(mcols(x)),
        " metadata column(s):\n")
    out <- S4Vectors:::makePrettyMatrixForCompactPrinting(
                           x, .make_naked_matrix_from_Spectrum2List)
    if (print.classinfo) {
        .COL2CLASS <- c(msLevel = "integer", rtime = "numeric",
                        precursorMz = "numeric", peaksCount = "integer")
        classinfo <- S4Vectors:::makeClassinfoRowForCompactPrinting(x,
                                                                    .COL2CLASS)
        out <- rbind(classinfo, out)
    }
    if (nrow(out) != 0L)
        rownames(out) <- paste0(margin, rownames(out))
    print(out, quote = FALSE, right = TRUE, max = length(out))
}

#' @importMethodsFrom S4Vectors showAsCell
#'
#' @importMethodsFrom MSnbase msLevel rtime precursorMz peaksCount
.make_naked_matrix_from_Spectrum2List <- function(x) {
    x_len <- length(x)
    mcls <- mcols(x, use.names = FALSE)
    x_mcls_len <- if (is.null(mcls)) 0L else ncol(mcls)
    res <- cbind(msLevel = as.character(unlist(lapply(x, msLevel))),
                 rtime = format(unlist(lapply(x, rtime)), digits = 6),
                 precursorMz = format(unlist(lapply(x, precursorMz)),
                                      digits = 6),
                 peaksCount = as.character(unlist(lapply(x, peaksCount))))
    if (!any(colnames(res) == "rtime"))
        res <- cbind(res, rtime = rep(NA, nrow(res)))
    if (!any(colnames(res) == "precursorMz"))
        res <- cbind(res, precursorMz = rep(NA, nrow(res)))
    res <- res[, c("msLevel", "rtime", "precursorMz", "peaksCount"),
               drop = FALSE]
    if (x_mcls_len > 0) {
        tmp <- do.call(data.frame, c(lapply(mcls, showAsCell),
                                     list(check.names = FALSE)))
        res <- cbind(res, `|` = rep.int("|", x_len), as.matrix(tmp))
    }
    res
}

#' @importFrom methods show
#'
#' @export
setMethod("show", "Spectrum2List", function(object) {
    .show_Spectrum2List(object, margin = "  ", print.classinfo = TRUE)
})

#' @rdname Spectrum2List
#'
#' @section Constructor:
#'
#' New [Spectrum2List] can be created with the `Spectrum2List(...)` function
#' where `...` can either be a single [Spectrum2] object, a `list` of
#' [Spectrum2] objects or a `data.frame`. If a `data.frame` is provided each
#' row is expected to represent one peak and the following required columns
#' have to be present:
#'
#' - spectrum_id: (`character`) defining the 
#' - mz: (`numeric`) the m/z value of each peak.
#' - intensity: (`numeric`) the peak's intensity.
#'
#' Optional columns listed below are used to fill the indicated slots of the
#' [Spectrum2] object:
#'
#' - polarity: `polarity`.
#' - ms_level: `msLevel`.
#' - rt: `rt`
#' - precursor_mz: `precursorMz`.
#' - precursor_charge: `precursorCharge`.
#' - precursor_intensity: `precursorIntensity`.
#' - collision_energy: `collisionEnergy`.
#' - acquisition_num: `acquisitionNum`.
#' - scan_index: `scanIndex`.
#' - from_file: `fromFile`.
#' - precursor_scan_num: `precScanNum`.
#' 
#' Any additional columns are added to the object's metadata columns `mcols`.
#'
#' @param ... [Spectrum2] object(s), a `list` of [Spectrum2] objects or a
#'     `data.frame`. See below for details.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @export
#'
#' @examples
#' 
#' ## Create from Spectrum2 objects
#' sp1 <- new("Spectrum2", mz = c(1, 2, 4), intensity = c(4, 5, 2))
#' sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5))
#'
#' spl <- Spectrum2List(sp1, sp2)
#' 
#' ## Build from data.frame
#' df <- data.frame(spectrum_id = c("b", "b", "b", "b", "a", "a", "a"),
#'     mz = c(1, 2, 3, 4, 1, 2, 4), intensity = c(5, 3, 2, 5, 4, 5, 2),
#'     polarity = c(1, 1, 1, 1, 0, 0, 0), compound_id = rep("cp_1", 7),
#'     stringsAsFactors = FALSE)
#'
#' spl <- Spectrum2List(df)
#' polarity(spl[[1]])
#' 
#' ## Additional annotations are provided in mcols
#' mcols(spl)
Spectrum2List <- function(...) {
    args <- list(...)
    if (length(args) == 1L && is.list(args[[1L]]))
        args <- args[[1L]]
    ## If we've got a data.frame check if it's in the expected format.
    if (is.data.frame(args)) {
        res <- .spectra2_from_df(args)
        args <- res$spectra
        mcls <- res$mcols
    } else mcls <- NULL
    new("Spectrum2List", listData = args, elementMetadata = mcls)
}

#' @description Create a list of `Spectrum2` objects from a `data.frame`.
#'
#' @param x `data.frame` with spectrum data.
#' 
#' @author Johannes Rainer
#'
#' @return `list` with elements `"spectra"` containing the `list` of
#'     `Spectrum2` objects and `mcols` with the metadata columns not used
#'     for the `Spectrum2` object creation.
#' @noRd
#'
#' @md
#'
#' @importFrom S4Vectors DataFrame
#' 
#' @examples
#'
#' df <- data.frame(spectrum_id = c("a", "a", "b"), mz = c(1, 2, 1),
#'     intensity = c(2, 3, 5), comp_id = "Z", polarity = 1,
#'     stringsAsFactors = FALSE)
#'
#' res <- CompoundDb:::.spectra2_from_df(df)
.spectra2_from_df <- function(x) {
    ## Check required columns.
    req_cols <- c("spectrum_id", "mz", "intensity")
    if (!all(req_cols %in% colnames(x)))
        stop("required columns 'spectrum_id', 'mz' and 'intensity' missing")
    supp_cols <- c(req_cols, "polarity", "ms_level", "rt", "precursor_mz",
                   "precursor_charge", "precursor_intensity",
                   "collision_energy")
    mz <- x$mz
    int <- x$intensity
    nvals <- lengths(split(x$mz, factor(x$spectrum_id,
                                          levels = unique(x$spectrum_id))))
    x <- unique(x[, !(colnames(x) %in% c("mz", "intensity"))])
    if (nrow(x) != length(unique(x$spectrum_id)))
        stop("Unexpected number of rows in data.frame")
    ## Process optional columns.
    ## polarity -> polarity
    if (length(colnames(x)) && any(colnames(x) == "polarity")) {
        polarity <- x$polarity
        x <- x[, colnames(x) != "polarity"]
    } else polarity <- integer()
    ## rt -> rt
    if (length(colnames(x)) && any(colnames(x) == "rt")) {
        rt <- x$rt
        x <- x[, colnames(x) != "rt"]
    } else rt <- numeric()
    ## ms_level -> msLevel
    if (length(colnames(x)) && any(colnames(x) == "ms_level")) {
        msLevel <- x$ms_level
        x <- x[, colnames(x) != "ms_level"]
    } else msLevel <- rep(2L, nrow(x))
    ## precursor_mz -> precursorMz
    if (length(colnames(x)) && any(colnames(x) == "precursor_mz")) {
        precursorMz <- x$precursor_mz
        x <- x[, colnames(x) != "precursor_mz"]
    } else precursorMz <- numeric()
    ## precursor_charge -> precursorCharge
    if (length(colnames(x)) && any(colnames(x) == "precursor_charge")) {
        precursorCharge <- x$precursor_charge
        x <- x[, colnames(x) != "precursor_charge"]
    } else precursorCharge <- integer()
    ## precursor_intensity -> precursorIntensity
    if (length(colnames(x)) && any(colnames(x) == "precursor_intensity")) {
        precursorIntensity <- x$precursor_intensity
        x <- x[, colnames(x) != "precursor_intensity"]
    } else precursorIntensity <- numeric()
    ## collision_energy -> collisionEnergy
    if (length(colnames(x)) && any(colnames(x) == "collision_energy")) {
        collisionEnergy <- x$collision_energy
        x <- x[, colnames(x) != "collision_energy"]
    } else collisionEnergy <- numeric()
    ## acquisition_num -> acquisitionNum
    if (length(colnames(x)) && any(colnames(x) == "acquisition_num")) {
        acquisitionNum <- x$acquisition_num
        x <- x[, colnames(x) != "acquisition_num"]
    } else acquisitionNum <- integer()
    ## scan_index -> scanIndex
    if (length(colnames(x)) && any(colnames(x) == "scan_index")) {
        scanIndex <- x$scan_index
        x <- x[, colnames(x) != "scan_index"]
    } else scanIndex <- integer()
    ## from_file -> fromFile
    if (length(colnames(x)) && any(colnames(x) == "from_file")) {
        fromFile <- x$from_file
        x <- x[, colnames(x) != "from_file"]
    } else fromFile <- integer()
    ## precursor_scan_num -> precScanNum
    if (length(colnames(x)) && any(colnames(x) == "precursor_scan_num")) {
        precScanNum <- x$precursor_scan_num
        x <- x[, colnames(x) != "precursor_scan_num"]
    } else precScanNum <- integer()
    
    ## Create the spectra
    spl <- MSnbase:::Spectra2_mz_sorted(peaksCount = nvals, rt = rt,
                                        acquisitionNum = acquisitionNum,
                                        scanIndex = scanIndex, mz = mz,
                                        intensity = int, fromFile = fromFile,
                                        polarity = polarity, msLevel = msLevel,
                                        precScanNum = precScanNum,
                                        precursorMz = precursorMz,
                                        precursorIntensity = precursorIntensity,
                                        precursorCharge = precursorCharge,
                                        collisionEnergy = collisionEnergy,
                                        nvalue = nvals)
    list(spectra = spl, mcols = DataFrame(x))
}
