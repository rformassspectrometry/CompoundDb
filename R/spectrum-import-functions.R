## Functions to import spectrum data.

#' @description
#'
#' Utility function to parse the data from a single spectrum xml file from
#' HMDB.
#'
#' @param x `character(1)` with the file path/name of the xml file.
#'
#' @param nonStop `logical(1)` whether content-related errors should be
#'     reported as a `warning`.
#'
#' @param collapsed `logical(1)` whether the returned `data.frame` should be
#'     *collapsed* or *expanded*. Collapsed means that the `data.frame` will
#'     have a single row and the m/z and intensity values are stored as a
#'     `list` in columns `"mz"` and `"intensity"`.
#' 
#' @return `data.frame`
#'
#' @author Johannes Rainer
#'
#' @return
#'
#' `data.frame` with as many rows as there are peaks and columns:
#' - external_spectrum_id (`character`): the HMDB-internal ID of the spectrum.
#' - compound_id (`character`): the HMDB ID the spectrum is associated with.
#' - polarity (`integer`): 0 for negative, 1 for positive, `NA` for not set.
#' - collision_energy (`numeric`): collision energy voltage.
#' - predicted (`logical`): whether the spectrum is predicted or experimentally
#'   verified.
#' - splash (`character`): the SPLASH key of the spectrum.
#' - instrument_type (`character`): the type of instrument on which the
#'   spectrum was measured.
#' - mz (`numeric` or `list`): m/z values of the spectrum.
#' - intensity (`numeric` or `list`): intensity of the spectrum.
#' 
#' @md
#'
#' @noRd
#'
#' @importFrom xml2 read_xml xml_text xml_find_first xml_find_all xml_double
.import_hmdb_ms_ms_spectrum <- function(x, nonStop = FALSE, collapsed = TRUE) {
    x_ml <- read_xml(x)
    id <- xml_text(xml_find_first(x_ml, "id"))
    cmp_id <- xml_text(xml_find_first(x_ml, "database-id"))
    if (id == "" || cmp_id == "") {
        msg <- paste0("Could not extract the HMDB ID from ", basename(x),
                      "! Is the file a spectrum xml file from HMDB?")
        if (nonStop) {
            warning(msg)
            return(data.frame())
        } else stop(msg)
    }
    plrty <- xml_text(xml_find_first(x_ml, "ionization-mode"))
    ## 0: negative, +1: positive, NA: not set.
    if (plrty == "")
        plrty <- NA_integer_
    else plrty <- ifelse(length(grep("pos", tolower(plrty))), yes = 1L, no = 0L)
    cev <- xml_double(xml_find_first(x_ml, "collision-energy-voltage"))
    prd <- xml_text(xml_find_first(x_ml, "predicted"))
    if (prd == "")
        prd <- NA
    else prd <- ifelse(prd == "false", yes = FALSE, no = TRUE)
    splsh <- xml_text(xml_find_first(x_ml, "splash-key"))
    itype <- xml_text(xml_find_first(x_ml, "instrument-type"))
    if (itype == "")
        itype <- NA_character_
    mz <- xml_double(xml_find_all(x_ml, "ms-ms-peaks/ms-ms-peak/mass-charge"))
    int <- xml_double(xml_find_all(x_ml, "ms-ms-peaks/ms-ms-peak/intensity"))
    if (!length(mz) | !length(int) | length(mz) != length(int)) {
        msg <- paste0("No mz and intensity values found in file ", basename(x))
        if (nonStop) {
            warning(msg)
            return(data.frame())
        } else stop(msg)
    }
    if (collapsed) {
        res <- data.frame(spectrum_id = id,
                          compound_id = cmp_id,
                          polarity = plrty,
                          collision_energy = cev,
                          predicted = prd,
                          splash = splsh,
                          instrument_type = itype,
                          stringsAsFactors = FALSE)
        res$mz <- list(mz)
        res$intensity <- list(int)
        res
    } else
        data.frame(spectrum_id = id,
                   compound_id = cmp_id,
                   polarity = plrty,
                   collision_energy = cev,
                   predicted = prd,
                   splash = splsh,
                   instrument_type = itype,
                   mz = mz,
                   intensity = int,
                   stringsAsFactors = FALSE)
}

#' @title Import MS/MS spectra from HMDB xml files
#'
#' @description
#'
#' `msms_spectra_hmdb` imports MS/MS spectra from corresponding xml files from
#' HMDB (http://www.hmdb.ca) and returns the data as a `data.frame`. HMDB
#' stores MS/MS spectrum data in xml files, one file per spectrum.
#'
#' Depending on the parameter `collapsed`, the returned `data.frame` is either
#' *collapsed*, meaning that each row represents data from one spectrum xml
#' file, or *expanded* with one row for each m/z and intensity pair for each
#' spectrum. Columns `"mz"` and `"intensity"` are of type `list` for
#' `collapsed = TRUE` and `numeric` for `collapsed = FALSE`.
#'
#' @note
#'
#' The HMDB xml files are supposed to be extracted from the downloaded zip file
#' into a folder and should not be renamed. The function identifies xml files
#' containing MS/MS spectra by their file name.
#'
#' The same spectrum ID can be associated with multiple compounds. Thus, the
#' function assignes an arbitrary ID (column `"spectrum_id"`) to values from
#' each file. The original ID of the spectrum in HMDB is provided in column
#' `"original_spectrum_id"`.
#' 
#' @param x `character(1)`: with the path to directory containing the xml files.
#'
#' @param collapsed `logical(1)` whether the returned `data.frame` should be
#'     *collapsed* or *expanded*. See description for more details.
#'
#' @return `data.frame` with as many rows as there are peaks and columns:
#' 
#' - spectrum_id (`character`): an arbitrary, unique ID identifying values
#'   from one xml file.
#' - original_spectrum_id (`character`): the HMDB-internal ID of the spectrum.
#' - compound_id (`character`): the HMDB compound ID the spectrum is associated
#'   with.
#' - polarity (`integer`): 0 for negative, 1 for positive, `NA` for not set.
#' - collision_energy (`numeric`): collision energy voltage.
#' - predicted (`logical`): whether the spectrum is predicted or experimentally
#'   verified.
#' - splash (`character`): the SPLASH (SPectraL hASH) key of the spectrum
#'   (Wohlgemuth 2016).
#' - instrument_type (`character`): the type of MS instrument on which the
#'   spectrum was measured.
#' - mz (`numeric` or `list` of `numeric`): m/z values of the spectrum.
#' - intensity (`numeric` or `list` of `numeric`): intensity of the spectrum.
#'
#' @md
#'
#' @author Johannes Rainer
#'
#' @export
#'
#' @seealso
#'
#' [Spectrum2List()] for converting the returned `data.frame` into
#' a [Spectrum2List] object (list of [Spectrum2] objects with annotations).
#'
#' [createCompDb()] for the function to create a [CompDb] database with
#' compound annotation and spectrum data.
#' 
#' @references
#'
#' Wohlgemuth G, Mehta SS, Mejia RF, Neumann S, Pedrosa D, Pluskal T,
#' Schymanski EL, Willighagen EL, Wilson M, Wishart DS, Arita M,
#' Dorrestein PC, Bandeira N, Wang M, Schulze T, Selak RM, Steinbeck C,
#' Nainala VC, Mistrik R, Nishioka T, Fiehn O. SPLASH, A hashed identifier for
#' mass spectra. Nature Biotechnology 2016 34(11):1099-1101
#' 
#' @examples
#'
#' ## Locate the folder within the package containing test xml files.
#' pth <- system.file("xml", package = "CompoundDb")
#'
#' ## List all files in that directory
#' dir(pth)
#'
#' ## Import spectrum data from HMDB MS/MS spectrum xml files in that directory
#' msms_spectra_hmdb(pth)
#'
#' ## Import the data as an *expanded* data frame, i.e. with a row for each
#' ## single m/z (intensity) value.
#' msms_spectra_hmdb(pth, collapsed = FALSE)
msms_spectra_hmdb <- function(x, collapsed = TRUE) {
    fls <- dir(x, pattern = "ms_ms_spectrum(.)+xml$", full.names = TRUE)
    if (!length(fls))
        stop("Unable to find any MS/MS spectrum xml files from HMDB in ",
             "folder ", x)
    message("Going to process ", length(fls), " xml files.")
    res <- do.call(rbind, lapply(fls, .import_hmdb_ms_ms_spectrum,
                                 nonStop = TRUE, collapsed = collapsed))
    ## Assign an arbitrary spectrum ID.
    message("Postprocessing data ...", appendLF = FALSE)
    colnames(res)[colnames(res) == "spectrum_id"] <- "original_spectrum_id"
    ids <- factor(paste0(res$compound_id, "-", res$original_spectrum_id))
    res$spectrum_id <- as.character(as.numeric(ids))
    message("OK")
    res
}

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
    x_new <- unique(x[, !colnames(x) %in% c("mz", "intensity")])
    x_new$mz <- list(x$mz)
    x_new$intensity <- list(x$intensity)
    x_new
}

#' @description
#'
#' Function to expand a collapsed spectrum `data.frame` (such as collapsed with
#' the `.collapse_spectrum_df` function).
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
    x_exp
}


## Function to import spectrum data from MoNa etc.
