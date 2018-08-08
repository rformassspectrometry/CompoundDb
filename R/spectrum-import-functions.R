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
#' @return `data.frame`
#'
#' @author Johannes Rainer
#'
#' @return
#'
#' `data.frame` with as many rows as there are peaks and columns:
#'
#' - spectrum_id (`character`): the HMDB-internal ID of the spectrum.
#' - compound_id (`character`): the HMDB ID the spectrum is associated with.
#' - polarity (`integer`): 0 for negative, 1 for positive, `NA` for not set.
#' - collision_energy (`numeric`): collision energy voltage.
#' - predicted (`logical`): whether the spectrum is predicted or experimentally
#'   verified.
#' - splash (`character`): the SPLASH key of the spectrum.
#' - instrument_type (`character`): the type of instrument on which the
#'   spectrum was measured.
#' - mz (`numeric`): m/z values of the spectrum.
#' - intensity (`numeric`): intensity of the spectrum.
#' 
#' @md
#'
#' @noRd
#'
#' @importFrom xml2 read_xml xml_text xml_find_first xml_find_all xml_double
.import_hmdb_ms_ms_spectrum <- function(x, nonStop = FALSE) {
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
    ## Return result.
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
#' HMDB (http://www.hmdb.ca). HMDB stores MS/MS spectrum data in xml files, one
#' file per spectrum.
#'
#' @note
#'
#' The HMDB xml files are supposed to be extracted from the downloaded zip file
#' into a folder and should not be renamed. The function identifies xml files
#' containing MS/MS spectra by their file name.
#' 
#' @param x `character(1)`: with the path to directory containing the xml files.
#'
#' @return `data.frame` with as many rows as there are peaks and columns:
#' 
#' - spectrum_id (`character`): the HMDB-internal ID of the spectrum.
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
#' - mz (`numeric`): m/z values of the spectrum.
#' - intensity (`numeric`): intensity of the spectrum.
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
msms_spectra_hmdb <- function(x) {
    fls <- dir(x, pattern = "ms_ms_spectrum(.)+xml$", full.names = TRUE)
    if (!length(fls))
        stop("Unable to find any MS/MS spectrum xml files from HMDB in ",
             "folder ", x)
    message("Going to process ", length(fls), " xml files.")
    ## for (i in 15383:length(fls)) {
    ##     tmp <- CompoundDb:::.import_hmdb_ms_ms_spectrum(fls[i])
    ## }
    do.call(rbind, lapply(fls, .import_hmdb_ms_ms_spectrum, nonStop = TRUE))
}

## Function to import spectrum data from MoNa etc.
