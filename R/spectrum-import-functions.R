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
    else plrty <- if(length(grep("pos", plrty, ignore.case = TRUE))) 1L else 0L
    cev <- xml_double(xml_find_first(x_ml, "collision-energy-voltage"))
    prd <- xml_text(xml_find_first(x_ml, "predicted"))
    if (prd == "")
        prd <- NA
    else prd <- prd != "false"
    splsh <- xml_text(xml_find_first(x_ml, "splash-key"))
    itype <- xml_text(xml_find_first(x_ml, "instrument-type"))
    if (itype == "")
        itype <- NA_character_
    instr <- xml_text(xml_find_first(x_ml, "notes"))
    if (instr == "")
        instr <- NA_character_
    else instr <- .extract_field_from_string(instr, "instrument=", "\n")
    mz <- xml_double(xml_find_all(x_ml, "ms-ms-peaks/ms-ms-peak/mass-charge"))
    int <- xml_double(xml_find_all(x_ml, "ms-ms-peaks/ms-ms-peak/intensity"))
    if (!length(mz) || !length(int) || length(mz) != length(int)) {
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
                          instrument = instr,
                          precursor_mz = NA_real_,
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
                   instrument = instr,
                   precursor_mz = NA_real_,
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
#' - spectrum_id (`integer`): an arbitrary, unique ID identifying values
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
#' - instrument (`character`): the MS instrument (not available for all spectra
#'   in HMDB).
#' - precursor_mz (`numeric`): not provided by HMDB and thus `NA`.
#' - mz (`numeric` or `list` of `numeric`): m/z values of the spectrum.
#' - intensity (`numeric` or `list` of `numeric`): intensity of the spectrum.
#'
#' @md
#'
#' @author Johannes Rainer
#'
#' @export
#'
#' @family spectrum data import functions.
#'
#' @seealso
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
    message("Postprocessing data ... ", appendLF = FALSE)
    colnames(res)[colnames(res) == "spectrum_id"] <- "original_spectrum_id"
    ids <- factor(paste0(res$compound_id, "-", res$original_spectrum_id))
    res$spectrum_id <- as.integer(ids)
    message("OK")
    res
}

#' @title Import MS/MS spectra from MoNa
#'
#' @description
#'
#' `msms_spectra_mona` imports MS/MS spectra from a MoNa (Massbank of North
#' America, http://mona.fiehnlab.ucdavis.edu/downloads) SDF file and returns
#' the data as a `data.frame`.
#'
#' Depending on the parameter `collapsed`, the returned `data.frame` is either
#' *collapsed*, meaning that each row represents data from one spectrum,
#' or *expanded* with one row for each m/z and intensity pair for each
#' spectrum. Columns `"mz"` and `"intensity"` are of type `list` for
#' `collapsed = TRUE` and `numeric` for `collapsed = FALSE`.
#'
#' @note
#'
#' The identifiers provided by MoNa are used as *spectrum_id*. Note also that
#' the MoNa data is not normalized in the sense that each spectrum is
#' associated to one compound and the compound data is partially redundant.
#' Also, MoNa does not provide a *splash* for a spectrum, hence the
#' corresponding column will only contain `NA`.
#'
#' @param x `character(1)`: with the path to directory containing the xml files.
#'
#' @param collapsed `logical(1)` whether the returned `data.frame` should be
#'     *collapsed* or *expanded*. See description for more details.
#'
#' @return `data.frame` with as many rows as there are peaks and columns:
#'
#' - spectrum_id (`integer`): an arbitrary, unique ID for each spectrum.
#' - original_spectrum_id (`character`): The ID from the spectrum as specified
#'   in the MoNa SDF.
#' - compound_id (`character`): the compound ID the spectrum is associated
#'   with.
#' - polarity (`integer`): 0 for negative, 1 for positive, `NA` for not set.
#' - collision_energy (`character`): collision energy voltage.
#' - predicted (`logical`): whether the spectrum is predicted or experimentally
#'   verified.
#' - splash (`character`): `NA` since SPLASH (SPectraL hASH) keys are not
#'   provided.
#' - instrument_type (`character`): the type of MS instrument on which the
#'   spectrum was measured.
#' - instrument (`character`): the MS instrument.
#' - precursor_mz (`numeric`): precursor m/z.
#' - adduct (`character`): ion formed from the precursor ion.
#' - ms_level (`integer`): stage of the sequential mass spectrometry (MSn).
#' - mz (`numeric` or `list` of `numeric`): m/z values of the spectrum.
#' - intensity (`numeric` or `list` of `numeric`): intensity of the spectrum.
#'
#' @author Johannes Rainer
#'
#' @export
#'
#' @family spectrum data import functions.
#'
#' @seealso
#'
#' [createCompDb()] for the function to create a [CompDb] database with
#' compound annotation and spectrum data.
#'
#' @examples
#'
#' ## Define the test file containing the data
#' fl <- system.file("sdf/MoNa_export-All_Spectra_sub.sdf.gz",
#'     package = "CompoundDb")
#' ## Import spectrum data from the SDF file with a subset of the MoNa data
#' msms_spectra_mona(fl)
#'
#' ## Import the data as an *expanded* data frame, i.e. with a row for each
#' ## single m/z (intensity) value.
#' msms_spectra_mona(fl, collapsed = FALSE)
msms_spectra_mona <- function(x, collapsed = TRUE) {
    spctra <- .import_mona_sdf(x, TRUE, FALSE)$msms_spectrum
    colnames(spctra)[colnames(spctra) == "precursor_type"] <- "adduct"
    spec_trim <- gsub(" ", "", spctra$spectrum_type, fixed = TRUE)
    is_mslevel <- grepl("^ms\\d+$", spec_trim, ignore.case = TRUE)
    if(any(is_mslevel)){
        spctra$ms_level <- rep(NA_integer_, length(spec_trim))
        spctra$ms_level[is_mslevel] <- as.integer(
            sub("ms", "", spec_trim[is_mslevel], ignore.case = TRUE))
    }
    if (collapsed) spctra
    else .expand_spectrum_df(spctra)
}

#' Extract spectrum data from a MoNa SDF file
#'
#' @author Johannes Rainer
#'
#' @noRd
.extract_spectra_mona_sdf <- function(x) {
    n <- nrow(x)
    plrty <- rep(NA_integer_, n)
    plrty[grep("P", x[, "ION MODE"])] <- 1L
    plrty[grep("N", x[, "ION MODE"])] <- 0L
    mzint <- lapply(strsplit(x[, "MASS SPECTRAL PEAKS"], " | __ "),
                    function(z) matrix(as.numeric(z), ncol = 2, byrow = TRUE))
    res <- data.frame(original_spectrum_id = x[, "ID"],
               compound_id = .compound_id_from_mona_sdf(x),
               polarity = plrty,
               collision_energy = x[, "COLLISION ENERGY"],
               predicted = NA,
               splash = NA_character_,
               instrument_type = x[, "INSTRUMENT TYPE"],
               instrument = x[, "INSTRUMENT"],
               precursor_mz = as.numeric(x[, "PRECURSOR M/Z"]),
               precursor_type = x[, "PRECURSOR TYPE"],
               spectrum_type = x[, "SPECTRUM TYPE"],
               spectrum_id = seq_len(nrow(x)),
               stringsAsFactors = FALSE)
    res$mz <- lapply(mzint, function(z) z[, 1])
    res$intensity <- lapply(mzint, function(z) z[, 2])
    res
}

#' Create compound IDs from a mona SDF. For now we create simply IDs from 1
#' to nrow(x), in future we might use this to *normalize* the matrix and
#' generate compound IDs e.g. on unique InChI keys or similar.
#'
#' @noRd
#'
#' @author Johannes Rainer
.compound_id_from_mona_sdf <- function(x, prefix = "CMP") {
    n <- nrow(x)
    sprintf(paste0(prefix, "%0", ceiling(log10(n + 1)), "d"), seq_len(n))
}
