library(CompoundDb)

#' @title Create a CompDb from MassBank data
#'
#' @description
#'
#' Create a `CompDb` database with the data from a `MassBank` database.
#'
#' @param con connection to a MassBank database.
#'
#' @param path `character(1)` defining the path where the `CompDb` SQLite
#'     database should be stored.
#'
#' @author Johannes Rainer
massbank_to_compdb <- function(con, path = ".") {
    message("Getting compound information ... ", appendLF = FALSE)
    compound_data <- dbGetQuery(con, "select * from ms_compound")
    name_data <- dbGetQuery(con, "select * from synonym")
    names <- split(name_data$synonym,
                   factor(name_data$compound_id,
                          levels = compound_data$compound_id))
    compound_data$name <- vapply(names, function(z) z[1], character(1))
    compound_data$synonyms <- names
    message("Done\nGetting MS2 spectra data ... ", appendLF = FALSE)
    spectra_data <- dbGetQuery(con, "select * from msms_spectrum")
    spectra_data$compound_id <- as.character(spectra_data$compound_id)
    spectra_data$precursor_mz <- as.numeric(spectra_data$precursor_mz_text)
    spectra_data$precursor_intensity <-
        as.numeric(spectra_data$precursor_intensity)
    spectra_data$original_spectrum_id <- spectra_data$spectrum_id
    spectra_data$spectrum_id <- seq_len(nrow(spectra_data))
    pol <- rep(NA_integer_, nrow(spectra_data))
    pol[grep("POS", spectra_data$polarity)] <- 1L
    pol[grep("NEG", spectra_data$polarity)] <- 0L
    spectra_data$polarity <- pol
    spectra_data$ms_level <- as.integer(sub("MS", "", spectra_data$ms_level))
    message("Done\nGetting MS2 peak data ... ", appendLF = FALSE)
    pks <- dbGetQuery(con, "select * from msms_spectrum_peak")
    mzs <- split(pks$mz, pks$spectrum_id)
    ints <- split(pks$intensity, pks$spectrum_id)
    ## Ensure they are correctly ordered.
    idx <- match(names(mzs), spectra_data$original_spectrum_id)
    spectra_data$mz <- unname(mzs[idx])
    spectra_data$intensity <- unname(ints[idx])
    message("Done")
    version <- dbGetQuery(con, "select * from LAST_UPDATE")
    metad <- make_metadata(source = "MassBank",
                           source_version = sub("\n", "", version$VERSION[1]),
                           source_date = format(version$LAST_UPDATE,
                                                format = "%Y-%m-%d"),
                           url = "https://massbank.eu/MassBank/")
    message("Creating database ... ", appendLF = FALSE)
    db <- createCompDb(compound_data, metadata = metad,
                       msms_spectra = spectra_data,
                       path = path)
    message("Done")
    db
}
