#' @title Extract compound data from a file in SDF format
#'
#' @description
#'
#' `compound_tbl_sdf` extracts basic compound annotations from a file in SDF
#' format (structure-data file). The function currently supports SDF files from:
#' + HMDB (Human Metabolome Database): http://www.hmdb.ca
#' + ChEBI (Chemical Entities of Biological Interest): http://ebi.ac.uk/chebi
#' + LMSD (LIPID MAPS Structure Database): http://www.lipidmaps.org
#' + PubChem: https://pubchem.ncbi.nlm.nih.gov/
#' + MoNa: http://mona.fiehnlab.ucdavis.edu/ (see notes below!)
#'
#' @details
#'
#' Column `"compound_name"` reports for HMDB files the `"GENERIC_NAME"`, for
#' ChEBI the `"ChEBI Name"`, for PubChem the `"PUBCHEM_IUPAC_TRADITIONAL_NAME"`,
#' and for Lipid Maps the `"COMMON_NAME"`, if that is
#' not available, the first of the compounds synonyms and, if that is also not
#' provided, the `"SYSTEMATIC_NAME"`.
#'
#' @note
#'
#' `compound_tbl_sdf` supports also to read/process gzipped files.
#'
#' MoNa SDF files organize the data by individual spectra (i.e. each element
#' is one spectrum) and individual compounds can not easily and consistently
#' defined (i.e. not all entries have an InChI ID or other means to uniquely
#' identify compounds). Thus, the function returns a highly redundant compount
#' table. Feedback on how to reduce this redundancy would be highly welcome!
#'
#' @param file `character(1)` with the name of the SDF file.
#'
#' @param collapse optional `character(1)` to be used to collapse multiple
#'     values in the columns `"synonyms"`. See examples for details.
#'
#' @return A [tibble::tibble] with general compound information (one row per
#' compound):
#' + `compound_id`: the ID of the compound.
#' + `compound_name`: the compound's name.
#' + `inchi`: the InChI of the compound.
#' + `inchikey`: the InChI key.
#' + `formula`: the chemical formula.
#' + `mass`: the compound's mass.
#' + `synonyms`: the compound's synonyms (aliases). This type of this column is
#'   by default a `list` to support multiple aliases per compound, unless
#'   argument `collapse` is provided, in which case multiple synonyms are pasted
#'   into a single element separated by the value of `collapse`.
#'
#' @family compound table creation functions
#'
#' @author Johannes Rainer and Jan Stanstrup
#'
#' @export
#'
#' @md
#'
#' @seealso [createCompDb()] for a function to create a SQLite-based compound
#'     database.
#'
#' @importFrom ChemmineR read.SDFset datablock datablock2ma
#'
#' @examples
#'
#' ## Read compound information from a subset of HMDB
#' fl <- system.file("sdf/HMDB_sub.sdf.gz", package = "CompoundDb")
#' cmps <- compound_tbl_sdf(fl)
#' cmps
#'
#' ## Column synonyms contains a list
#' cmps$synonyms
#'
#' ## If we provide the optional argument collapse, multiple entries will be
#' ## collapsed.
#' cmps <- compound_tbl_sdf(fl, collapse = "|")
#' cmps
#' cmps$synonyms
compound_tbl_sdf <- function(file, collapse) {
    if (missing(file))
        stop("Please provide the file name using 'file'")
    if (!file.exists(file))
        stop("Can not fine file ", file)
    res <- .simple_extract_compounds_sdf(
        datablock2ma(datablock(read.SDFset(file))))
    if (!missing(collapse)) {
        ## collapse elements from lists.
        res$synonyms <- vapply(res$synonyms, paste0, collapse = collapse,
                               FUN.VALUE = "character")
    }
    res
}

#' @title Extract compound data from LipidBlast
#'
#' @description
#'
#' `compound_tbl_lipidblast` extracts basic comopund annotations from a
#' LipidBlast file in (json format) downloaded from
#' http://mona.fiehnlab.ucdavis.edu/downloads
#'
#' @param file `character(1)` with the name of the file name.
#'
#' @param collapse optional `character(1)` to be used to collapse multiple
#'     values in the columns `"synonyms"`. See examples for details.
#'
#' @return A [tibble::tibble] with general compound information (one row per
#' compound):
#' + `compound_id`: the ID of the compound.
#' + `compound_name`: the compound's name.
#' + `inchi`: the inchi of the compound.
#' + `formula`: the chemical formula.
#' + `mass`: the compound's mass.
#' + `synonyms`: the compound's synonyms (aliases). This type of this column is
#'   by default a `list` to support multiple aliases per compound, unless
#'   argument `collapse` is provided, in which case multiple synonyms are pasted
#'   into a single element separated by the value of `collapse`.
#'
#' @family compound table creation functions
#'
#' @author Johannes Rainer and Jan Stanstrup
#'
#' @export
#'
#' @md
#'
#' @examples
#'
#' ## Read compound information from a subset of HMDB
#' fl <- system.file("json/MoNa-LipidBlast_sub.json", package = "CompoundDb")
#' cmps <- compound_tbl_lipidblast(fl)
#' cmps
compound_tbl_lipidblast <- function(file, collapse) {
    if (missing(file))
        stop("Please provide the file name using 'file'")
    if (!file.exists(file))
        stop("Can not fine file ", file)
    res <- .import_lipidblast(file)
    if (!missing(collapse)) {
        ## collapse elements from lists.
        res$synonyms <- vapply(res$synonyms, paste0, collapse = collapse,
                               FUN.VALUE = "character")
    }
    res
}

#' @description
#'
#' Internal function to extract compound information from a file in SDF format.
#'
#' @param x what is returned by datablock2ma(datablock(read.SDFset)).
#'
#' @return A [tibble::tibble] with columns `"compound_id"`, `"compound_name"`,
#'     `"inchi"`, `"formula"`, `"mass"`.
#'
#' @importFrom tibble data_frame
#'
#' @md
#'
#' @author Johannes Rainer
#'
#' @noRd
.simple_extract_compounds_sdf <- function(x) {
    source_db <- .guess_sdf_source(colnames(x))
    if (is.null(source_db))
        stop("The SDF file is not supported. Supported are SDF files from ",
             "HMDB, ChEBI and LipidMaps.")
    colmap <- get(paste0(".", source_db, "_colmap"))
    sep <- get(paste0(".", source_db, "_separator"))
    ## Fix missing COMMON_NAME entries in LipidMaps (see issue #1)
    nms <- x[, colmap["name"]]
    syns <- strsplit(x[, colmap["synonyms"]], split = sep)
    if (source_db == "lipidmaps") {
        nas <- is.na(nms)
        if (any(nas))
            nms[nas] <- vapply(syns[nas], `[[`, 1, FUN.VALUE = "character")
        nas <- is.na(nms)
        if (any(nas))
            nms[nas] <- x[nas, "SYSTEMATIC_NAME"]
    }
    res <- data_frame(compound_id = x[, colmap["id"]],
                      compound_name = nms,
                      inchi = x[, colmap["inchi"]],
                      inchi_key = x[, colmap["inchi_key"]],
                      formula = x[, colmap["formula"]],
                      mass = as.numeric(x[, colmap["mass"]]),
                      synonyms = syns
                      )
    if (source_db == "mona") {
        warning("MoNa data can currently not be normalized and the ",
                "compound table contains thus highly redundant data.",
                call. = FALSE)
        ## Eventually reduce and normalize.
        inchis <- .extract_field_from_string(x[, "COMMENT"], "InChI")
        nona <- !is.na(inchis)
        inchis[nona] <- paste0("InChI", inchis[nona])
        res$inchi <- inchis
        res$compound_id <- .compound_id_from_mona_sdf(x)
        res
    } else res
}

#' @description
#'
#' Based on the provided `colnames` guess whether the file is from HMDB,
#' ChEBI or LipidBlast.
#'
#' @param x `character` with the column names of the data table.
#'
#' @return `character(1)` with the name of the resource or `NULL`.
#'
#' @md
#'
#' @author Johannes Rainer
#'
#' @noRd
.guess_sdf_source <- function(x) {
    if (any(x == "HMDB_ID"))
        return("hmdb")
    if (any(x == "ChEBI ID"))
        return("chebi")
    if (any(x == "LM_ID"))
        return("lipidmaps")
    if (any(x == "PUBCHEM_COMPOUND_CID"))
        return("pubchem")
    if (any(x == "ID")) {
        message("Guessing we've got a SDF file from MoNa.")
        return("mona")
    }
    NULL
}


.hmdb_colmap <- c(id = "HMDB_ID",
                  name = "GENERIC_NAME",
                  inchi = "INCHI_IDENTIFIER",
                  inchi_key = "INCHI_KEY",
                  formula = "FORMULA",
                  mass = "EXACT_MASS",
                  synonyms = "SYNONYMS"
                  )
.hmdb_separator <- "; "
.chebi_colmap <- c(id = "ChEBI ID",
                   name = "ChEBI Name",
                   inchi = "InChI",
                   inchi_key = "InChIKey",
                   formula = "Formulae",
                   mass = "Monoisotopic Mass",
                   synonyms = "Synonyms"
                   )
.chebi_separator <- " __ "
.lipidmaps_colmap <- c(id = "LM_ID",
                       name = "COMMON_NAME",
                       inchi = "INCHI",
                       inchi_key = "INCHI_KEY",
                       formula = "FORMULA",
                       mass = "EXACT_MASS",
                       synonyms = "SYNONYMS"
                       )
.lipidmaps_separator <- "; "
.pubchem_colmap <- c(id = "PUBCHEM_COMPOUND_CID",
                     name = "PUBCHEM_IUPAC_TRADITIONAL_NAME",
                     inchi = "PUBCHEM_IUPAC_INCHI",
                     inchi_key = "PUBCHEM_IUPAC_INCHIKEY",
                     formula = "PUBCHEM_MOLECULAR_FORMULA",
                     mass = "PUBCHEM_EXACT_MASS",
                     synonyms = "PUBCHEM_IUPAC_TRADITIONAL_NAME"
                       # Others:
                                        # PUBCHEM_IUPAC_SYSTEMATIC_NAME,
                                        # PUBCHEM_IUPAC_CAS_NAME,
                                        # PUBCHEM_IUPAC_OPENEYE_NAME,
                                        # PUBCHEM_IUPAC_NAME
                     )
.pubchem_separator <- "; " # there seems to be none
.mona_colmap <- c(id = "ID",
                  name = "NAME",
                  inchi = "INCHIKEY",
                  inchi_key = "INCHIKEY",
                  formula = "FORMULA",
                  mass = "EXACT MASS",
                  synonyms = "SYNONYMS")
.mona_separator <- " __ "

#' @description Import compound information from a LipidBlast file in json
#'     format.
#'
#' @note This is a modified version from Jan's generate_lipidblast_tbl that
#'     extracts the mass also from the json and does not calculate it.
#'
#' @author Jan Stanstrup and Johannes Rainer
#'
#' @importFrom jsonlite read_json
#' @importFrom dplyr bind_rows
#'
#' @md
#'
#' @noRd
.import_lipidblast <- function(file) {
    lipidb <- read_json(file)

    parse_element <- function(x) {
        id <- x$id
        cmp <- x$compound[[1]]
        ## get the name(s) -> name + aliases
        nms <- vapply(cmp$names, `[[`, "name", FUN.VALUE = "character")
        mass <- unlist(lapply(cmp$metaData, function(z) {
            if (z$name == "total exact mass")
                z$value
        }))
        if (is.null(mass))
            mass <- NA_character_
        frml <- unlist(lapply(cmp$metaData, function(z) {
            if (z$name == "molecular formula")
                z$value
        }))
        if (is.null(frml))
            mass <- NA_character_
        list(
            compound_id = x$id,
            compound_name = nms[1],
            inchi = cmp$inchi,
            inchi_key = NA_character_,
            formula = frml,
            mass = mass,
            synonyms = nms[-1]
        )
    }

    res <- lapply(lipidb, parse_element)
    bind_rows(res)
}

#' @title Create a CompDb database
#'
#' @description
#'
#' `createCompDb` creates a `SQLite`-based [`CompDb`] object/database
#' from a compound resource provided as a `data.frame` or `tbl`. Alternatively,
#' the name(s) of the file(s) from which the annotation should be extracted can
#' be provided. Supported are SDF files (such as those from the
#' *Human Metabolome Database* HMDB) that can be read using the
#' [compound_tbl_sdf()] or LipidBlast files (see [compound_tbl_lipidblast()].
#'
#' An additional `data.frame` providing metadata information including the data
#' source, date, version and organism is mandatory.
#'
#' Optionally MS/MS (MS2) spectra for compounds can be also stored in the
#' database. Currently only MS/MS spectra from HMDB are supported. These can
#' be downloaded in XML format from HMDB (http://www.hmdb.ca), loaded with
#' the [msms_spectra_hmdb()] or [msms_spectra_mona()] function and passed to
#' the function with the `msms_spectra` argument. See [msms_spectra_hmdb()] or
#' [msms_spectra_mona()] for information on the expected columns and format.
#'
#' Required columns for the `data.frame` providing the compound information (
#' parameter `x`) are:
#' + `"id"`: the ID of the compound (e.g. an HMDB ID).
#' + `"name"`: the compound's name.
#' + `"inchi"`: the inchi of the compound.
#' + `"formula"`: the chemical formula.
#' + `"mass"`: the compound's mass.
#' + `"synonyms"`: additional synonyms/aliases for the compound. Should be
#'   either a single character or a list of values for each compound.
#'
#' See e.g. [compound_tbl_sdf()] or [compound_tbl_lipidblast()] for functions
#' creating such compound tables.
#'
#' The metadata `data.frame` is supposed to have two columns named `"name"` and
#' `"value"` providing the following minimal information as key-value pairs
#' (see `make_metadata` for a unitlity function to create such a `data.frame):
#' + `"source"`: the source from which the data was retrieved (e.g. `"HMDB"`).
#' + `"url"`: the url from which the original data was retrieved.
#' + `"source_version"`: the version from the original data source
#'   (e.g. `"v4"`).
#' + `"source_date"`: the date when the original data source was generated.
#' + `"organism"`: the organism. Should be in the form `"Hsapiens"` or
#'   `"Mmusculus"`.
#'
#' @details
#'
#' Metadata information is also used to create the file name for the database
#' file. The name starts with `"CompDb"`, followed by the organism, the
#' data source and its version. A compound database file for HMDB version 4
#' with human metabolites will thus be named: `"CompDb.Hsapiens.HMDB.v4"`.
#'
#' A single `CompDb` database is created from multiple SDF files (e.g. for
#' *PubChem*) if all the file names are provided with parameter `x`. Parallel
#' processing is currently not enabled because SQLite does not support it yet
#' natively.
#'
#' @param x For `createCompDb`: `data.frame` or `tbl` with the compound
#'     annotations or `character` with the file name(s) from which the compound
#'     annotations should be retrieved. See description for details.
#'
#'     For `createCompDbPackage`: `character(1)` with the file name of the
#'     `CompDb` SQLite file (created by `createCompDb`).
#'
#' @param metadata For `createCompDb`: `data.frame` with metadata
#'     information. See description for details.
#'
#' @param msms_spectra For `createCompDb`: `data.frame` with MS/MS spectrum
#'     data. See [msms_spectra_hmdb()] for the expected format and a function
#'     to import such data from spectrum xml files from HMDB.
#'
#' @param path `character(1)` with the path to the directory where the database
#'     file or package folder should be written. Defaults to the current
#'     directory.
#'
#' @return For `createCompDb`: a `character(1)` with the database name
#'     (invisibly).
#'
#' @importFrom DBI dbDriver dbWriteTable dbExecute dbDisconnect
#' @importFrom RSQLite dbConnect
#' @importFrom dplyr bind_cols bind_rows
#'
#' @export
#'
#' @author Johannes Rainer
#'
#' @seealso
#'
#' [compound_tbl_sdf()] and [compound_tbl_lipidblast()] for functions
#' to extract compound annotations from files in SDF format, or files from
#' LipidBlast.
#'
#' [import_mona_sdf()] to import both the compound and spectrum data from a
#' SDF file from MoNa (Massbank of North America) in one call.
#'
#' [msms_spectra_hmdb()] and [msms_spectra_mona()] for functions to import
#' MS/MS spectrum data from xml files from HMDB or an SDF file from MoNa.
#'
#' [CompDb()] for how to use a compound database.
#'
#' @md
#'
#' @examples
#'
#' ## Read compounds for a HMDB subset
#' fl <- system.file("sdf/HMDB_sub.sdf.gz", package = "CompoundDb")
#' cmps <- compound_tbl_sdf(fl)
#'
#' ## Create a metadata data.frame for the compounds.
#' metad <- data.frame(name = c("source", "url", "source_version",
#'     "source_date", "organism"), value = c("HMDB", "http://www.hmdb.ca",
#'     "v4", "2017-08-27", "Hsapiens"))
#'
#' ## Alternatively use the make_metadata helper function
#' metad <- make_metadata(source = "HMDB", source_version = "v4",
#'     source_date = "2017-08", organism = "Hsapiens",
#'     url = "http://www.hmdb.ca")
#' ## Create a SQLite database in the temporary folder
#' db_f <- createCompDb(cmps, metadata = metad, path = tempdir())
#'
#' ## The database can be loaded and accessed with a CompDb object
#' db <- CompDb(db_f)
#' db
#'
#' ## Create a database for HMDB that includes also MS/MS spectrum data
#' metad2 <- make_metadata(source = "HMDB_with_spectra", source_version = "v4",
#'     source_date = "2017-08", organism = "Hsapiens",
#'     url = "http://www.hmdb.ca")
#'
#' ## Import spectrum information from selected MS/MS xml files from HMDB
#' ## that are provided in the package
#' xml_path <- system.file("xml", package = "CompoundDb")
#' spctra <- msms_spectra_hmdb(xml_path)
#'
#' ## Create a SQLite database in the temporary folder
#' db_f2 <- createCompDb(cmps, metadata = metad2, msms_spectra = spctra,
#'     path = tempdir())
#'
#' ## The database can be loaded and accessed with a CompDb object
#' db2 <- CompDb(db_f2)
#' db2
#'
#' ## Does the database contain MS/MS spectrum data?
#' hasMsMsSpectra(db2)
#'
#' ## Create a database for a ChEBI subset providing the file name of the
#' ## corresponding SDF file
#' metad <- make_metadata(source = "ChEBI_sub", source_version = "2",
#'     source_date = NA, organism = "Hsapiens", url = "www.ebi.ac.uk/chebi")
#' db_f <- createCompDb(system.file("sdf/ChEBI_sub.sdf.gz",
#'     package = "CompoundDb"), metadata = metad, path = tempdir())
#' db <- CompDb(db_f)
#' db
#'
#' compounds(db)
#'
#' ## connect to the database and query it's tables using RSQlite
#' library(RSQLite)
#' con <- dbConnect(dbDriver("SQLite"), db_f)
#'
#' dbGetQuery(con, "select * from metadata")
#' dbGetQuery(con, "select * from compound")
#'
#' ## To create a CompDb R-package we could simply use the
#' ## createCompDbPackage function on the SQLite database file name.
createCompDb <- function(x, metadata, msms_spectra, path = ".") {
    .valid_metadata(metadata)
    db_file <- paste0(path, "/", .db_file_from_metadata(metadata), ".sqlite")
    con <- dbConnect(dbDriver("SQLite"), dbname = db_file)
    ## Add additional metadata info
    metadata$name <- as.character(metadata$name)
    metadata$value <- as.character(metadata$value)
    metadata <- bind_rows(metadata, c(name = "db_creation_date",
                                      value = date()))
    metadata <- bind_rows(metadata, c(name = "supporting_package",
                                      value = "CompoundDb"))
    metadata <- bind_rows(metadata, c(name = "supporting_object",
                                      value = "CompDb"))
    dbWriteTable(con, name = "metadata", metadata, row.names = FALSE)
    ## compound table
    if (is.character(x)) {
        lapply(x, function(z) {
            message("Import data from ", basename(z), " ...", appendLF = FALSE)
            if (.is_sdf_filename(z)) {
                tbl <- compound_tbl_sdf(z)
            } else {
                ## Assume this is a json file... might simply be wrong.
                tbl <- compound_tbl_lipidblast(z)
            }
            message("OK")
            .valid_compound(tbl, db = FALSE)
            tbl_syn <- bind_cols(compound_id = rep(tbl$compound_id,
                                                   lengths(tbl$synonyms)),
                                 synonym = unlist(tbl$synonyms))
            tbl <- tbl[, colnames(tbl) != "synonyms"]
            dbWriteTable(con, name = "synonym", tbl_syn, row.names = FALSE,
                         append = TRUE)
            dbWriteTable(con, name = "compound", row.names = FALSE,
                         tbl[, colnames(tbl) != "synonyms"], append = TRUE)
        })
    } else {
        .valid_compound(x, db = FALSE)
        x_synonym <- bind_cols(compound_id = rep(x$compound_id,
                                                 lengths(x$synonyms)),
                               synonym = unlist(x$synonyms))
        dbWriteTable(con, name = "synonym", x_synonym, row.names = FALSE)
        dbWriteTable(con, name = "compound", x[, colnames(x) != "synonyms"],
                     row.names = FALSE)
    }
    ## Creating indices
    dbExecute(con, "create index compound_id_idx on compound (compound_id)")
    dbExecute(con, "create index compound_name_idx on compound (compound_name)")
    ## Process spectra.
    if (!missing(msms_spectra) && is.data.frame(msms_spectra)) {
        comp_ids <- unique(x$compound_id)
        if (!all(msms_spectra$compound_id %in% comp_ids))
            stop("All compound identifiers in 'msms_spectra' have to be ",
                 "present also in 'x'")
        msms_spectra <- msms_spectra[msms_spectra$compound_id %in% comp_ids, ]
        .insert_msms_spectra_blob(con, msms_spectra)
    }
    dbDisconnect(con)
    invisible(db_file)
}

#' Function to insert MS/MS spectrum data into a single database table with
#' one row per spectrum and m/z and intensity vectors stored as BLOB.
#'
#' @param con database connection
#'
#' @param x `data.frame` with the spectrum data, such as returned by
#' [msms_spectra_hmdb()]
#'
#' @noRd
#'
#' @author Johannes Rainer
.insert_msms_spectra_blob <- function(con, x, ...) {
    if (is.numeric(x$mz))
        x <- .collapse_spectrum_df(x)
    x <- .add_mz_range_column(x)
    x$collision_energy <- as.character(x$collision_energy)
    .valid_msms_spectrum(x, blob = TRUE)
    x$mz <- lapply(x$mz, base::serialize, NULL)
    x$intensity <- lapply(x$intensity, base::serialize, NULL)
    dbWriteTable(con, name = "msms_spectrum", x, row.names = FALSE, ...)
    dbExecute(con, "create index msms_cidb_idx on msms_spectrum (compound_id)")
}

#' Add mz_min and mz_max columns to data.frame x
#'
#' @param x `data.frame` with MS MS spectrum data such as returned by
#'     [msms_spectra_hmdb()] **and** eventually collapsed with
#'     [.collapse_spectrum_df].
#'
#' @noRd
#'
#' @author Johannes Rainer
.add_mz_range_column <- function(x) {
    if (!any(colnames(x) == "mz"))
        stop("Required column 'mz' not found")
    if (is.list(x$mz)) {
        mzr <- do.call(rbind, lapply(x$mz, range))
        x$msms_mz_range_min <- mzr[, 1]
        x$msms_mz_range_max <- mzr[, 2]
        x
    } else {
        fct <- factor(x$spectrum_id, levels = unique(x$spectrum_id))
        mzs <- split(x$mz, f = fct)
        mzr <- do.call(rbind, lapply(mzs, range, na.rm = TRUE))
        x$msms_mz_range_min <- unsplit(mzr[, 1], fct)
        x$msms_mz_range_max <- unsplit(mzr[, 2], fct)
        x
    }
}

#' Function to insert MS/MS spectrum data into two database tables:
#' msms_spectrum_peak with the individual m/z and intensity values and
#' msms_spectrum_metadata with the spectrum annotations. Both tables are linked
#' via the spectrum_id.
#'
#' @param con database connection
#'
#' @param x `data.frame` with the spectrum data, such as returned by
#' [msms_spectra_hmdb]
#'
#' @noRd
#'
#' @author Johannes Rainer
.insert_msms_spectra <- function(con, x) {
    if (is.list(x$mz))
        x <- .expand_spectrum_df(x)
    x$collision_energy <- as.character(x$collision_energy)
    .valid_msms_spectrum(x, blob = FALSE)
    msms_spectrum_peak <- x[, c("spectrum_id", "mz", "intensity")]
    dbWriteTable(con, name = "msms_spectrum_peak", msms_spectrum_peak,
                 row.names = FALSE)
    x <- .add_mz_range_column(x)
    msms_spectrum_metadata <- unique(x[, !(colnames(x) %in%
                                           c("mz", "intensity"))])
    dbWriteTable(con, name = "msms_spectrum_metadata",
                 msms_spectrum_metadata, row.names = FALSE)
    dbExecute(con, "create index msms_id_idx on msms_spectrum_peak (spectrum_id)")
    dbExecute(con, "create index msms_mid_idx on msms_spectrum_metadata (spectrum_id)")
    dbExecute(con, "create index msms_cid_idx on msms_spectrum_metadata (compound_id)")
}


.required_metadata_keys <- c("source", "url", "source_version", "source_date",
                             "organism")
.required_compound_db_columns <- c("compound_id", "compound_name", "inchi",
                                   "inchi_key", "formula", "mass")
.required_compound_columns <- c(.required_compound_db_columns, "synonyms")

.required_msms_spectrum_columns <- c(spectrum_id = "character",
                                     compound_id = "character",
                                     polarity = "integer",
                                     collision_energy = "character",
                                     predicted = "logical",
                                     splash = "character",
                                     instrument_type = "character",
                                     instrument = "character",
                                     precursor_mz = "numeric",
                                     mz = "numeric",
                                     intensity = "numeric"
                                     )

.required_msms_spectrum_columns_blob <- c(spectrum_id = "character",
                                          compound_id = "character",
                                          polarity = "integer",
                                          collision_energy = "character",
                                          predicted = "logical",
                                          splash = "character",
                                          instrument_type = "character",
                                          instrument = "character",
                                          precursor_mz = "numeric",
                                          mz = "list",
                                          intensity = "list"
                                          )

.valid_msms_spectrum <- function(x, blob = TRUE, error = TRUE) {
    coltypes <- unlist(lapply(x, class))
    msg <- NULL
    if (blob)
        req_cols <- .required_msms_spectrum_columns_blob
    else req_cols <- .required_msms_spectrum_columns
    not_found <- which(!(names(req_cols) %in% names(coltypes)))
    if (length(not_found))
        msg <- c(msg, paste0("Required column(s) ",
                             paste0("'", names(req_cols)[not_found], "'"),
                             " not found in 'msms_spectra'"))
    common_cols <- union(names(coltypes), names(req_cols))
    wrong_type <- which(coltypes[common_cols] != req_cols[common_cols])
    if (length(wrong_type))
        msg <- c(msg,
                 paste0("One or more columns contain data of the wrong type: ",
                        paste(names(coltypes[names(wrong_type)]),
                              coltypes[names(wrong_type)], sep = ":",
                              collapse = ", "),
                        ". Please refer to the documentation of createCompDb ",
                        "for the correct data types."))
    if (length(msg))
        if (error) stop(msg)
        else msg
    else TRUE
}

#' @description Create the database file name from the metadata `data.frame`.
#'     The function checks also for the presence of all required metadata fields
#'     ensuring that these are also not `NA` or `NULL`.
#'
#' @md
#'
#' @noRd
.db_file_from_metadata <- function(x) {
    paste0("CompDb.", x$value[x$name == "organism"], ".",
           x$value[x$name == "source"], ".",
           x$value[x$name == "source_version"])
}

#' @description Check the metadata `data.frame` for required columns.
#'
#' @md
#'
#' @noRd
.valid_metadata <- function(metadata, error = TRUE) {
    txt <- character()
    if (!is.data.frame(metadata))
        txt <- c(txt, "'metadata' is expected to be a data.frame")
    if (all(c("name", "value") %in% colnames(metadata))) {
        keys <- metadata$name
        vals <- metadata$value
        got_it <- .required_metadata_keys %in% keys
        if (!all(got_it))
            txt <- c(txt, paste0("required fields ",
                                 paste0(.required_metadata_keys[!got_it],
                                        collapse = ", "),
                                 " not found in metadata data.frame"))
        vals <- vals[keys %in% .required_metadata_keys]
    } else {
        txt <- c(txt, paste0("metadata data.frame needs to have columns ",
                             "named 'name' and 'value'"))
    }
    if (length(txt))
        if (error)
            stop(paste(txt, collapse = "\n"))
        else txt
    else TRUE
}

#' @description Check that the compounds table contains all required data.
#'
#' @param db `logical(1)` whether validity should be checked on the internal
#'     database table instead of the input file.
#' @md
#'
#' @noRd
.valid_compound <- function(x, error = TRUE, db = TRUE) {
    txt <- character()
    if (!is.data.frame(x))
        txt <- c(txt, "'x' is supposed to be a data.frame")
    .req_cols <- .required_compound_db_columns
    if (!db)
        .req_cols <- .required_compound_columns
    got_it <- .req_cols %in% colnames(x)
    if (!all(got_it)) {
        txt <- c(txt, paste0("Miss required columns: ",
                             paste0(.required_compound_columns[!got_it],
                                    collapse = ", ")))
    } else {
        if (!is.numeric(x$mass))
            txt <- c(txt, "Column 'mass' should be numeric")
    }
    if (db) {
        ## Do not allow more columns than expected!
        is_ok <- colnames(x) %in% .req_cols
        if (any(!is_ok)) {
            txt <- c(txt, paste0("Column(s) ", paste(colnames(x)[!is_ok]),
                                 " are not expected"))
        }
    }
    if (length(txt))
        if (error)
            stop(paste(txt, collapse = "\n"))
        else txt
    else TRUE
}

#' @description
#'
#' `createCompDbPackage` creates an R data package with the data from a
#' [`CompDb`] object.
#'
#' @importFrom Biobase createPackage
#'
#' @param version For `createCompDbPackage`: `character(1)` with the version
#'     of the package (ideally in the format `"x.y.z"`).
#'
#' @param maintainer For `createCompDbPackage`: `character(1)` with the
#'     name and email address of the package maintainer (in the form
#'     `"First Last <first.last@provider.com>"`.
#'
#' @param author For `createCompDbPackage`: `character(1)` with the name
#'     of the package author.
#'
#' @param license For `createCompDbPackage`: `character(1)` with the
#'     license of the package respectively the originating provider.
#'
#' @export
#'
#' @md
#'
#' @rdname createCompDb
createCompDbPackage <- function(x, version, maintainer, author,
                                    path = ".", license = "Artistic-2.0") {
    if (missing(x) | missing(version) | missing(maintainer) | missing(author))
        stop("'x', 'version', 'maintainer' and 'author' are required")
    if (!is.character(x))
        stop("'x' is supposed to be the file name of the CompDb")
    cdb <- CompDb(x)
    metad <- .metadata(cdb)
    pkg_name <- .db_file_from_metadata(metad)
    m_source <- .metadata_value(cdb, "source")
    m_source_version <- .metadata_value(cdb, "source_version")
    m_source_date <- .metadata_value(cdb, "source_date")
    m_organism <- .metadata_value(cdb, "organism")
    m_source_url <- .metadata_value(cdb, "url")
    template_path <- system.file("pkg-template", package = "CompoundDb")
    symvals <- list(
        PKGTITLE = paste0(m_source, " compound annotation package"),
        PKGDESCRIPTION = paste0("Exposes a CompDb compound annotation ",
                                "databases with annotations retrieved from ",
                                m_source, "."),
        PKGVERSION  = version,
        AUTHOR = author,
        MAINTAINER = maintainer,
        LIC = license,
        ORGANISM = m_organism,
        SPECIES = m_organism,
        SOURCE = m_source,
        SOURCEVERSION = as.character(m_source_version),
        RELEASEDATE = m_source_date,
        SOURCEURL =  m_source_url,
        DBOBJNAME = pkg_name
    )
    createPackage(pkgname = pkg_name, destinationDir = path,
                  originDir = template_path, symbolValues = symvals)
    sqlite_path <- file.path(path, pkg_name, "inst", "extdata")
    dir.create(sqlite_path, showWarnings = FALSE, recursive = TRUE)
    file.copy(x, to = file.path(sqlite_path, paste0(pkg_name, ".sqlite")))
    invisible(file.path(path, pkg_name))
}

#' @description Simple test function to evaluate whether a file is an SDF file
#'     based on its file name.
#'
#' @param x `character(1)` the file name to be tested.
#'
#' @noRd
#'
#' @md
#'
#' @author Johannes Rainer
.is_sdf_filename <- function(x) {
    grepl(".sdf($|.gz$)", x, ignore.case = TRUE)
}

#' @description `make_metadata` helps generating a metadata `data.frame` in the
#'     correct format expected by the `createCompDb` function. The function
#'     returns a `data.frame`.
#'
#' @param source For `make_metadata`: `character(1)` with the name of the
#'     resource that provided the compound annotation.
#'
#' @param url For `make_metadata`: `character(1)` with the url to the original
#'     resource.
#'
#' @param source_version For `make_metadata`: `character(1)` with the version
#'     of the original resource providing the annotation.
#'
#' @param source_date For `make_metadata`: `character(1)` with the date of the
#'     resource's release.
#'
#' @param organism For `make_metadata`: `character(1)` with the name of the
#'     organism. This should be in the format `"Hsapiens"` for human,
#'     `"Mmusculus"` for mouse etc.
#'
#' @export
#'
#' @md
#'
#' @rdname createCompDb
make_metadata <- function(source, url, source_version, source_date, organism) {
    if (any(c(missing(source), missing(url), missing(source_version),
              missing(source_date), missing(organism))))
        stop("Arguments 'source', 'url', 'source_version', 'source_date', ",
             "'organism' are required")
    if (any(c(is.null(source), is.null(url), is.null(source_version),
              is.null(source_date), is.null(organism))))
        stop("'NULL' is not allowed for 'source', 'url', 'source_version', ",
             "'organism'")
    data.frame(name = c("source", "url", "source_version", "source_date",
                        "organism"),
               value = as.character(c(source, url, source_version,
                                      source_date, organism)),
               stringsAsFactors = FALSE)
}

#' @title Import compound and spectrum information from MoNa
#'
#' `import_mona_sdf` allows to import compound and spectrum information from an
#' SDF file from MoNa (Massbank of North America
#' http://mona.fiehnlab.ucdavis.edu/). This function is a convenience function
#' using the [compound_tbl_sdf()] and [msms_spectra_mona()] functions for data
#' import but avoiding to read the SDF files twice.
#'
#' @note
#'
#' MoNa SDF files organize the data by individual spectra (i.e. each element
#' is one spectrum) and individual compounds can not easily and consistently
#' defined (i.e. not all entries have an InChI ID or other means to uniquely
#' identify compounds). Thus, the function returns a highly redundant compount
#' table. Feedback on how to reduce this redundancy would be highly welcome!
#'
#' @param x `character(1)` being the SDF file name.
#'
#' @param nonStop `logical(1)` wheter file content specific errors should
#'     only reported as warnings and not break the full import process.
#'
#' @author Johannes Rainer
#'
#' @return A `list` with elements `"compound"` and `"msms_spectrum"` containing
#'     data.frames with compound and MS/MS spectrum data, respectively.
#'
#' @export
#'
#' @md
#'
#' @seealso
#'
#' [compound_tbl_sdf()] to read only the compound information.
#'
#' [msms_spectra_mona()] to read only the spectrum data.
#'
#' @examples
#'
#' ## Define the test file containing a small subset from MoNa
#' fl <- system.file("sdf/MoNa_export-All_Spectra_sub.sdf.gz",
#'     package = "CompoundDb")
#'
#' ## Import the data
#' res <- import_mona_sdf(fl)
import_mona_sdf <- function(x, nonStop = TRUE) {
    message("Reading SDF file ... ", appendLF = FALSE)
    sdfs <- datablock2ma(datablock(read.SDFset(x, skipErrors = nonStop)))
    if (!any(colnames(sdfs) == "MASS SPECTRAL PEAKS"))
        stop("The provided file does not contain \"MASS SPECTRAL PEAKS\" ",
             "elements. Is this an SDF file from MoNa?")
    message("OK")
    message("Extracting compound information ... ", appendLF = FALSE)
    suppressMessages(cmps <- .simple_extract_compounds_sdf(sdfs))
    message("OK")
    message("Extracting spectrum information ... ", appendLF = FALSE)
    spctra <- .extract_spectra_mona_sdf(sdfs)
    message("OK")
    list(compound = cmps, msms_spectrum = spctra)
}
