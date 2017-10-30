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
#'
#' @details
#'
#' Column `"compound_name"` reports for HMDB files the `"GENERIC_NAME"`, for
#' ChEBI the `"ChEBI Name"`, for PubChem the `"PUBCHEM_IUPAC_TRADITIONAL_NAME"`,
#' and for Lipid Maps the `"COMMON_NAME"`, if that is
#' not available, the first of the compounds synonyms and, if that is also not
#' provided, the `"SYSTEMATIC_NAME"`.
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
#' @seealso [createCompDb()] for a function to create a SQLite-based compound
#'     database.
#' 
#' @examples
#'
#' ## Read compound information from a subset of HMDB
#' fl <- system.file("sdf/HMDB_sub.sdf", package = "CompoundDb")
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
    res <- .simple_import_compounds_sdf(file)
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
#' @param x `character(1)` with the name of the file.
#'
#' @return A [tibble::tibble] with columns `"compound_id"`, `"compound_name"`,
#'     `"inchi"`, `"formula"`, `"mass"`.
#' 
#' @importFrom ChemmineR read.SDFset datablock datablock2ma
#' @importFrom tibble data_frame
#' 
#' @md
#'
#' @author Johannes Rainer
#' 
#' @noRd
.simple_import_compounds_sdf <- function(x) {
    full_mat <- datablock2ma(datablock(read.SDFset(x)))
    source_db <- .guess_sdf_source(colnames(full_mat))
    if (is.null(source_db))
        stop("The SDF file is not supported. Supported are SDF files from ",
             "HMDB, ChEBI and LipidMaps.")
    colmap <- get(paste0(".", source_db, "_colmap"))
    sep <- get(paste0(".", source_db, "_separator"))
    ## Fix missing COMMON_NAME entries in LipidMaps (see issue #1)
    nms <- full_mat[, colmap["name"]]
    syns <- strsplit(full_mat[, colmap["synonyms"]], split = sep)
    if (source_db == "lipidmaps") {
        nas <- is.na(nms)
        if (any(nas))
            nms[nas] <- vapply(syns[nas], `[[`, 1, FUN.VALUE = "character")
        nas <- is.na(nms)
        if (any(nas))
            nms[nas] <- full_mat[nas, "SYSTEMATIC_NAME"]
    }
    data_frame(compound_id = full_mat[, colmap["id"]],
               compound_name = nms,
               inchi = full_mat[, colmap["inchi"]],
               formula = full_mat[, colmap["formula"]],
               mass = as.numeric(full_mat[, colmap["mass"]]),
               synonyms = syns
               )
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
    
    NULL
}


.hmdb_colmap <- c(id = "HMDB_ID",
                  name = "GENERIC_NAME",
                  inchi = "INCHI_IDENTIFIER",
                  formula = "FORMULA",
                  mass = "EXACT_MASS",
                  synonyms = "SYNONYMS"
                  )
.hmdb_separator <- "; "
.chebi_colmap <- c(id = "ChEBI ID",
                   name = "ChEBI Name",
                   inchi = "InChI",
                   formula = "Formulae",
                   mass = "Monoisotopic Mass",
                   synonyms = "Synonyms"
                   )
.chebi_separator <- " __ "
.lipidmaps_colmap <- c(id = "LM_ID",
                       name = "COMMON_NAME",
                       inchi = "INCHI",
                       formula = "FORMULA",
                       mass = "EXACT_MASS",
                       synonyms = "SYNONYMS"
                       )
.lipidmaps_separator <- "; "
.pubchem_colmap <- c(id = "PUBCHEM_COMPOUND_CID",
                       name = "PUBCHEM_IUPAC_TRADITIONAL_NAME",
                       inchi = "PUBCHEM_IUPAC_INCHI",
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
#' the name of the file from which the annotation should be extracted can be
#' provided. Supported are e.g. SDF files that can be processed using the
#' [compound_tbl_sdf()] or LipidBlast files (see [compound_tbl_lipidblast()].
#' 
#' An additional `data.frame` providing metadata information is mandatory.
#' Required columns for the `data.frame` providing the compound information are:
#' + `"id"`: the ID of the compound (e.g. an HMDB ID).
#' + `"name"`: the compound's name.
#' + `"inchi"`: the inchi of the compound.
#' + `"formula"`: the chemical formula.
#' + `"mass"`: the compound's mass.
#'
#' See e.g. [compound_tbl_sdf()] or [compound_tbl_lipidblast()] for functions
#' creating such compound tables.
#' 
#' The metadata `data.frame` is supposed to have two columns named `"name"` and
#' `"value"` providing the following minimal information as key-value pairs:
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
#' @param x For `createCompDb`: `data.frame` or `tbl` with the compound
#'     annotations. See description for details.
#'
#'     For `createCompDbPackage`: `character(1)` with the file name of the
#'     `CompDb` SQLite file (created by `createCompDb`).
#'
#' @param metadata For `createCompDb`: `data.frame` with metadata
#'     information. See description for details.
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
#' @seealso [compound_tbl_sdf()] and [compound_tbl_lipidblast()] for functions
#'     to extract compound annotations from files in SDF format, or files from
#'     LipidBlast.
#'     [CompDb()] for how to use a compound database.
#' 
#' @md
#'
#' @examples
#'
#' ## Read compounds for a HMDB subset
#' fl <- system.file("sdf/HMDB_sub.sdf", package = "CompoundDb")
#' cmps <- compound_tbl_sdf(fl)
#'
#' ## Create a metadata data.frame for the compounds.
#' metad <- data.frame(name = c("source", "url", "source_version",
#'     "source_date", "organism"), value = c("HMDB", "http://www.hmdb.ca",
#'     "v4", "2017-08-27", "Hsapiens"))
#'
#' ## Create a SQLite database in the temporary folder
#' db_f <- createCompDb(cmps, metadata = metad, path = tempdir())
#'
#' ## The database can be loaded and accessed with a CompDb object
#' db <- CompDb(db_f)
#' db
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
createCompDb <- function(x, metadata, path = ".") {
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
        ## Process each file iteratively.
        ## 1) Determine file type, SDF, json.
        ## 2) write to the table.
    } else {
        .valid_compound(x)
        x_synonym <- bind_cols(compound_id = rep(x$compound_id,
                                                 lengths(x$synonyms)),
                               synonym = unlist(x$synonyms))
        x[,] <- x[, colnames(x) != "synonyms"]
        dbWriteTable(con, name = "synonym", x_synonym, row.names = FALSE)
        dbWriteTable(con, name = "compound", x, row.names = FALSE)
    }
    ## Creating indices
    dbExecute(con, "create index compound_id_idx on compound (compound_id)")
    dbExecute(con, "create index compound_name_idx on compound (compound_name)")
    dbDisconnect(con)
    invisible(db_file)
}

.required_metadata_keys <- c("source", "url", "source_version", "source_date",
                             "organism")
.required_compound_columns <- c("compound_id", "compound_name", "inchi",
                                "formula", "mass", "synonyms")

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
        if (length(vals))
            if (any(is.na(vals)) | any(length(vals) == 0))
                txt <- c(txt, paste0("values for metadata data.frame fields ",
                                     "should not be empty or NA"))
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
#' @md
#' 
#' @noRd
.valid_compound <- function(x, error = TRUE) {
    txt <- character()
    if (!is.data.frame(x))
        txt <- c(txt, "'x' is supposed to be a data.frame")
    got_it <- .required_compound_columns %in% colnames(x)
    if (!all(got_it)) {
        txt <- c(txt, paste0("Miss required columns: ",
                             paste0(.required_compound_columns[!got_it],
                                    collapse = ", ")))
    } else {
        if (!is.numeric(x$mass))
            txt <- c(txt, "Column 'mass' should be numeric")
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
    invisible(TRUE)
}
