#' @title Create a CompoundDb database
#'
#' @description
#' 
#' `createCompoundDb` creates a `SQLite`-based [`CompoundDb`] object/database
#' from a compound resource provided as a `data.frame` or `tbl`.
#' An additional `data.frame` providing metadata information is mandatory.
#' Required columns for the `data.frame` providing the compound information are:
#' + `"id"`: the ID of the compound (e.g. an HMDB ID).
#' + `"name"`: the compound's name.
#' + `"inchi"`: the inchi of the compound.
#' + `"formula"`: the chemical formula.
#' + `"mass"`: the compound's mass.
#'
#' See e.g. [generate_hmdb_tbl()] or [generate_lipidblast_tbl()] for functions
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
#' file. The name starts with `"CompoundDb"`, followed by the organism, the
#' data source and its version. A compound database file for HMDB version 4
#' with human metabolites will thus be named: `"CompoundDb.Hsapiens.HMDB.v4"`.
#' 
#' @param x For `createCompoundDb`: `data.frame` or `tbl` with the compound
#'     annotations. See description for details.
#'
#'     For `createCompoundDbPackage`: `character(1)` with the file name of the
#'     `CompoundDb` SQLite file (created by `createCompoundDb`).
#'
#' @param metadata For `createCompoundDb`: `data.frame` with metadata
#'     information. See description for details.
#'
#' @param path `character(1)` with the path to the directory where the database
#'     file or package folder should be written. Defaults to the current
#'     directory.
#'
#' @return For `createCompoundDb`: a `character(1)` with the database name
#'     (invisibly).
#' 
#' @importFrom DBI dbDriver dbWriteTable dbExecute dbDisconnect
#' @importFrom RSQLite dbConnect
#'
#' @export
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' ## Read compounds for a HMDB subset
#' fl <- system.file("extdata/hmdb/hmdb_sub.xml", package = "PeakABro")
#' cmps <- generate_hmdb_tbl(fl)
#'
#' ## Create a metadata data.frame for the compounds.
#' metad <- data.frame(name = c("source", "url", "source_version",
#'     "source_date", "organism"), value = c("HMDB", "http://www.hmdb.ca",
#'     "v4", "2017-08-27", "Hsapiens"))
#'
#' ## Create a SQLite database in the temporary folder
#' db_f <- createCompoundDb(cmps, metadata = metad, path = tempdir())
#'
#' ## connect to the database and query it's tables using RSQlite
#' library(RSQLite)
#' con <- dbConnect(dbDriver("SQLite"), db_f)
#'
#' dbGetQuery(con, "select * from metadata")
#' dbGetQuery(con, "select * from compound")
#'
#' ## To create a CompoundDb R-package we could simply use the
#' ## createCompoundDb function on the SQLite database file name.
createCompoundDb <- function(x, metadata, path = ".") {
    .valid_metadata(metadata)
    .valid_compound(x)
    db_file <- paste0(path, "/", .db_file_from_metadata(metadata), ".sqlite")
    con <- dbConnect(dbDriver("SQLite"), dbname = db_file)
    ## Add additional metadata info
    metadata <- rbind(metadata, c("db_creation_date", date()))
    metadata <- rbind(metadata, c("supporting_package", "PeakABro"))
    metadata <- rbind(metadata, c("supporting_object", "CompoundDb"))
    dbWriteTable(con, name = "metadata", metadata, row.names = FALSE)
    dbWriteTable(con, name = "compound", x, row.names = FALSE)
    ## Creating indices
    dbExecute(con, paste0("create index compound_id_idx on compound (id)"))
    dbExecute(con, paste0("create index compound_name_idx on compound (name)"))
    dbDisconnect(con)
    invisible(db_file)
}

.required_metadata_keys <- c("source", "url", "source_version", "source_date",
                             "organism")
.required_compound_columns <- c("id", "name", "inchi", "formula", "mass")

#' @description Create the database file name from the metadata `data.frame`.
#'     The function checks also for the presence of all required metadata fields
#'     ensuring that these are also not `NA` or `NULL`.
#'
#' @md
#' 
#' @noRd
.db_file_from_metadata <- function(x) {
    paste0("CompoundDb.", x$value[x$name == "organism"], ".",
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
#' `createCompoundDbPackage` creates an R data package with the data from a
#' [`CompoundDb`] object.
#'
#' @importFrom Biobase createPackage
#'
#' @param version For `createCompoundDbPackage`: `character(1)` with the version
#'     of the package (ideally in the format `"x.y.z"`).
#'
#' @param maintainer For `createCompoundDbPackage`: `character(1)` with the
#'     name and email address of the package maintainer (in the form
#'     `"First Last <first.last@provider.com>"`.
#'
#' @param author For `createCompoundDbPackage`: `character(1)` with the name
#'     of the package author.
#'
#' @param license For `createCompoundDbPackage`: `character(1)` with the
#'     license of the package respectively the originating provider.
#' 
#' @export
#'
#' @md
#' 
#' @rdname createCompoundDb
createCompoundDbPackage <- function(x, version, maintainer, author,
                                    path = ".", license = "Artistic-2.0") {
    if (missing(x) | missing(version) | missing(maintainer) | missing(author))
        stop("'x', 'version', 'maintainer' and 'author' are required")
    if (!is.character(x))
        stop("'x' is supposed to be the file name of the CompoundDb")
    cdb <- CompoundDb(x)
    metad <- .metadata(cdb)
    pkg_name <- .db_file_from_metadata(metad)
    m_source <- .metadata_value(cdb, "source")
    m_source_version <- .metadata_value(cdb, "source_version")
    m_source_date <- .metadata_value(cdb, "source_date")
    m_organism <- .metadata_value(cdb, "organism")
    m_source_url <- .metadata_value(cdb, "url")
    template_path <- system.file("pkg-template", package = "PeakABro")
    symvals <- list(
        PKGTITLE = paste0(m_source, " compound annotation package"),
        PKGDESCRIPTION = paste0("Exposes a CompoundDb compound annotation ",
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

#' @description
#'
#' Internal function to extract compound information from a file in SDF format.
#'
#' @param x `character(1)` with the name of the file.
#'
#' @return A `data.frame` with columns `"compound_id"`, `"compound_name"`,
#'     `"inchi"`, `"formula"`, `"mass"`.
#' 
#' @importFrom ChemmineR read.SDFset datablock datablock2ma
#' 
#' @md
#'
#' @author Johannes Rainer
#' 
#' @noRd
.simple_import_compounds_sdf <- function(x) {
    dblock <- datablock(read.SDFset(x))
    full_mat <- datablock2ma(dblock)
    data.frame(compound_id = full_mat[, "DATABASE_ID"],
               compound_name = full_mat[, "GENERIC_NAME"],
               inchi = full_mat[, "INCHI_IDENTIFIER"],
               formula = full_mat[, "FORMULA"],
               mass = as.numeric(full_mat[, "EXACT_MASS"]),
               stringsAsFactors = FALSE)
}
