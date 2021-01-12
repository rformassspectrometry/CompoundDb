#' @include createCompDbPackage.R

#' @name CompDb
#'
#' @title Simple compound (metabolite) databases
#'
#' @aliases CompDb-class show dbconn,CompDb-method show,CompDb-method compoundVariables
#'
#' @description
#'
#' `CompDb` objects provide access to general (metabolite) compound
#' annotations along with *metadata* information such as the annotation's
#' source, date and release version. The data is stored internally in a
#' database (usually an SQLite database).
#'
#' @details
#'
#' `CompDb` objects should be created using the constructor function
#' `CompDb` providing the name of the (SQLite) database file providing
#' the compound annotation data.
#'
#' @section Retrieve annotations from the database:
#'
#' Annotations/compound informations can be retrieved from a `CompDb` database
#' with the `compounds` and `Spectra` functions:
#'
#' - `compounds` extracts compound data from the `CompDb` object. In contrast
#'   to `src_compdb` it returns the actual data as a `data.frame` (if
#'   `return.type = "data.frame"`) or a [tibble::tibble()] (if
#'   `return.type = "tibble"`). A `compounds` call will always return all
#'   elements from the *ms_compound* table (unless a `filter` is used).
#'
#' - `Spectra` extract spectra from the database and returns them as a
#'   [Spectra()] object from the `Spectra` package. Additional annotations
#'   requested with the `columns` parameter are added as additional spectra
#'   variables.
#'
#' @section General functions:
#'
#' - `CompDb`: connect to a compound database.
#'
#' - `compoundVariables`: returns all available columns/database fields for
#'   compounds.
#'
#' - `dbconn`: returns the connection (of type `DBIConnection`) to the database.
#'
#' - `metadata`: returns general meta data of the compound database.
#'
#' - `spectraVariables`: returns all spectra variables (i.e. columns) available
#'   in the `CompDb`.
#'
#' - `src_compdb` provides access to the `CompDb`'s database *via*
#'   the functionality from the `dplyr`/`dbplyr` package.
#'
#' - `supportedFilters`: provides an overview of the filters that can be
#'   applied on a `CompDb` object to extract only specific data from the
#'   database.
#'
#' - `tables`: returns a named `list` (names being table names) with
#'   the fields/columns from each table in the database.
#'
#'
#' @section Filtering the database:
#'
#' Data access methods such as `compounds` and `Spectra` allow to filter the
#' results using specific filter classes and expressions. Filtering uses the
#' concepts from Bioconductor's `AnnotationFilter` package. All information
#' for a certain compound with the ID `"HMDB0000001"` can for example be
#' retrieved by passing the filter expression
#' `filter = ~ compound_id == "HMDB0000001"` to the `compounds` function.
#'
#' Use the [supportedFilters] function on the [CompDb] object to get a list of
#' all supported filters. See also examples below or the usage vignette for
#' details.
#'
#' @param columns For `compounds`, `Spectra`: `character` with the names of the
#'     database columns that should be retrieved. Use `compoundVariables` and/or
#'     `spectraVariables` for a list of available column names.
#'
#' @param filter For `compounds` and `Spectra`: filter expression or
#'     [AnnotationFilter()] defining a filter to be used to retrieve specific
#'     elements from the database.
#'
#' @param flags flags passed to the SQLite database connection.
#'     See [SQLite()]. Defaults to read-only, i.e. `RSQLite::SQLITE_RO`.
#'
#' @param includeId for `compoundVariables`: `logical(1)` whether the comound
#'     ID (column `"compound_id"`) should be included in the result. The
#'     default is `includeIds = FALSE`.
#'
#' @param object For all methods: a `CompDb` object.
#'
#' @param return.type For `compounds`: either `"data.frame"` or `"tibble"` to
#'     return the result as a [data.frame()] or [tibble()], respectively.
#'
#' @param x For `CompDb`: `character(1)` with the file name of the SQLite
#'     compound database. Alternatively it is possible to provide the connection
#'     to the database with parameter `x`.
#'
#'     For all other methods: a `CompDb` object.
#'
#' @param ... additional arguments. Currently not used.
#'
#' @author Johannes Rainer
#'
#' @seealso
#'
#' [createCompDb()] for the function to create a SQLite compound database.
#'
#' [CompoundIdFilter()] for filters that can be used on the `CompDb` database.
#'
#' @examples
#'
#' ## We load a small compound test database based on MassBank which is
#' ## distributed with this package.
#' cdb <- CompDb(system.file("sql/CompDb.MassBank.sql", package = "CompoundDb"))
#' cdb
#'
#' ## Get general metadata information from the database, such as originating
#' ## source and version:
#' metadata(cdb)
#'
#' ## List all available compound annotations/fields
#' compoundVariables(cdb)
#'
#' ## Extract a data.frame with these annotations for all compounds
#' compounds(cdb)
#'
#' ## Add also the synonyms (aliases) for the compounds. This will cause the
#' ## tables compound and synonym to be joined. The elements of the compound_id
#' ## and name are now no longer unique
#' res <- compounds(cdb, columns = c("name", "synonym"))
#' head(res)
#'
#' ## List all database tables and their columns
#' tables(cdb)
#'
#' ## Any of these columns can be used in the `compounds` call to retrieve
#' ## the specific annotations. The corresponding database tables will then be
#' ## joined together
#' compounds(cdb, columns = c("formula", "publication"))
#'
#' ## Create a Spectra object with all MS/MS spectra from the database.
#' sps <- Spectra(cdb)
#' sps
#'
#' ## Extract spectra for a specific compound.
#' sps <- Spectra(cdb, filter = ~ name == "Mellein")
#' sps
#'
#' ## List all available annotations for MS/MS spectra
#' spectraVariables(sps)
#'
#' ## Get access to the m/z values of these
#' mz(sps)
#'
#' library(Spectra)
#' ## Plot the first spectrum
#' plotSpectra(sps[1])
#'
#'
#' #########
#' ## Filtering the database
#' ##
#' ## Get all compounds with an exact mass between 310 and 320
#' res <- compounds(cdb, filter = ~ exactmass > 310 & exactmass < 320)
#' res
#'
#' ## Get all compounds that have an H14 in their formula.
#' res <- compounds(cdb, filter = FormulaFilter("H14", "contains"))
#' res
#'
#' #########
#' ## Using CompDb with the *tidyverse*
#' ##
#' ## Using return.type = "tibble" the result will be returned as a "tibble"
#' compounds(cdb, return.type = "tibble")
#'
#' ## Use the CompDb in a dplyr setup
#' library(dplyr)
#' src_cmp <- src_compdb(cdb)
#' src_cmp
#'
#' ## Get a tbl for the ms_compound table
#' cmp_tbl <- tbl(src_cmp, "ms_compound")
#'
#' ## Extract the id, name and inchi
#' cmp_tbl %>% select(compound_id, name, inchi) %>% collect()
NULL

#' @importFrom methods new
#'
#' @exportClass CompDb
.CompDb <- setClass("CompDb",
                    slots = c(dbcon = "DBIConnection",
                              .properties = "list"),
                    prototype = list(.properties = list(),
                                     dbcon = NULL))

#' @importFrom methods validObject
setValidity("CompDb", function(object) {
    if (!is.null(object@dbcon))
        .validCompDb(object@dbcon)
    else TRUE
})

#' @importFrom DBI dbListTables dbGetQuery
.validCompDb <- function(x) {
    txt <- character()
    tables <- dbListTables(x)
    required_tables <- c("ms_compound", "metadata")
    got <- required_tables %in% tables
    if (!all(got))
        txt <- c(txt, paste0("Required tables ", paste0(required_tables[!got]),
                             "not found in the database"))
    ## Check table columns.
    comps <- dbGetQuery(x, "select * from ms_compound limit 3")
    res <- .valid_compound(comps, error = FALSE)
    if (is.character(res))
        txt <- c(txt, res)
    metad <- .metadata(x)
    res <- .valid_metadata(metad, error = FALSE)
    if (is.character(res))
        txt <- c(txt, res)
    if (length(grep("msms", tables))) {
        if (!all(c("msms_spectrum", "msms_spectrum_peak") %in% tables))
            txt <- c(txt, paste0("Required tables msms_spectrum and ",
                                 "msms_spectrum_peak not found in",
                                 " the database"))
        res <- dbGetQuery(x, "select * from msms_spectrum limit 3")
        if (!any(colnames(res) == "spectrum_id"))
            stop("Required column 'spectrum_id' not found in table ",
                 "msms_spectrum")
        res2 <- dbGetQuery(x, "select * from msms_spectrum_peak limit 3")
        if (!all(c("spectrum_id", "mz", "intensity", "peak_id") %in%
                 colnames(res2)))
            stop("Required columns 'spectrum_id', 'mz', 'intensity' and ",
                 "'peak_id' not found in table msms_spectrum_peak")
        res <- dbGetQuery(
            x, paste0("select * from msms_spectrum join msms_spectrum_peak on ",
                      "(msms_spectrum.spectrum_id=msms_spectrum_peak.spectrum_",
                      "id) where compound_id = '", res$compound_id[1], "'"))
        res$predicted <- as.logical(res$predicted)
        res <- .valid_msms_spectrum(res, error = FALSE)
        if (is.character(res))
            txt <- c(txt, res)
        compound_cmp_id <- dbGetQuery(
            x, "select compound_id from ms_compound")[,1]
        spectrum_cmp_id <- dbGetQuery(x, paste0("select compound_id from ",
                                                "msms_spectrum"))[, 1]
        if (!all(spectrum_cmp_id %in% compound_cmp_id))
            txt <- c(txt, paste0("Not all compound ids in the msms_spectrum",
                                 " table are also in the compound table"))
    }
    if (length(txt)) txt else TRUE
}

#' @export
#'
#' @rdname CompDb
CompDb <- function(x, flags = RSQLite::SQLITE_RO) {
    if (missing(x))
        stop("Argument 'x' is required")
    if (is.character(x)) {
        ## Assume it's the file name of the SQLite database
        x <- dbConnect(dbDriver("SQLite"), dbname = x,
                         flags = flags)
    }
    if (is(x, "DBIConnection")) {
        res <- .validCompDb(x)
        if (is.character(res))
            stop(res)
        cdb <- .CompDb(dbcon = x)
        ## fetch all tables and all columns for all tables.
        tbl_nms <- dbListTables(x)
        tbls <- lapply(tbl_nms, function(z) {
            colnames(dbGetQuery(x, paste0("select * from ", z, " limit 1")))
        })
        names(tbls) <- tbl_nms
        cdb@.properties$tables <- tbls
        return(cdb)
    }
    stop("Can not create a 'CompDb' from 'x' of type '", class(x), "'.")
}

#' @importFrom methods is
.metadata <- function(x) {
    if (!is(x, "DBIConnection"))
        x <- .dbconn(x)
    dbGetQuery(x, "select * from metadata")
}

.metadata_value <- function(x, key) {
    metad <- .metadata(x)
    metad[metad$name == key, "value"]
}

.dbconn <- function(x) {
    x@dbcon
}

.has_msms_spectra <- function(x) {
    ## all(c("msms_spectrum_peak", "msms_spectrum_metadata") %in%
    ##     names(.tables(x)))
    any(names(.tables(x)) == "msms_spectrum")
}

#' @description `hasMsMsSpectra` returns `TRUE` if MS/MS spectrum data is
#'     available in the database and `FALSE` otherwise.
#'
#' @export
#'
#' @rdname CompDb
hasMsMsSpectra <- function(x) {
    .has_msms_spectra(x)
}

#' @importFrom dbplyr src_dbi
#'
#' @export
#'
#' @rdname CompDb
src_compdb <- function(x) {
    if (!is(x, "CompDb"))
        stop("'x' is supposed to be a 'CompDb' object")
    src_dbi(.dbconn(x), auto_disconnect = FALSE)
}

#' @description Get a list of all tables and their columns.
#'
#' @param x `CompDb` object.
#'
#' @param name optional `character` to return the table/columns for specified
#'     tables.
#'
#' @param metadata `logical(1)` whether the metadata should be returned too.
#'
#' @noRd
.tables <- function(x, name, metadata = FALSE) {
    tbls <- .get_property(x, "tables")
    if (!missing(name))
        tbls <- tbls[name]
    if (!metadata)
        tbls <- tbls[names(tbls) != "metadata"]
    tbls
}

#' @export
#'
#' @rdname CompDb
tables <- function(x) {
    .tables(x)
}

.get_property <- function(x, name) {
    x@.properties[[name]]
}
