#' @include createCompDbPackage.R

#' @name CompDb
#'
#' @title Simple compound (metabolite) databases
#'
#' @aliases CompDb-class show dbconn,CompDb-method show,CompDb-method
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
#' with the `compounds` and `spectra` functions:
#'
#' - `compounds` extracts compound data from the `CompDb` object. In contrast
#'   to `src_compdb` it returns the actual data as a `data.frame` (if
#'   `return.type = "data.frame"`) or a [tibble::tibble()] (if
#'   `return.type = "tibble"`). A `compounds` call will always return all
#'   elements from the *compound* table (unless a `filter` is used). Also, the
#'   result `data.frame` will always contain the compound identifier in column
#'   `"compound_id"`.
#'
#' - `spectra` extract spectra from the database and returns them as a
#'   [Spectra()] object. Additional annotations requested with the
#'   `columns` parameter will be added as metadata columns.
#'
#' @section Filtering the database:
#'
#' Data access methods such as `compounds` and `spectra` allow to filter the
#' results using specific filter classes and expressions. Filtering uses the
#' concepts from Bioconductor's `AnnotationFilter` package. All information
#' for a certain compound with the ID `"HMDB0000001"` can for example be
#' retrieved by passing the filter expression
#' `filter = ~ compound_id == "HMDB0000001"` to the `compounds` function.
#'
#' @usage
#' show(object)
#'
#' @param object For all methods: a `CompDb` object.
#'
#' @param x For `CompDb`: `character(1)` with the file name of the SQLite
#'     compound database. Alternatively it is possible to provide the connection
#'     to the database with parameter `x`.
#'
#'     For all other methods: a `CompDb` object.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @seealso
#'
#' [createCompDb()] for the function to create a SQLite compound database.
#'
#' [CompoundIdFilter()] for filters that can be used on the `CompDb` database.
#'
#' @examples
#'
#' ## Create a small CompDb from a provided HMDB subset
#' cmps <- compound_tbl_sdf(system.file("sdf/HMDB_sub.sdf",
#'     package = "CompoundDb"))
#' metad <- data.frame(name = c("source", "url", "source_version",
#'     "source_date", "organism"),
#'     value = c("sub_HMDB", "http://www.hmdb.ca", "4", "2017", "Hsapiens"),
#'     stringsAsFactors = FALSE)
#'
#' ## Load also MS/MS spectra from HMDB xml files
#' xml_path <- system.file("xml", package = "CompoundDb")
#' spctra <- msms_spectra_hmdb(xml_path)
#'
#' ## Create the SQLite database:
#' db_file <- createCompDb(cmps, metadata = metad, msms_spectra = spctra,
#'     path = tempdir())
#'
#' ## Create a CompDb object
#' cmp_db <- CompDb(db_file)
#' cmp_db
#'
#' ## List all tables in the database and their columns
#' tables(cmp_db)
#'
#' ## Extract a data.frame with the id, name and inchi of all compounds
#' compounds(cmp_db, columns = c("compound_id", "compound_name", "inchi"))
#'
#' ## Add also the synonyms (aliases) for the compounds. This will cause the
#' ## tables compound and synonym to be joined. The elements of the compound_id
#' ## and compound_name are now no longer unique
#' res <- compounds(cmp_db, columns = c("compound_id", "compound_name", "synonym"))
#' head(res)
#'
#' ## Extract spectra for a specific HMDB compound.
#' sps <- spectra(cmp_db, filter = ~ compound_id == "HMDB0000001")
#' sps
#'
#' ## Using return.type = "tibble" the result will be returned as a "tibble"
#' compounds(cmp_db, return.type = "tibble")
#'
#' ## Use the CompDb in a dplyr setup
#' library(dplyr)
#' src_cmp <- src_compdb(cmp_db)
#' src_cmp
#'
#' ## Get a tbl for the compound table
#' cmp_tbl <- tbl(src_cmp, "compound")
#'
#' ## Extract the id, name and inchi
#' cmp_tbl %>% select(compound_id, compound_name, inchi) %>% collect()
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
    required_tables <- c("compound", "metadata")
    got <- required_tables %in% tables
    if (!all(got))
        txt <- c(txt, paste0("Required tables ", paste0(required_tables[!got]),
                             "not found in the database"))
    ## Check table columns.
    comps <- dbGetQuery(x, "select * from compound limit 3")
    res <- .valid_compound(comps, error = FALSE)
    if (is.character(res))
        txt <- c(txt, res)
    metad <- .metadata(x)
    res <- .valid_metadata(metad, error = FALSE)
    if (is.character(res))
        txt <- c(txt, res)
    if (length(grep("msms", tables))) {
        ## BLOB
        if (!any(tables == "msms_spectrum"))
            txt <- c(txt, paste0("Required table msms_spectrum not found in",
                                 " the database"))
        compound_cmp_id <- dbGetQuery(x, "select compound_id from compound")[,1]
        spectrum_cmp_id <- dbGetQuery(x, paste0("select compound_id from ",
                                                "msms_spectrum"))[, 1]
        ## Uncomment below if we switch back to individual m/z value storing
        ## req_tables <- c("msms_spectrum_metadata", "msms_spectrum_peak")
        ## got <- req_tables %in% tables
        ## if (!all(got))
        ##     txt <- c(txt, paste0("Required tables ", paste0(req_tables[!got]),
        ##                          "not found in the database"))
        ## compound_cmp_id <- dbGetQuery(x, "select compound_id from compound")[,1]
        ## spectrum_cmp_id <- dbGetQuery(
        ##     x, "select compound_id from msms_spectrum_metadata")[, 1]
        if (!all(spectrum_cmp_id %in% compound_cmp_id))
            txt <- c(txt, paste0("Not all compound ids in the msms_spectrum",
                                 " table are also in the compound table"))
    }
    if (length(txt)) txt else TRUE
}

#' @description `CompDb` *constructs* a `CompDb` object by connecting
#'     to the provided database file.
#'
#' @md
#'
#' @export
CompDb <- function(x) {
    if (missing(x))
        stop("Argument 'x' is required")
    if (is.character(x)) {
        ## Assume it's the file name of the SQLite database, open it read only
        x <- dbConnect(dbDriver("SQLite"), dbname = x,
                         flags = RSQLite::SQLITE_RO)
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

.hasSpectra <- function(x) {
    ## all(c("msms_spectrum_peak", "msms_spectrum_metadata") %in%
    ##     names(.tables(x)))
    any(names(.tables(x)) == "msms_spectrum")
}

#' @description `hasSpectra` returns `TRUE` if MS/MS spectrum data is available
#'     in the database and `FALSE` otherwise.
#'
#' @export
#'
#' @rdname CompDb
#'
#' @md
hasSpectra <- function(x) {
    .hasSpectra(x)
}

#' @param columns For `compounds`, `spectra`: `character` with the names of the
#'     database columns that should be retrieved. Use [tables()] for a list of
#'     available column names.
#'
#' @param filter For `compounds`: not yet supported.
#'
#' @param return.type For `compounds`: `character` defining the type/class of
#'     the return object. Can be either `"data.frame"` (default) or
#'     `"tibble"`.
#'     For `spectra`: either `"Spectra"` (default), `"data.frame"` or
#'     `"tibble"`.
#'
#' @importFrom tibble as_tibble
#'
#' @export
#'
#' @rdname CompDb
#'
#' @md
compounds <- function(x, columns, filter, return.type = "data.frame") {
    if (!is(x, "CompDb"))
        stop("'x' is supposed to be a 'CompDb' object")
    match.arg(return.type, c("data.frame", "tibble"))
    if (missing(columns))
        columns <- .tables(x, "compound")[[1]]
    if (!any(columns == "compound_id"))
        columns <- c("compound_id", columns)
    res <- .fetch_data(x, columns = columns, filter = filter,
                       start_from = "compound")
    if (return.type == "tibble")
        as_tibble(res)
    else res
}

#' @description
#'
#' `src_compdb` provides access to the `CompDb`'s database *via*
#' the functionality from the `dplyr`/`dbplyr` package.
#'
#' @importFrom dbplyr src_dbi
#'
#' @export
#'
#' @rdname CompDb
#'
#' @md
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
#' @md
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

#' @description `tables` returns a named `list` (names being table names) with
#'     the fields/columns from each table in the database.
#'
#' @export
#'
#' @rdname CompDb
#'
#' @md
tables <- function(x) {
    .tables(x)
}

.get_property <- function(x, name) {
    x@.properties[[name]]
}
