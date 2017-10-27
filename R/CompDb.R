#' @include createCompoundDbPackage.R query-engine.R

#' @name CompoundDb
#'
#' @title Simple compound (metabolite) databases
#'
#' @aliases CompoundDb-class show,CompoundDb-method dbconn,CompoundDb-method
#'     show dbconn
#' 
#' @description
#'
#' `CompoundDb` objects provide access to general (metabolite) compound
#' annotations along with *metadata* information such as the annotation's
#' source, date and release version. The data is stored internally in a
#' database (usually an SQLite database).
#'
#' @details
#'
#' `CompoundDb` objects should be created using the constructor function
#' `CompoundDb` providing the name of the (SQLite) database file providing
#' the compound annotation data.
#'
#' @usage
#' show(object)
#' 
#' @param object For all methods: a `CompoundDb` object.
#'
#' @param x For `CompoundDb`: `character(1)` with the file name of the SQLite
#'     compound database. Alternatively it is possible to provide the connection
#'     to the database with parameter `x`.
#'
#'     For all other methods: a `CompoundDb` object.
#' 
#' @author Johannes Rainer
#'
#' @md
#'
#' @seealso [createCompoundDb()] for the function to create a SQLite compound
#'     database.
#'
#' @examples
#'
#' ## Create a small CompoundDb from a provided HMDB subset
#' cmps <- generate_hmdb_tbl(system.file("extdata/hmdb/hmdb_sub.xml",
#'     package = "PeakABro"))
#' metad <- data.frame(name = c("source", "url", "source_version",
#'     "source_date", "organism"),
#'     value = c("sub_HMDB", "http://www.hmdb.ca", "4", "2017", "Hsapiens"),
#'     stringsAsFactors = FALSE)
#' ## Create the SQLite database:
#' db_file <- createCompoundDb(cmps, metadata = metad, path = tempdir())
#'
#' ## Create a CompoundDb object
#' cmp_db <- CompoundDb(db_file)
#' cmp_db
#'
#' ## List all tables in the database and their columns
#' tables(cmp_db)
#'
#' ## Extract a data.frame with the id, name and inchi of all compounds
#' compounds(cmp_db, columns = c("id", "name", "inchi"))
#'
#' ## Use the CompoundDb in a dplyr setup
#' library(dplyr)
#' src_cmp <- src_compdb(cmp_db)
#' src_cmp
#'
#' ## Get a tbl for the compound table
#' cmp_tbl <- tbl(src_cmp, "compound")
#'
#' ## Extract the id, name and inchi
#' cmp_tbl %>% select(id, name, inchi) %>% collect()
NULL

#' @exportClass CompoundDb
.CompoundDb <- setClass("CompoundDb",
                        slots = c(dbcon = "DBIConnection",
                                  .properties = "list"),
                        prototype = list(.properties = list(),
                                         dbcon = NULL))

#' @importFrom methods validObject
setValidity("CompoundDb", function(object) {
    if (!is.null(object@dbcon))
        .validCompoundDb(object@dbcon)
    else TRUE
})

#' @importFrom DBI dbListTables dbGetQuery
.validCompoundDb <- function(x) {
    txt <- character()
    tables <- dbListTables(x)
    required_tables <- c("compound", "metadata")
    got <- required_tables %in% tables
    if (!all(got))
        stop("Required tables ", paste0(required_tables[!got]), "not found",
             " in the database")
    ## Check table columns.
    comps <- dbGetQuery(x, "select * from compound limit 3")
    res <- .valid_compound(comps, error = FALSE)
    if (is.character(res))
        txt <- c(txt, res)
    metad <- .metadata(x)
    res <- .valid_metadata(metad, error = FALSE)
    if (is.character(res))
        txt <- c(txt, res)
    if (length(txt)) txt else TRUE
}

#' @description `CompoundDb` *constructs* a `CompoundDb` object by connecting
#'     to the provided database file.
#'
#' @md
#' 
#' @export
CompoundDb <- function(x) {
    if (missing(x))
        stop("Argument 'x' is required")
    if (is.character(x)) {
        ## Assume it's the file name of the SQLite database, open it read only
        x <- dbConnect(dbDriver("SQLite"), dbname = x,
                         flags = RSQLite::SQLITE_RO)
    }
    if (is(x, "DBIConnection")) {
        res <- .validCompoundDb(x)
        if (is.character(res))
            stop(res)
        cdb <- .CompoundDb(dbcon = x)
        ## fetch all tables and all columns for all tables.
        tbl_nms <- dbListTables(x)
        tbls <- lapply(tbl_nms, function(z) {
            colnames(dbGetQuery(x, paste0("select * from ", z, " limit 1")))
        })
        names(tbls) <- tbl_nms
        cdb@.properties$tables <- tbls
        return(cdb)
    }
    stop("Can not create a 'CompoundDb' from 'x' of type '", class(x), "'.")
}

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

#' @description
#'
#' `compounds` extracts compound data from the `CompoundDb` object. In contrast
#' to `src_compdb` it returns the actual data as a `data.frame` (if
#' `return.type = "data.frame"`) or a [tibble::tibble()] (if
#' `return.type = "tibble"`).
#'
#' @param columns For `compounds`: `character` with the names of the database
#'     columns that should be retrieved. Use `tables` for a list of available
#'     column names.
#'
#' @param filter For `compounds`: not yet supported.
#'
#' @param return.type For `compounds`: `character` defining the type/class of
#'     the return object. Can be either `"data.frame"` (default) or
#'     `"tibble"`.
#'
#' @importFrom tibble as_tibble
#'
#' @export
#'
#' @rdname CompoundDb
#' 
#' @md
compounds <- function(x, columns, filter, return.type = "data.frame") {
    if (!is(x, "CompoundDb"))
        stop("'x' is supposed to be a 'CompoundDb' object")
    match.arg(return.type, c("data.frame", "tibble"))
    if (missing(columns))
        columns <- .tables(x, "compound")[[1]]
    res <- dbGetQuery(.dbconn(x),
                      .build_query_CompoundDb(x, columns = columns,
                                              filter = filter))
    if (return.type == "tibble")
        as_tibble(res)
    else res
}

#' @description
#'
#' `src_compdb` provides access to the `CompoundDb`'s database *via*
#' the functionality from the `dplyr`/`dbplyr` package.
#'
#' @importFrom dbplyr src_dbi
#'
#' @export
#' 
#' @rdname CompoundDb
#' 
#' @md
src_compdb <- function(x) {
    if (!is(x, "CompoundDb"))
        stop("'x' is supposed to be a 'CompoundDb' object")
    src_dbi(.dbconn(x), auto_disconnect = FALSE)
}

#' @description Get a list of all tables and their columns.
#'
#' @param x `CompoundDb` object.
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
#' @rdname CompoundDb
#'
#' @md
tables <- function(x) {
    .tables(x)
}

.get_property <- function(x, name) {
    x@.properties[[name]]
}
