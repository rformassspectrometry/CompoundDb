#' @include createCompDbPackage.R

#' @name CompDb
#'
#' @import BiocGenerics
#'
#' @title Simple compound (metabolite) databases
#'
#' @aliases CompDb-class show dbconn,CompDb-method show,CompDb-method
#'
#' @aliases compoundVariables insertSpectra deleteSpectra mass2mz
#'
#' @aliases mass2mz,ANY-method insertCompound deleteCompound
#'
#' @aliases deleteCompound,IonDb-method
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
#' `CompDb()` providing the name of the (SQLite) database file providing
#' the compound annotation data.
#'
#' @section Retrieve annotations from the database:
#'
#' Annotations/compound informations can be retrieved from a `CompDb` database
#' with the `compounds()` and `Spectra()` functions:
#'
#' - `compounds()` extracts compound data from the `CompDb` object. In contrast
#'   to `src_compdb` it returns the actual data as a `data.frame` (if
#'   `return.type = "data.frame"`) or a [tibble::tibble()] (if
#'   `return.type = "tibble"`). A `compounds()` call will always return all
#'   elements from the *ms_compound* table (unless a `filter` is used).
#'
#' - `Spectra()` extract spectra from the database and returns them as a
#'   [Spectra()] object from the *Spectra* package. Additional annotations
#'   requested with the `columns` parameter are added as additional spectra
#'   variables.
#'
#' @section General functions:
#'
#' - `CompDb()`: connect to a compound database.
#'
#' - `compoundVariables()`: returns all available columns/database fields for
#'   compounds.
#'
#' - `copyCompDb()`: allows to copy the content from a CompDb to another
#'   database. Parameter `x` is supposed to be either a `CompDb` or a database
#'   connection from which the data should be copied and `y` a connection to
#'   a database to which it should be copied.
#'
#' - `dbconn()`: returns the connection (of type `DBIConnection`) to the
#'   database.
#'
#' - `metadata()`: returns general meta data of the compound database.
#'
#' - `spectraVariables()`: returns all spectra variables (i.e. columns)
#'   available in the `CompDb`.
#'
#' - `src_compdb()` provides access to the `CompDb`'s database *via*
#'   the functionality from the `dplyr`/`dbplyr` package.
#'
#' - `supportedFilters()`: provides an overview of the filters that can be
#'   applied on a `CompDb` object to extract only specific data from the
#'   database.
#'
#' - `tables()`: returns a named `list` (names being table names) with
#'   the fields/columns from each table in the database.
#'
#' - `mass2mz()`: calculates a table of the m/z values for each compound based
#'   on the provided set of adduct(s). Adduct definitions can be provided with
#'   parameter `adduct`. See [MetaboCoreUtils::mass2mz()] for more details.
#'   Parameter `name` defines the database table column that should be used as
#'   `rownames` of the returned `matrix`. By default `name = "formula"`, m/z
#'   values are calculated for each unique formula in the `CompDb` `x`.
#'
#' @section Adding and removing data from a database:
#'
#' Note that inserting and deleting data requires read-write access to the
#' database. Databases returned by `CompDb` are by default *read-only*. To get
#' write access `CompDb` should be called with parameter
#' `flags = RSQLite::SQLITE_RW`.
#'
#' - `insertCompound()`: adds additional compound(s) to a `CompDb`. The
#'   compound(s) to be added can be specified with parameter `compounds` that
#'   is expected to be a `data.frame` with columns `"compound_id"`, `"name"`,
#'   `"inchi"`, `"inchikey"`, `"formula"`, `"exactmass"`.
#'   Column `"exactmass"` is expected to contain numeric values, all other
#'   columns `character`. Missing values are allowed for all columns except
#'   `"compound_id"`. An optional column `"synonyms"` can be used to provide
#'   alternative names for the compound. This column can contain a single
#'   `character` by row, or a `list` with multiple `character` (names) per
#'   row/compound (see examples below for details). By setting parameter
#'   `addColumns = TRUE` any additional columns in `compound` will be added to
#'   the database table. The default is `addColumns = FALSE`. The function
#'   returns the `CompDb` with the compounds added.
#'   See also [createCompDb()] for more information and details on expected
#'   compound data and the examples below for general usage.
#'
#' - `deleteCompound()`: removes specified compounds from the `CompDb` database.
#'   The IDs of the compounds that should be deleted need to be provided with
#'   parameter `ids`. To include compound IDs in the output of a `compounds()`
#'   call `"compound_id"` should be added to the `columns` parameter. By
#'   default an error is thrown if for some of the specified compounds also MS2
#'   spectra are present in the database. To force deletion of the compounds
#'   along with all associated MS2 spectra use `recursive = TRUE`. See examples
#'   below for details. The function returns the updated `CompDb` database.
#'
#' - `insertSpectra()`: adds further spectra to the database.
#'   The method always adds all the spectra specified through the `spectra`
#'   parameter and does not check if they are already in the database. Note that
#'   the input spectra must have the variable `compound_id` and only `Spectra`
#'   whose `compound_id` values are also in `compounds(object, "compound_id")`
#'   can be added. Parameter `columns` defines which spectra variables from the
#'   `spectra` should be inserted into the database. By default, all spectra
#'   variables are added but it is strongly suggested to specifically select
#'   (meaningful) spectra variables that should be stored in the database.
#'   Note that a spectra variable `"compound_id"` is mandatory.
#'   If needed, the function adds additional columns to the `msms_spectrum`
#'   database table. The function returns the updated `CompDb` object.
#'
#' - `deleteSpectra()`: deletes specified spectra from the database. The IDs of
#'   the spectra to be deleted need to be provided with parameter `ids`.
#'
#' @section Filtering the database:
#'
#' Data access methods such as `compounds()` and `Spectra` allow to filter the
#' results using specific filter classes and expressions. Filtering uses the
#' concepts from Bioconductor's `AnnotationFilter` package. All information
#' for a certain compound with the ID `"HMDB0000001"` can for example be
#' retrieved by passing the filter expression
#' `filter = ~ compound_id == "HMDB0000001"` to the `compounds` function.
#'
#' Use the [supportedFilters()] function on the [CompDb] object to get a list
#' of all supported filters. See also examples below or the usage vignette for
#' details.
#'
#' @param addColumns For `insertCompound()`: `logical(1)` whether all (extra)
#'     columns in parameter `compounds` should be stored also in the database
#'     table. The default is `addColumns = FALSE`.
#'
#' @param columns For `compounds()`, `Spectra`: `character` with the names of
#'     the database columns that should be retrieved. Use `compoundVariables()`
#'     and/or `spectraVariables()` for a list of available column names.
#'     For `insertSpectra()`: columns (spectra variables) that should be
#'     inserted into the database (to avoid inserting all variables).
#'
#' @param compounds For `insertCompound()`: `data.frame` with compound data to
#'     be inserted into a `CompDb` database. See function description for
#'     details.
#'
#' @param filter For `compounds()` and `Spectra()`: filter expression or
#'     [AnnotationFilter()] defining a filter to be used to retrieve specific
#'     elements from the database.
#'
#' @param flags flags passed to the SQLite database connection.
#'     See [SQLite()]. Defaults to read-only, i.e. `RSQLite::SQLITE_RO`.
#'
#' @param ids For `deleteSpectra()`: `integer()`
#'     specifying the IDs of the spectra to delete. IDs in `ids` that are
#'     not associated to any spectra in the `CompDb` object are ignored.
#'     For `deleteCompound`: `character()` with the compound IDs to be deleted.
#'
#' @param includeId for `compoundVariables()`: `logical(1)` whether the comound
#'     ID (column `"compound_id"`) should be included in the result. The
#'     default is `includeIds = FALSE`.
#'
#' @param name For `mass2mz()`: `character(1)`. Defines the `CompDb` column
#'     that will be used to name/identify the returned m/z values. By default
#'     (`name = "formula"`) m/z values for all unique molecular formulas are
#'     calculated and these are used as `rownames` for the returned `matrix`.
#'     With `name = "compound_id"` the adduct m/z for all compounds (even those
#'     with equal formulas) are calculated and returned.
#'
#' @param object For all methods: a `CompDb` object.
#'
#' @param recursive For `deleteCompound()`: `logical(1)` whether also MS2
#'     spectra associated with the compounds should be deleted.
#'
#' @param return.type For `compounds()`: either `"data.frame"` or `"tibble"` to
#'     return the result as a [data.frame()] or [tibble()], respectively.
#'
#' @param x For `CompDb()`: `character(1)` with the file name of the SQLite
#'     compound database. Alternatively it is possible to provide the
#'     connection to the database with parameter `x`. For `copyCompDb()`:
#'     either a `CompDb` or a database connection.
#'
#'     For all other methods: a `CompDb` object.
#'
#' @param y For `copyCompDb()`: connection to a database to which the content
#'     should be copied.
#'
#' @param spectra For `insertSpectra()`: `Spectra` object containing the
#'     spectra to be added to the `IonDb` database.
#'
#' @param ... additional arguments. Currently not used.
#'
#' @return See description of the respective function.
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
#' ## Note that the `compounds` function will by default always return a
#' ## data frame of **unique** entries for the specified columns. Including
#' ## also the `"compound_id"` to the requested columns will ensure that all
#' ## data is returned from the tables.
#' compounds(cdb, columns = c("compound_id", compoundVariables(cdb)))
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
#' ## Calculating m/z values for the exact masses of unique chemical formulas
#' ## in the database:
#' mass2mz(cdb, adduct = c("[M+H]+", "[M+Na]+"))
#'
#' ## By using `name = "compound_id"` the calculation will be performed for
#' ## each unique compound ID instead (resulting in potentially redundant
#' ## results)
#' mass2mz(cdb, adduct = c("[M+H]+", "[M+Na]+"), name = "compound_id")
#'
#' ## Create a Spectra object with all MS/MS spectra from the database.
#' library(Spectra)
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
#'
#' ########
#' ## Creating an empty CompDb and sequentially adding content
#' ##
#' ## Create an empty CompDb and store the database in a temporary file
#' cdb <- emptyCompDb(tempfile())
#' cdb
#'
#' ## Define a data.frame with some compounds to add
#' cmp <- data.frame(
#'     compound_id = c(1, 2),
#'     name = c("Caffeine", "Glucose"),
#'     formula = c("C8H10N4O2", "C6H12O6"),
#'     exactmass = c(194.080375584, 180.063388116))
#'
#' ## We can also add multiple synonyms for each compound
#' cmp$synonyms <- list(c("Cafeina", "Koffein"), "D Glucose")
#' cmp
#'
#' ## These compounds can be added to the empty database with insertCompound
#' cdb <- insertCompound(cdb, compounds = cmp)
#' compounds(cdb)
#'
#' ## insertCompound would also allow to add additional columns/annotations to
#' ## the database. Below we define a new compound adding an additional column
#' ## hmdb_id
#' cmp <- data.frame(
#'     compound_id = 3,
#'     name = "Alpha-Lactose",
#'     formula = "C12H22O11",
#'     exactmass = 342.116211546,
#'     hmdb_id = "HMDB0000186")
#'
#' ## To add additional columns we need to set addColumns = TRUE
#' cdb <- insertCompound(cdb, compounds = cmp, addColumns = TRUE)
#' cdb
#' compounds(cdb)
#'
#' ######
#' ## Deleting selected compounds from a database
#' ##
#' ## Compounds can be deleted with the deleteCompound function providing the
#' ## IDs of the compounds that should be deleted. IDs of compounds in the
#' ## database can be retrieved by adding "compound_id" to the columns parameter
#' ## of the compounds function:
#' compounds(cdb, columns = c("compound_id", "name"))
#'
#' ## Compounds can be deleted with the deleteCompound function. Below we delete
#' ## the compounds with the IDs "1" and "3" from the database
#' cdb <- deleteCompound(cdb, ids = c("1", "3"))
#' compounds(cdb)
#'
#' ## If also MS2 spectra associated with any of these two compounds an error
#' ## would be thrown. Setting the parameter `recursive = TRUE` in the
#' ## `deleteCompound` call would delete the compounds along with their MS2
#' ## spectra.
NULL

setClassUnion("DBIConnectionOrNULL", c("DBIConnection", "NULL"))

#' @importFrom methods new
#'
#' @exportClass CompDb
.CompDb <- setClass("CompDb",
                    slots = c(dbcon = "DBIConnectionOrNULL",
                              .properties = "list",
                              dbname= "character",
                              dbflags = "integer"),
                    prototype = list(.properties = list(),
                                     dbcon = NULL,
                                     dbname = character(),
                                     dbflags = 1L))

#' @importFrom methods validObject
setValidity("CompDb", function(object) {
    con <- .dbconn(object)
    if (!is.null(con)) {
        if (length(.dbname(object)))
            on.exit(dbDisconnect(con))
        .validCompDb(con)
    } else TRUE
})

#' @importFrom DBI dbListTables dbGetQuery dbIsValid
.validCompDb <- function(x) {
    if (!dbIsValid(x))
        return("Database connection not available or closed.")
    tables <- dbListTables(x)
    required_tables <- c("ms_compound", "metadata")
    got <- required_tables %in% tables
    if (!all(got))
        return(paste0("Required tables ",
                      paste0("'", required_tables[!got], "'", collapse = ", "),
                      " not found in the database"))
    ## Check table columns.
    comps <- dbGetQuery(x, "select * from ms_compound limit 3")
    res <- .valid_compound(comps, error = FALSE)
    if (is.character(res))
        return(res)
    metad <- .metadata(x)
    res <- .valid_metadata(metad, error = FALSE)
    if (is.character(res))
        return(res)
    if (length(grep("msms", tables))) {
        if (!all(c("msms_spectrum", "msms_spectrum_peak") %in% tables))
            return(paste0("Required tables msms_spectrum and ",
                          "msms_spectrum_peak not found in the database"))
        res <- dbGetQuery(x, "select * from msms_spectrum limit 3")
        if (!any(colnames(res) == "spectrum_id"))
            return(paste0("Required column 'spectrum_id' not found in table ",
                          "msms_spectrum"))
        res2 <- dbGetQuery(x, "select * from msms_spectrum_peak limit 3")
        if (!all(c("spectrum_id", "mz", "intensity", "peak_id") %in%
                 colnames(res2)))
            return(paste0("Required columns 'spectrum_id', 'mz', 'intensity'",
                          " and 'peak_id' not found in table ",
                          "msms_spectrum_peak"))
        res <- dbGetQuery(
            x, paste0("select * from msms_spectrum join msms_spectrum_peak on ",
                      "(msms_spectrum.spectrum_id=msms_spectrum_peak.spectrum_",
                      "id) where compound_id = '", res$compound_id[1], "'"))
        res$predicted <- as.logical(res$predicted)
        res <- .valid_msms_spectrum(res, error = FALSE)
        if (is.character(res))
            return(res)
        compound_cmp_id <- dbGetQuery(
            x, "select compound_id from ms_compound")[, 1L]
        spectrum_cmp_id <- dbGetQuery(x, paste0("select compound_id from ",
                                                "msms_spectrum"))[, 1L]
        if (!all(spectrum_cmp_id %in% compound_cmp_id))
            return(paste0("Not all compound ids in the msms_spectrum",
                          " table are also in the compound table"))
    }
    TRUE
}

#' @export
#'
#' @importFrom RSQLite SQLITE_RO
#'
#' @rdname CompDb
CompDb <- function(x, flags = SQLITE_RO) {
    if (missing(x))
        stop("Argument 'x' is required and should be either a connection to ",
             "the database or, for SQLite, the database file.")
    if (is.character(x))
        return(.initialize_compdb(.CompDb(dbname = x, dbflags = flags)))
    if (is(x, "DBIConnection"))
        return(.initialize_compdb(.CompDb(dbcon = x, dbflags = flags)))
    stop("'x' should be either a connection to a database or a character ",
         "specifying the (SQLite) database file.")
}

.initialize_compdb <- function(x) {
    con <- .dbconn(x)
    if (length(.dbname(x)) && !is.null(con))
        on.exit(dbDisconnect(con))
    res <- .validCompDb(con)
    if (is.character(res))
        stop(res)
    ## fetch all tables and all columns for all tables.
    tbl_nms <- dbListTables(con)
    tbls <- lapply(tbl_nms, function(z) {
        colnames(dbGetQuery(con, paste0("select * from ", z, " limit 1")))
    })
    names(tbls) <- tbl_nms
    x@.properties$tables <- tbls
    x
}

#' @importFrom methods is
.metadata <- function(x) {
    if (!is(x, "DBIConnection")) {
        n <- .dbname(x)
        x <- .dbconn(x)
        if (length(n) && !is.null(x))
            on.exit(dbDisconnect(x))
    }
    dbGetQuery(x, "select * from metadata")
}

.metadata_value <- function(x, key) {
    metad <- .metadata(x)
    metad[metad$name == key, "value"]
}

#' @importFrom methods .hasSlot
.dbflags <- function(x) {
    if (.hasSlot(x, "dbflags"))
        x@dbflags
    else 1L
}

.dbconn <- function(x) {
    if (length(.dbname(x)))
        dbConnect(dbDriver("SQLite"), dbname = x@dbname, flags = .dbflags(x))
    else x@dbcon
}

.dbname <- function(x) {
    if (.hasSlot(x, "dbname"))
        x@dbname
    else character()
}

.has_msms_spectra <- function(x) {
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


#' @export
#'
#' @rdname CompDb
copyCompDb <- function(x, y) {
    if (inherits(x, "CompDb")) {
        n <- .dbname(x)
        x <- .dbconn(x)
        if (length(n) && !is.null(x))
            on.exit(dbDisconnect(x))
    }
    .copy_compdb(x, y)
}

.require_spectra <- function() {
    requireNamespace("Spectra", quietly = TRUE)
}
