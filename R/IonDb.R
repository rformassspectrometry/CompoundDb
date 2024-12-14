#' @include CompDb.R

#' @name IonDb
#'
#' @title IonDb: compound database with additional ion information
#'
#' @aliases IonDb-class show,IonDb-method ionVariables ions insertIon deleteIon
#'
#' @description
#'
#' `IonDb` objects extends `CompDb` by allowing to store also information about
#' measured ions to a [CompDb()] database. This information includes the type
#' (adduct) of the ion, it's measured (or expected) retention time for a certain
#' LC-MS setup and its mass-to-charge ratio.
#'
#' As suggested use case, users might create (or download) a `CompDb` (SQLite)
#' database e.g. containing compound (and eventually MS/MS spectra) annotations
#' from public databases such as the Human Metabolome Database (HMDB) or
#' MassBank. To store now measured ions (e.g. of lab-internal standards) for a
#' certain LC-MS setup, such a `CompDb` can then be converted to an `IonDb`
#' using the `IonDb()` constructor function. Ions can be subsequently added
#' using the `insertIon()` function. In general, it is suggested to create
#' one `IonDb` database for one specific LC-MS setup. Such an `IonDb`
#' database can then be used to match experimental m/z and retention times
#' against ions defined in the database (using the functionality of the
#' [MetaboAnnotation](https://rformassspectrometry.github.io/MetaboAnnotation)
#' package).
#'
#' @section Creation of `IonDb` objects/databases:
#'
#' - A new `IonDb` database can be created and initialized with data from an
#'   existing `CompDb` database by passing either the database connection
#'   (e.g. an `SQLiteConnection`) or the file path of a (to be created) SQLite
#'   database with parameter `x` to the `IonDb()` function and the `CompDb`
#'   object with parameter `cdb`. Optional parameter `ions` allows insert in
#'   addition ion definitions (which can also be added later using
#'   `insertIon()` function calls).
#'
#' - An existing `CompDb` can be converted to an `IonDb` by passing the
#'   [CompDb()] object with parameter `x` to the `IonDb` function. Optional
#'   parameter `ions` allows to provide a `data.frame` with ion definitions to
#'   be inserted in to the database (which can also be added later using
#'   `insertIon()` function calls). Note that this fails if the database
#'   connection for the `CompDb` is read-only.
#'
#' - Previously created `IonDb` databases can be loaded by passing either the
#'   database connection (e.g. an `SQLiteConnection`) or the file path of the
#'   (SQLite) database with parameter `x` to the `IonDb()` function.
#'
#' @section Retrieve annotations and ion information from the database:
#'
#' Annotations/compound informations can be retrieved from a `IonDb` in the
#' same way as thay are extracted from a `CompDb`. In addition, the function
#' `ions()` allows to retrieve the specific ion information from the database.
#' It returns the actual data as a `data.frame` (if
#' `return.type = "data.frame"`) or a [tibble::tibble()]
#' (if `return.type = "tibble"`). An `ions()` call will always
#' return all elements from the *ms_ion* table (unless a `filter` is used).
#'
#' @section General functions (beside those inherited from `CompDb`):
#'
#' - `IonDb()`: connect to or create a compound/ion database.
#'
#' - `ionVariables()`: returns all available columns/database fields for ions.
#'
#' @section Adding and removing data from a database:
#'
#' `IonDb` inherits the `insertCompound()`, `insertSpectra()`,
#' `deleteCompound()` and `deleteSpectra()` functions from [CompDb()].
#' In addition, `IonDb` defines the functions:
#'
#' - `insertIon()`: adds ions to the `IonDb` object. Note that `insertIon()`
#'   always adds all the ions specified through the `ions` parameter
#'   and does not check if they are already in the database. To add columns
#'   present in the submitted `data.frame` to the database table set
#'   `addColumns = TRUE` (default is `addColumns = FALSE`).
#'
#' - `deleteIon()`: deletes ions from the `IonDb` object by specifying
#'    their IDs.
#'
#' @section Filtering the database:
#'
#' Like `compounds()` and `Spectra()` also `ions()` allows to filter the
#' results using specific filter classes and expressions. Filtering uses the
#' concepts from Bioconductor's *AnnotationFilter* package. All information
#' for a certain compound with the ID `"1"` can for example be
#' retrieved by passing the filter expression `filter = ~ ion_id == 1` to
#' the `ions()` function.
#'
#' Use the [AnnotationFilter::supportedFilters()] function on the `IonDb`
#' object to get a list of all supported filters. See also examples below
#' or the usage vignette for details.
#'
#' @param addColumns For `insertIons()`: `logical(1)` whether columns being
#'     present in the submitted `data.frame` but not in the database table
#'     should be added to the database's ion table.
#'
#' @param cdb For `IonDb()`: `CompDb` object from which data should be
#'     transferred to the `IonDb` database.
#'
#' @param columns For `ions()`: `character` with the names of the database
#'     columns that should be retrieved. Use `ionVariables` for a list
#'     of available column names.
#'
#' @param filter For `ions()`: filter expression or
#'     [AnnotationFilter::AnnotationFilter()] defining a filter to be used
#'     to retrieve specific elements from the database.
#'
#' @param flags For `IonDb()`: optional `integer(1)` defining the flags for
#'     the SQLite database connection. Only used if `x` is a `character()`.
#'
#' @param ids For `deleteIon()`: `character()` or (alternatively `integer()`)
#'     specifying the IDs of the ions to delete. IDs in `ids` that are
#'     not associated to any ion in the `IonDb` object are ignored.
#'
#' @param includeId For `ionVariables()`: `logical(1)` whether the ion
#'     ID (column `"ion_id"`) should be included in the result. The
#'     default is `includeId = FALSE`.
#'
#' @param ions for `insertIon()` and `IonDb()`: `data.frame` with ion
#'     definitions to be added to the `IonDb` database. Columns `"compound_id"`
#'     (`character()`), `"ion_adduct"` (`character()`), `"ion_mz"`
#'     (`numeric()`) and `"ion_rt"` (`numeric()`) are mandatory (but, with the
#'     exception of `"compound_id"`, can contain `NA`).
#'
#' @param object For all methods: a `IonDb` object.
#'
#' @param return.type For `ions()`: either `"data.frame"` or `"tibble"` to
#'     return the result as a [data.frame()] or [tibble::tibble()],
#'     respectively.
#'
#' @param x For `IonDb()`: database connection or `character(1)` with the file
#'     name of the SQLite database where the `IonDb` data will be stored or a
#'     [CompDb()] object that should be converted into an `IonDb` object.
#'
#'     For all other methods: an `IonDb` object.
#'
#' @param .DBNAME `character(1)` defining the SQLite database file. This is
#'     an internal parameter not intended to be used/provided by the user.
#'
#' @param ... additional arguments. Currently not used.
#'
#' @return See description of the respective function.
#'
#' @author Andrea Vicini, Johannes Rainer
#'
#' @examples
#'
#' # We load a small compound test database based on MassBank which is
#' # distributed with this package.
#' cdb <- CompDb(system.file("sql/CompDb.MassBank.sql", package = "CompoundDb"))
#' cdb
#'
#' # We next want to convert this CompDb into an IonDb, but the original CompDb
#' # database is read only, thus we have to provide the name (or connection) of
#' # an other database to transfer all the data from the CompDb to that.
#' idb <- IonDb(paste0(tempdir(), "/idb_ex.db"), cdb)
#' idb
#'
#' # It is also possible to load a previously created IonDb passing only the
#' # connection to the database.
#' idb2 <- IonDb(paste0(tempdir(), "/idb_ex.db"))
#'
#' # Ion definitions can be added to the database with the `insertIon` function
#' # providing a `data.frame` with ion definition. This `data.frame` is expected
#' # to provide the IDs of the compounds, an adduct name/definition and the
#' # (experimentally determined) m/z and retention time of the ion. To list
#' # compound IDs from the CompDb database:
#' head(compounds(cdb, c("compound_id", "name")))
#'
#' ions = data.frame(compound_id = c("1", "1", "2", "3", "6", "35"),
#'                   ion_adduct = c("[M+H]+", "[M+Na]+", "[M+Na]+",
#'                                  "[M+Na]+", "[M+2H]2+", "[M+H-NH3]+"),
#'                   ion_mz = c(179.0703, 201.0522, 201.0522,
#'                              201.0522, 253.66982, 312.0390),
#'                   ion_rt = 1:6)
#'
#' # Inserting ion definitions.
#' idb <- insertIon(idb, ions)
#' idb
#'
#' ions(idb, columns = c("name", "formula", "ion_adduct", "ion_mz", "ion_rt"))
#'
#' ## List all available ion variables
#' ionVariables(idb)
#'
#' ## Extract a data.frame with ion variables for all ions
#' ions(idb)
#'
#' ## List all database tables and their columns
#' tables(idb)
#'
#' ## Filtering the database
#' ##
#' ## Get all ions with an m/z between 200 and 300
#' res <- ions(idb, filter = ~ ion_mz > 200 & ion_mz < 300)
#' res
#'
#' ## Get all ions that have a H in their adduct definition.
#' res <- ions(idb, filter = IonAdductFilter("H", "contains"))
#' res
NULL

#' @importFrom methods new
#'
#' @exportClass IonDb
.IonDb <- setClass("IonDb",
                   contains = "CompDb")

#' @importFrom methods validObject
setValidity("IonDb", function(object) {
    con <- .dbconn(object)
    if (length(.dbname(object)) && !is.null(con))
        on.exit(dbDisconnect(con))
    if (!is.null(con))
        .validIonDb(con)
    else TRUE
})

#' @importFrom DBI dbListTables dbGetQuery
.validIonDb <- function(x) {
    txt <- character()
    if (!dbIsValid(x))
        return("Database connection not available or closed.")
    if(!("ms_ion" %in% dbListTables(x)))
        txt <- c(txt, "Required table 'ms_ion' not found not found in database")
    else {
       ions <- dbGetQuery(x, "select * from ms_ion limit 3")
       res <- .valid_ion(ions, error = FALSE)
       if (is.character(res))
           txt <- c(txt, res)
    }
    if (length(txt)) txt else TRUE
}

#' @description Check that the ions table contains all required data.
#'
#' @param db `logical(1)` whether validity should be checked on the internal
#'     database table instead of the input file.
#' @md
#'
#' @noRd
.valid_ion <- function (ions, error = TRUE){
    txt <- .valid_data_frame_columns(ions, "ions", .required_ion_columns)
    if (!length(txt)) {
        if (!is.character(ions$compound_id))
          txt <- c(txt, "Column 'compound_id' should be of type character")
        if (!is.character(ions$ion_adduct))
            txt <- c(txt, "Column 'ion_adduct' should be of type character")
        if (!is.numeric(ions$ion_mz))
            txt <- c(txt, "Column 'ion_mz' should be numeric")
        if (!is.numeric(ions$ion_rt))
            txt <- c(txt, "Column 'ion_rt' should be numeric")
        if (nrow(ions) && (any(is.na(ions$compound_id)) ||
                           any(ions$compound_id == "")))
            txt <- c(txt, "No missing values in column 'compound_id' allowed")
    }
    .throw_error(txt, error = error)
}

#' @importMethodsFrom DBI dbDataType
#'
#' @noRd
.create_ion_table <- function(con, ions = data.frame(ion_id = integer(),
                                                     compound_id = character(),
                                                     ion_adduct = character(),
                                                     ion_mz = numeric(),
                                                     ion_rt = numeric())) {
    dbtype <- dbDataType(con, ions)
    sql <- paste0("CREATE TABLE ms_ion (",
                  paste(names(dbtype), dbtype, collapse = ","), ");")
    suppressWarnings(res <- dbExecute(con, sql))
    res <- dbExecute(
        con, "create index ms_ion_compound_id_idx on ms_ion (compound_id)")
}

#' Copy all tables from a compdb to another database creating also all indices.
#'
#' @noRd
.copy_compdb <- function(x, y) {
    tbls <- dbListTables(x)
    lapply(tbls, function(tbl)
        dbWriteTable(y, tbl, dbGetQuery(x, paste0("select * from ", tbl))))
    dbExecute(y, "create index compound_id_idx on ms_compound (compound_id)")
    dbExecute(y, "create index compound_name_idx on ms_compound (name)")
    if (any(tbls == "msms_spectrum")) {
        dbExecute(y, "create index msms_mid_idx on msms_spectrum (spectrum_id)")
        dbExecute(y, "create index msms_cid_idx on msms_spectrum (compound_id)")
    }
    if (any(tbls == "msms_spectrum_peak")) {
        dbExecute(
            y, "create index msms_id_idx on msms_spectrum_peak (spectrum_id)")
        dbExecute(
            y, "create index msms_pid_idx on msms_spectrum_peak (peak_id)")
    }
}
