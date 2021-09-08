#' @include createCompDbPackage.R

#' @name IonDb
#'
#' @title Simple database containing compound (metabolite) and ion information 
#'
#' @aliases IonDb-method show,IonDb-method ionVariables
#'
#' @description
#' 
#' `IonDb` objects extends `CompDb` by adding ions measured on a certain LC-MS 
#' setup. The data is stored internally in a database (usually an SQLite 
#' database).
#'
#' @details
#'
#' `IonDb` objects should be created using the constructor function
#' `IonDb` providing a database connection, a `CompDb` object and a table with 
#' ions (the content of the `CompDb` object and ion table is copied to
#' the database pointed by the connection). Alternatively it is possible to 
#' directly pass just the database connection to an already existing database 
#' (with all the required tables) and load it as a `IonDb`.
#'
#' @section Retrieve annotations and ion information from the database:
#'
#' Annotations/compound informations can be retrieved from a `IonDb` in the same
#' way as thay are extracted from a `CompDb`. In addition to that the function 
#' `ions` allows to retrieve ion information from the `IonDb` object. It returns 
#' the actual data as a `data.frame` (if `return.type = "data.frame"`) or a 
#' [tibble::tibble()] (if `return.type = "tibble"`). A `ions` call will always 
#' return all elements from the *ms_ion* table (unless a `filter` is used).
#'
#' @section General functions (beside those inherited from `CompDb`):
#'
#' - `IonDb`: connect to a compound/ion database.
#'
#' - `ionVariables`: returns all available columns/database fields for ions.
#' 
#' - `insertIon`: allows to add further ions to the `IonDb` object.
#'
#'
#' @section Filtering the database:
#'
#' Like `compounds` and `Spectra` also `ions` allows to filter the
#' results using specific filter classes and expressions. Filtering uses the
#' concepts from Bioconductor's `AnnotationFilter` package. All information
#' for a certain compound with the ID `"1"` can for example be
#' retrieved by passing the filter expression `filter = ~ ion_id == "1"` to 
#' the `ions` function.
#'
#' Use the [supportedFilters] function on the [IonDb] object to get a list of
#' all supported filters. See also examples below or the usage vignette for
#' details.
#'
#' @param columns For `ions`: `character` with the names of the database 
#'     columns that should be retrieved. Use `ionVariables` for a list 
#'     of available column names.
#'
#' @param filter For `ions`: filter expression or [AnnotationFilter()] defining 
#'     a filter to be used to retrieve specific elements from the database.
#'
#' @param includeId for `ionVariables`: `logical(1)` whether the ion
#'     ID (column `"ms_ion_id"`) should be included in the result. The
#'     default is `includeId = FALSE`.
#'
#' @param object For all methods: a `IonDb` object.
#'
#' @param return.type For `ions`: either `"data.frame"` or `"tibble"` to
#'     return the result as a [data.frame()] or [tibble()], respectively.
#'
#' @param x For `IonDb`: database connection (should be to empty database?) or 
#'     `character(1)` with the file name of the SQLite database where the 
#'     `IonDb` data will be stored.
#'
#'     For all other methods: a `IonDb` object.
#' 
#' @param cdb for `IonDb`: `CompDb` object used to construct a `IonDb` object 
#'     (the content of `cdb` is copied to the `IonDb` object).
#' 
#' @param ions for `IonDb`: table to be added to the `IonDb` object (if missing
#'     an empty table is added). `ions` is required to have the columns 
#'     "compound_id", "ion_adduct", "ion_mz" and "ion_rt". 
#'
#' @param ... additional arguments. Currently not used.
#'
#' @author Andrea Vicini, Johannes Rainer
#'
#' @examples
#'
# # To be added.
NULL

#' @importFrom methods new
#'
#' @exportClass IonDb
.IonDb <- setClass("IonDb",
                   contains = "CompDb")

#' @importFrom methods validObject
setValidity("IonDb", function(object) {
    if (!is.null(object@dbcon))
        .validIonDb(object@dbcon)
    else TRUE
})

#' @importFrom DBI dbListTables dbGetQuery
.validIonDb <- function(x) {
    txt <- character()
    # maybe this is already called automatically by the validation of base class?
    #txt <- c(txt, .validCompDb(x))            
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

#' @export
#'
#' @rdname IonDb
IonDb <- function(x, cdb, ions) {
    if (missing(x))        
        stop("Argument 'x' is required")
    if (is.character(x)) {
        ## Assume it's the file name of the SQLite database
        x <- dbConnect(dbDriver("SQLite"), dbname = x)
    }
    if (!is(x, "DBIConnection"))
        stop("Can not create a 'IonDb' from 'x' of type '",
             class(x), "'.")
    if (!missing(cdb)) {
        # Not sure about the following condition. Should we also handle in some 
        # ways connections to databases with tables in it?
        if (length(dbListTables(x)))
            stop(paste0("'x' must be a connection to a empty database when", 
                        " 'cdb' is provided"))
        con_cdb <- dbconn(cdb)
        sapply(dbListTables(con_cdb), function(tbl)
            dbWriteTable(x, tbl, dbGetQuery(con_cdb,
                                                paste0("select * from ", tbl))))
        if(missing(ions)) {
            ions <- data.frame(compound_id = character(0), adduct = character(0),
                               mz = numeric(0), rt = numeric(0),
                               ion_id = character(0))
            dbWriteTable(x, "ms_ion", ions)
        } else {
            txt <- character(0)
            res <- .valid_ion(ions)
            if(is.character(res))
                txt <- c(txt, res)
            if (!all(ions$compound_id %in%
                     dbGetQuery(con_cdb,
                                "select compound_id from ms_compound")[, 1]))
                txt <- c(txt, paste0("All values of 'compound_id' column of ", 
                                     "'ions' must be in 'compound_id' column ",
                                     "of 'ms_compound' table of 'cdb'"))
            if (length(txt))
                stop(paste(txt, collapse = "\n"))
            #ions$ion_id <- seq_len(nrow(ions))
            ions <- data.frame(ion_id = seq_len(nrow(ions)), ions)
            dbWriteTable(x, "ms_ion", ions, overwrite = TRUE)
        }
    }
    res <- .validCompDb(x)
    if (is.character(res))
        stop(res)
    res <- .validIonDb(x)
    if (is.character(res))
        stop(res)
    
    idb <- .IonDb(dbcon = x)
    tbl_nms <- dbListTables(x)
    tbls <- lapply(tbl_nms, function(z) {
        colnames(dbGetQuery(x, paste0("select * from ", z, " limit 1")))
    })
    names(tbls) <- tbl_nms
    idb@.properties$tables <- tbls
    return(idb)
}

