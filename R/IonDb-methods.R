#' @include IonDb.R

#' @importFrom methods show
#'
#' @export
setMethod("show", "IonDb", function(object) {
    callNextMethod()
    con <- .dbconn(object)
    if (!is.null(con)) {
        ion_nr <- dbGetQuery(con, paste0("select count(distinct ion_id) ",
                                         "from ms_ion"))
        cat(" ion count:", ion_nr[1, 1], "\n")
    }
})

#' @export
#'
#' @rdname IonDb
setMethod("ionVariables", "IonDb", function(object, includeId = FALSE, ...) {
    if (length(.tables(object))) {
        cn <- .tables(object)$ms_ion
        if (includeId)
            cn
        else cn[!cn %in% "ion_id"]
    } else character()
})

#' @importFrom tibble as_tibble
#' 
#' @importMethodsFrom ProtGenerics ions
#'
#' @export
#'
#' @rdname IonDb
setMethod("ions", "IonDb", function(object,
                                    columns = ionVariables(object),
                                    filter,
                                    return.type = c("data.frame",
                                                    "tibble"), ...) {
    return.type <- match.arg(return.type)
    if (length(columns))
        res <- .fetch_data(object, columns = columns, filter = filter,
                           start_from = "ms_ion")
    else res <- data.frame()
    if (return.type == "tibble")
        as_tibble(res)
    else res
})

#' @importFrom DBI dbAppendTable dbGetQuery
#' 
#' @export
#' 
#' @rdname IonDb
setMethod("insertIon", "IonDb", function(object, ions)
{
    .valid_ion(ions, error = TRUE)
    dbcon <- dbConnect(dbDriver("SQLite"), dbname = dbconn(object)@dbname)
    #dbcon <- dbconn(object)
    if (!all(ions$compound_id %in%
             dbGetQuery(dbcon, "select compound_id from ms_compound")[, 1]))
        stop(paste0("All values of 'compound_id' column of 'ions' must be",
                    " in 'compound_id' column of 'ms_compound' table of 'cdb'"))
    ions$ion_id <- seq_len(nrow(ions)) +
        dbGetQuery(dbcon, "select count(distinct ion_id) from ms_ion")[1, 1]
    dbAppendTable(dbcon, "ms_ion", ions)
    dbDisconnect(dbcon)
})
