#' @include IonDb.R

#' @importFrom methods show
#'
#' @export
setMethod("show", "IonDb", function(object) {
    cat("class:", class(object), "\n")
    con <- .dbconn(object)
    if (!is.null(con)) {
        cat(" data source:", .metadata_value(con, "source"), "\n")
        cat(" version:", .metadata_value(con, "source_version"), "\n")
        cat(" organism:", .metadata_value(con, "organism"), "\n")
        cmp_nr <- dbGetQuery(con, paste0("select count(distinct compound_id) ",
                                         "from ms_compound"))
        cat(" compound count:", cmp_nr[1, 1], "\n")
        ion_nr <- dbGetQuery(con, paste0("select count(distinct ion_id) ",
                                         "from ms_ion"))
        cat(" ion count:", ion_nr[1, 1], "\n")
        if (.has_msms_spectra(object)) {
            spctra <- dbGetQuery(con, paste0("select count(distinct spectrum_",
                                             "id) from msms_spectrum"))
            cat(" MS/MS spectra count:", spctra[1, 1], "\n")
        }
    } else cat(" no database connection available\n")
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
    # should I check here that the newly inserted ions are not already in the 
    # object, right? or we can add two different ions with the same id?
    ions$ion_id <- seq_len(nrow(ions)) +
        dbGetQuery(dbcon, "select count(distinct ion_id) from ms_ion")[1, 1]
    dbAppendTable(dbcon, "ms_ion", ions)
    dbDisconnect(dbcon)
})
