#' @include CompoundDb.R

#' @description `dbconn` returns the connection (`DBIConnection`) to the
#'     database.
#'
#' @importMethodsFrom BiocGenerics dbconn
#'
#' @export
#'
#' @md
#' 
#' @rdname CompoundDb
setMethod("dbconn", "CompoundDb", function(x) {
    .dbconn(x)
})

#' @importFrom methods show
#'
#' @md
#' 
#' @export
setMethod("show", "CompoundDb", function(object) {
    cat("class:", class(object), "\n")
    con <- .dbconn(object)
    if (!is.null(con)) {
        cat(" data source:", .metadata_value(con, "source"), "\n")
        cat(" version:", .metadata_value(con, "source_version"), "\n")
        cat(" organism:", .metadata_value(con, "organism"), "\n")
        cmp_nr <- dbGetQuery(con, "select count(distinct id) from compound")
        cat(" compound count:", cmp_nr[1, 1], "\n")
    } else {
        cat(" no database connection available\n")
    }
})

# organism
