## Simple query engine to define a query for a e.g. CompoundDb based on
## provided table and columns

#' @description
#'
#' Utility function to create a SQL query for a `CompoundDb` database
#' given the provided column names and filter.
#'
#' @details
#'
#' + Check first the parameter `columns` to see if valid column names were
#'   provided.
#' + Based on the columns, get the table name from which the data should be
#'   retrieved (or which tables should be joined if more than one).
#' 
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.build_query_CompoundDb <- function(x, columns, filter) {
    if (missing(x))
        stop("'x' is required")
    if (missing(columns))
        stop("'columns' is required")
    if (!missing(filter))
        stop("Not implemented yet")
    tbls <- .tables(x)
    col_ok <- columns %in% unique(unlist(tbls, use.names = FALSE))
    if (!all(col_ok))
        stop("Columns ", paste0(columns[!col_ok], collapse = ", "),
             " are not present in the database. Use 'tables' to list ",
             "all tables and their columns.")
    columns_tbl <- .reduce_tables(tbls, columns)
    if (length(columns_tbl) > 1) {
        stop("Joining of tables not yet implemented")
    }
    paste0(.select(unlist(.prefix_columns(columns_tbl), use.names = FALSE)),
           .from(names(columns_tbl)),
           .where(filter))
}

#' @description
#'
#' Create the *select* part of the SQL query based on the provided columns.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.select <- function(columns) {
    if (missing(columns))
        stop("No columns provided")
    paste0("select ", paste0(columns, collapse = ","))
}

#' @description
#'
#' Create the *from* part of the SQL query based on the provided table names.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.from <- function(tables) {
    q <- ""
    if (length(tables) > 1) {
        stop("Joining of tables is not yet implemented")
    } else
        q <- tables
    paste0(" from ", q)
}

#' @description
#'
#' Create the *where* condition for the SQL query based on the provided filter.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.where <- function(filter) {
    if (!missing(filter)) {
        stop("Filtering is not yet implemented")
    }
    NULL
}

#' @description
#' 
#' Helper function that reduces the provided `list` of table columns to
#' contain only the provided columns. In addition the function ensures that
#' each element in `columns` is only present once across all tables: if a
#' column is present in more than one table, it is retained in the table with
#' most columns and removed from all other tables. See examples below for
#' details.
#'
#' @param tables `list` of `character`. Names of the list are the table names,
#'     the elements its column names.
#'
#' @param columns `character` with the column names that should be retained.
#' 
#' @return `list` being a subset of `tables` that contains only the `columns`.
#' 
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
#'
#' @examples
#'
#' tabs <- list(gene = c("gene_name", "gene_id", "redundant_field"),
#'     tx = c("tx_id", "gene_id", "tx_name", "redundant_field", "tx_start"),
#'     exon = c("exon_id", "redundant_field", "tx_id"))
#'
#' .reduce_tables(tabs, columns = c("tx_id", "gene_id"))
#' .reduce_tables(tabs, columns = c("gene_id", "gene_name"))
#' .reduce_tables(tabs, columns = c("gene_name", "exon_id", "redundant_field",
#'     "tx_id"))
.reduce_tables <- function(tables, columns) {
    ## Get to known in which tables the columns are.
    columns_tbl <- lapply(tables, function(z){
        z[z %in% columns]
    })
    columns_tbl <- columns_tbl[order(lengths(columns_tbl), decreasing = TRUE)]
    ## If we have redundancy, i.e. the same-named column in multiple tables,
    ## remove it from the table with fewer columns.
    tmp_columns <- columns
    for (i in 1:length(columns_tbl)) {
        got_them <- tmp_columns %in% columns_tbl[[i]]
        columns_tbl[[i]] <- tmp_columns[got_them]
        tmp_columns <- tmp_columns[!got_them]
    }
    columns_tbl[lengths(columns_tbl) > 0]
}

.prefix_columns <- function(x) {
    mapply(names(x), x, FUN = function(y, z) paste0(y, ".", z),
                  SIMPLIFY = FALSE)
}


