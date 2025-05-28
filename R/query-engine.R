## Simple query engine to define a query for a e.g. CompDb based on
## provided table and columns

#' @description
#'
#' Utility function to create a SQL query for a `CompDb` database
#' given the provided column names and filter.
#'
#' @param x `CompDb`
#'
#' @details
#'
#' + Check first the parameter `columns` to see if valid column names were
#'   provided.
#' + Based on the columns, get the table name from which the data should be
#'   retrieved (or which tables should be joined if more than one).
#' + `start_from` allows to specify from which table we want to start the join
#'   query. That's important if we don't have all features in all tables
#'   annotated. Example is that we don't have MS/MS spectra from all compounds!
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.build_query_CompDb <- function(x, columns, filter, order, start_from) {
    if (missing(x))
        stop("'x' is required")
    if (missing(columns))
        stop("'columns' is required")
    tbls <- .tables(x)
    col_miss <- columns[!columns %in% unique(unlist(tbls, use.names = FALSE))]
    if (nchar(msg <- paste(col_miss, collapse = ", ")))
        stop("Columns ", msg,
             " are not present in the database. Use 'tables' to list ",
             "all tables and their columns.")
    ## Depending on 'filter' we might have to add some more tables/columns!
    if (!missing(filter)) {
        filter <- .process_filter(filter, x)
        columns_flts <- .field(filter)
        columns <- unique(c(columns, columns_flts))
    }
    ## By default we return also the filter columns!
    columns_tbl <- .reduce_tables_start_from(tbls, columns, start_from)
    paste0(.select(unlist(.prefix_columns(columns_tbl), use.names = FALSE)),
           .from(names(columns_tbl)),
           .where(filter, columns_tbl), .order(order))
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
.select <- function(columns, distinct = TRUE) {
    if (missing(columns))
        stop("No columns provided")
    if (distinct)
        dst <- "distinct "
    else dst <- ""
    paste0("select ", dst, paste0(columns, collapse = ","))
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
.from <- function(tables, joins = .JOINS) {
    q <- .join_tables(tables, joins = joins)
    paste0(" from ", q)
}

.joins <- function(x = new("CompDb")) {
    if (length(x@.properties[["joins"]]))
        x@.properties[["joins"]]
    else .JOINS
}

#' Initialize that in zzz.R and allow adding joins.
#'
#' @noRd
.JOINS <- rbind(
    c("ms_compound", "synonym",
      "on (ms_compound.compound_id=synonym.compound_id)",
      "left outer join"),
    c("ms_ion", "ms_compound",
      "on (ms_ion.compound_id=ms_compound.compound_id)",
      "left outer join"),
    c("ms_compound", "msms_spectrum",
      "on (ms_compound.compound_id=msms_spectrum.compound_id)",
      "left outer join"),
    c("ms_ion", "msms_spectrum",
      "on (ms_ion.compound_id=msms_spectrum.compound_id)",
      "left outer join"),
    c("msms_spectrum", "synonym",
      "on (msms_spectrum.compound_id=synonym.compound_id)",
      "left outer join"),
    c("msms_spectrum", "msms_spectrum_peak",
      "on (msms_spectrum.spectrum_id=msms_spectrum_peak.spectrum_id)",
      "left outer join"),
    c("ms_compound", "experiment",
      "on (ms_compound.expid=experiment.expid)",
      "left outer join")
)

#' Function to convert the "join definition" into a node list
#'
#' @noRd
.table_to_graph <- function(x) {
    u <- unique(as.vector(x))
    res <- lapply(u, function(z) {
        v <- as.vector(x)[x[, 1L] %in% z | x[, 2L] %in% z]
        unique(v[!v %in% z])
    })
    names(res) <- u
    res
}

#' @description Joins two tables from the database.
#'
#' @param x `character` with the names of the tables to be joined.
#'
#' @param joins `character` `matrix` with the definition of the joins that
#'     can be done between database tables.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.join_tables <- function(x, joins = .JOINS){
    x <- unique(x)
    if (length(x) < 2L)
        return(x)
    x <- .add_join_tables(x)
    q <- x[1]
    tbls_used <- x[1]
    tbls <- x[-1]
    idxs <- c(1, 2)
    while(length(tbls)) {
        got_it <- which((joins[, 1] %in% tbls_used & joins[, 2] %in% tbls) |
                        (joins[, 2] %in% tbls_used & joins[, 1] %in% tbls))
        join <- joins[got_it[1], ]
        if (length(got_it)) {
            new_tbl <- join[idxs][!(join[idxs] %in% tbls_used)]
            q <- paste(q, join[4], new_tbl, join[3])
            tbls_used <- c(tbls_used, new_tbl)
            tbls <- tbls[tbls != new_tbl]
        }
        else stop("Tables ", paste(tbls_used, collapse = ", ") ,
                  " can not be joined with ", paste(tbls, collapse = ", "))
    }
    q
}

#' @description
#'
#' Helper function to add additional tables required to join the two provided
#' tables.
#'
#' @param x `character` with the names of the tables to be joined. Has to be
#'     of length > 1.
#'
#' @param join `character` `matrix` with the definition of the tables that
#'     can be joined.
#'
#' @return `character` with all tables needed to join the tables in `x`
#'     (contain `x` plus eventually required additional tables).
#'
#' @md
#'
#' @author Johannes Rainer
#'
#' @importFrom utils combn
#'
#' @noRd
.add_join_tables <- function(x, join = .JOINS) {
    g <- .table_to_graph(join[, 1:2, drop = FALSE])
    cmb <- combn(x, 2, simplify = FALSE)
    res <- unlist(lapply(
        cmb, function(z) .shortest_path(g, start = z[1L], end = z[2L])),
        use.names = FALSE)
    if (!length(res))
        stop("Unable to join tables ", paste0("'", x, "'", collapse = ", "),
             call. = FALSE)
    unique(res)
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
.where <- function(filter, column_list = list()) {
    if (!missing(filter)) {
        if (length(column_list)) {
            nms <- unlist(column_list, use.names = FALSE)
            column_list <- unlist(.prefix_columns(column_list),
                                  use.names = FALSE)
            names(column_list) <- nms
        }
        paste0(" where ", .where_filter(filter, column_list))
    } else NULL
}

#' @description
#'
#' Add a order statement. Thus far we are not testing/checking for correctness.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
.order <- function(order) {
    if (!missing(order))
        paste0(" order by ", order)
    else NULL
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
#' @details
#'
#' This function works well for databases where each entity is annotated/linked
#' to each other. For databases where e.g. not all features in one table a
#' have entries in table b this function might lead to unexpected/unwanted
#' results since table joins are always performed as left (outer joins). In
#' such cases the [.reduce_tables_start_from()] function should be used instead.
#'
#' @param tables `list` of `character`. Names of the list are the table names,
#'     the elements its column names.
#'
#' @param columns `character` with the column names that should be retained.
#'
#' @param start_with `character(1)` with the name of the table from which
#'     the query should start.
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
#' CompoundDb:::.reduce_tables(tabs, columns = c("tx_id", "gene_id"))
#' CompoundDb:::.reduce_tables(tabs, columns = c("gene_id", "gene_name"))
#' CompoundDb:::.reduce_tables(tabs, columns = c("gene_name", "exon_id",
#'    "redundant_field", "tx_id"))
.reduce_tables <- function(tables, columns) {
    ## Get to known in which tables the columns are.
    columns_tbl <- lapply(tables, function(z){
        z[z %in% columns]
    })
    columns_tbl <- columns_tbl[order(lengths(columns_tbl), decreasing = TRUE)]
    ## If we have redundancy, i.e. the same-named column in multiple tables,
    ## remove it from the table with fewer columns.
    tmp_columns <- columns
    for (i in seq_along(columns_tbl)) {
        got_them <- tmp_columns %in% columns_tbl[[i]]
        columns_tbl[[i]] <- tmp_columns[got_them]
        tmp_columns <- tmp_columns[!got_them]
    }
    columns_tbl[lengths(columns_tbl) > 0]
}

#' @description
#'
#' This function is similar to [.reduce_tables()] function but ensures in
#' addition that the table defined by `start_from` is included in the result
#' set and is listed as first table. The function thus can be used for database
#' queries that should contain elements from a certain database table and join
#' other tables to that. See examples below for more details.
#'
#' @note
#'
#' `columns` has to contain at least one column from the database table
#' `start_from`, otherwise `start_from` can not be included resulting in a
#' warning.
#'
#' @param tables see [.reduce_tables()]
#'
#' @param columns see [.reduce_tables()]
#'
#' @param start_from `character(1)` defining the database table name that
#'     should be listed as first table in the result.
#'
#' @md
#'
#' @noRd
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' ## Define database tables with some redundand fields.
#' tabs <- list(
#'     compound = c("compound_id", "name", "red_field"),
#'     spectrum = c("spectrum_id", "compound_id"),
#'     other_tab = c("compound_id", "red_field"))
#'
#' .reduce_tables_start_from(tabs, c("compound_id"))
#' .reduce_tables_start_from(tabs, c("compound_id", "red_field"))
#' .reduce_tables_start_from(tabs, c("compound_id", "red_field"),
#'     start_from = "other_tab")
#' .reduce_tables_start_from(tabs, c("compound_id", "red_field"),
#'     start_from = "spectrum")
#' .reduce_tables_start_from(tabs, c("name", "red_field"),
#'     start_from = "spectrum")
#' .reduce_tables_start_from(tabs, c("name", "red_field"),
#'     start_from = "spectrum_bla")
.reduce_tables_start_from <- function(tables, columns, start_from) {
    tbls <- .reduce_tables(tables, columns)
    if (!missing(start_from)) {
        if (!any(names(tables) == start_from))
            stop("Table '", start_from, "' not known")
        start_from_clms <- intersect(tables[[start_from]], unlist(tbls))
        if (length(start_from_clms)) {
            tbls <- lapply(tbls, function(z) {
                z[!(z %in% start_from_clms)]
            })
            tbls <- tbls[lengths(tbls) > 0]
            tbls <- c(list(start_from_clms), tbls)
            names(tbls)[1] <- start_from
        } else
            warning("No column from table '", start_from, "' in 'columns' ",
                    .call = FALSE)
    }
    tbls
}

.prefix_columns <- function(x) {
    mapply(names(x), x, FUN = function(y, z) paste0(y, ".", z),
                  SIMPLIFY = FALSE)
}

#' Main interface function to retrieve data from the database. Performs the
#' SQL call, gets data, formats data etc.
#'
#' @author Johannes Rainer
#'
#' @noRd
.fetch_data <- function(x, columns, filter, start_from, order) {
    ## If any column is mz or intensity we have to add also spectrum_id, other
    ## wise it's not possible to assign them correctly
    if (any(columns %in% c("mz", "intensity")) & !any(columns == "spectrum_id"))
        columns <- c(columns, "spectrum_id")
    con <- .dbconn(x)
    if (length(.dbname(x)))
        on.exit(dbDisconnect(con))
    res <- dbGetQuery(con, .build_query_CompDb(x, columns = columns,
                                               filter = filter,
                                               start_from = start_from,
                                               order = order))
    if (any(columns == "predicted"))
        res$predicted <- as.logical(res$predicted)
    .deserialize_mz_intensity(res)
}

#' Deserialize m/z and intensity values stored as BLOB in the database.
#'
#' @author Johannes Rainer
#'
#' @noRd
.deserialize_mz_intensity <- function(x) {
    if (nrow(x)) {
        if (is.raw(x$mz[[1]]))
            x$mz <- lapply(x$mz, unserialize)
        if (is.raw(x$intensity[[1]]))
            x$intensity <- lapply(x$intensity, unserialize)
    }
    x
}

#' Counting the number of nodes in a path
#'
#' @noRd
.path_nodes <- function(x = NULL) {
    if (is.null(x)) return(Inf)
    length(x)
}

#' Get the shortest path between two nodes in a graph.
#'
#' @param graph `list` of connected nodes. Names of the `list` are node names
#'     and the elements are the names of the nodes the current node is directly
#'     connected with.
#'
#' @param start `character(1)` with the name of the start node.
#'
#' @param end `character(1)` with the name of the end node.
#'
#' @param path `character` with the path already travelled. Should be empty
#'     for the first call. This will be used in the recursive call to keep
#'     track of the path.
#'
#' @noRd
.shortest_path <- function(graph, start, end, path = c()) {
    if (is.null(graph[[start]])) return(NULL)
    path <- c(path, start)
    ## end of recursion - return path ended in start node
    if (start == end) return(path)
    shortest <- NULL
    for (node in graph[[start]]) {
        if (!node %in% path) {
            newpath <- .shortest_path(graph, node, end, path)
            if (.path_nodes(newpath) < .path_nodes(shortest))
                shortest <- newpath
        }
    }
    shortest
}
