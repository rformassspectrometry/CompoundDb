setGeneric("hasMz", function(object, mz, ...) standardGeneric("hasMz"))

setGeneric("annotateMz", function(object, compounds, ...)
    standardGeneric("annotateMz"))

#' @importClassesFrom tibble tbl_df
#'
#' @noRd
setClassUnion("DataFrameOrEquivalent", c("data.frame", "DataFrame", "tbl_df"))

setClassUnion("numericOrDataFrameOrEquivalent",
              c("numeric", "DataFrameOrEquivalent"))
