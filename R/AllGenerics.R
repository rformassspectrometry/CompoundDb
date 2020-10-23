setGeneric("annotateMz", function(object, compounds, ...)
    standardGeneric("annotateMz"))

setGeneric("compoundVariables", function(object, ...)
    standardGeneric("compoundVariables"))

#' @importClassesFrom tibble tbl_df
#'
#' @noRd
setClassUnion("DataFrameOrEquivalent", c("data.frame", "DataFrame", "tbl_df"))

setClassUnion("numericOrDataFrameOrEquivalent",
              c("numeric", "DataFrameOrEquivalent"))
