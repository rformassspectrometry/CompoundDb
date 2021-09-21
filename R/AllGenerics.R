setGeneric("annotateMz", function(object, compounds, ...)
    standardGeneric("annotateMz"))

setGeneric("compoundVariables", function(object, ...)
    standardGeneric("compoundVariables"))

setGeneric("ionVariables", function(object, ...)
  standardGeneric("ionVariables"))

setGeneric("insertIon", function(object, ions, ...)
  standardGeneric("insertIon"))

setGeneric("IonDb", function(x, cdb, ...)
    standardGeneric("IonDb"))

#' @importClassesFrom tibble tbl_df
#'
#' @noRd
setClassUnion("DataFrameOrEquivalent", c("data.frame", "DataFrame", "tbl_df"))

setClassUnion("numericOrDataFrameOrEquivalent",
              c("numeric", "DataFrameOrEquivalent"))
