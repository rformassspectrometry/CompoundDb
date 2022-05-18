setGeneric("annotateMz", function(object, compounds, ...)
    standardGeneric("annotateMz"))

setGeneric("compoundVariables", function(object, ...)
    standardGeneric("compoundVariables"))

setGeneric("ionVariables", function(object, ...)
  standardGeneric("ionVariables"))

setGeneric("insertCompound", function(object, compounds, ...)
  standardGeneric("insertCompound"))

setGeneric("insertIon", function(object, ions, ...)
  standardGeneric("insertIon"))

setGeneric("insertSpectra", function(object, spectra, ...)
  standardGeneric("insertSpectra"))

setGeneric("deleteCompound", function(object, ...)
  standardGeneric("deleteCompound"))

setGeneric("deleteIon", function(object, ...)
  standardGeneric("deleteIon"))

setGeneric("deleteSpectra", function(object, ...)
  standardGeneric("deleteSpectra"))

setGeneric("IonDb", function(x, cdb, ...)
    standardGeneric("IonDb"))

setGeneric("mass2mz", function(x, ...)
  standardGeneric("mass2mz"))

#' @importClassesFrom tibble tbl_df
#'
#' @noRd
setClassUnion("DataFrameOrEquivalent", c("data.frame", "DataFrame", "tbl_df"))

setClassUnion("numericOrDataFrameOrEquivalent",
              c("numeric", "DataFrameOrEquivalent"))
