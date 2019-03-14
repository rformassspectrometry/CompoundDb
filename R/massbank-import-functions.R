
#' Import a single massbank file.
#'
#' @param x `character(1)` with the file name of the Massbank txt file. The
#'     file will be read with `readLines`.
#'
#' @return A single-row `tibble` or a 0-row tibble if none of the fields were
#'     found in the file.
#'
#' @noRd
#'
#' @importFrom tibble as_tibble
#'
#' @author Johannes Rainer
.import_massbank_file <- function(x, collapsed = TRUE) {
    lns <- readLines(x)
    ## Import the data, return as a data.frame
    res <- mapply(.massbank_fields$field, .massbank_fields$data_type,
                  FUN = function(field, type) {
                      vls <- .massbank_extract_field(field = field, lns)
                      if (type == "list")
                          list(vls)
                      else as(vls, type)
    }, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    names(res) <- .massbank_fields$column
    res <- as_tibble(res)
    if (all(is.na(res)))
        res[integer(), ]
    else res
}

#' Extract a field from a Massbank txt file.
#'
#' @param field `character(1)` with the name of the field
#'
#' @param x `character` e.g. representing the content of an txt file read with
#'     `readLines`.
#'
#' @return `character` with the value for the field. Can be of length > 1 if
#'     multiple rows match the field. Returns `NA_character_` if the field
#'     is not present.
#'
#' @author Johannes Rainer
#'
#' @noRd
.massbank_extract_field <- function(field = "CH$NAME: ", x) {
    idx <- grep(field, x, fixed = TRUE)
    if (length(idx))
        sub(field, "", x[idx], fixed = TRUE)
    else NA_character_
}

.massbank_fields <- data.frame(
    field = c("CH$NAME: ",
              "CH$FORMULA: ",
              "CH$EXACT_MASS: ",
              "CH$IUPAC: ",
              "CH$LINK: INCHIKEY ",
              "CH$LINK: CHEBI ",
              "CH$LINK: KEGG ",
              "CH$LINK: PUBCHEM ",
              "CH$LINK: CHEMSPIDER "),
    column = c("synonyms",
               "formula",
               "mass",
               "inchi",
               "inchikey",
               "chebi_id",
               "kegg_id",
               "pubchem_id",
               "chemspider_id"),
    data_type = c("list",
                  "character",
                  "numeric",
                  "character",
                  "character",
                  "character",
                  "character",
                  "character",
                  "character"),
    stringsAsFactors = FALSE
)
