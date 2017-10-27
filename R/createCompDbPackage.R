
#' @description
#'
#' Internal function to extract compound information from a file in SDF format.
#'
#' @param x `character(1)` with the name of the file.
#'
#' @return A [tibble::tibble] with columns `"compound_id"`, `"compound_name"`,
#'     `"inchi"`, `"formula"`, `"mass"`.
#' 
#' @importFrom ChemmineR read.SDFset datablock datablock2ma
#' @importFrom tibble data_frame
#' 
#' @md
#'
#' @author Johannes Rainer
#' 
#' @noRd
.simple_import_compounds_sdf <- function(x) {
    full_mat <- datablock2ma(datablock(read.SDFset(x)))
    source_db <- .guess_sdf_source(colnames(full_mat))
    if (is.null(source_db))
        stop("The SDF file is not supported. Supported are SDF files from ",
             "HMDB, ChEBI and LipidMaps.")
    colmap <- get(paste0(".", source_db, "_colmap"))
    data_frame(compound_id = full_mat[, colmap["id"]],
               compound_name = full_mat[, colmap["name"]],
               inchi = full_mat[, colmap["inchi"]],
               formula = full_mat[, colmap["formula"]],
               mass = as.numeric(full_mat[, colmap["mass"]])
               )
}

#' @description
#'
#' Based on the provided `colnames` guess whether the file is from HMDB,
#' ChEBI or LipidBlast.
#'
#' @param x `character` with the column names of the data table.
#'
#' @return `character(1)` with the name of the resource or `NULL`.
#' 
#' @md
#'
#' @author Johannes Rainer
#' 
#' @noRd
.guess_sdf_source <- function(x) {
    if (any(x == "HMDB_ID"))
        return("hmdb")
    if (any(x == "ChEBI ID"))
        return("chebi")
    if (any(x == "LM_ID"))
        return("lipidmaps")
    NULL
}


.hmdb_colmap <- c(id = "HMDB_ID",
                  name = "GENERIC_NAME",
                  inchi = "INCHI_IDENTIFIER",
                  formula = "FORMULA",
                  mass = "EXACT_MASS"
                  )
.chebi_colmap <- c(id = "ChEBI ID",
                   name = "ChEBI Name",
                   inchi = "InChI",
                   formula = "Formulae",
                   mass = "Monoisotopic Mass"
                   )
.lipidmaps_colmap <- c(id = "LM_ID",
                       name = "COMMON_NAME",
                       inchi = "INCHI",
                       formula = "FORMULA",
                       mass = "EXACT_MASS"
                       )
