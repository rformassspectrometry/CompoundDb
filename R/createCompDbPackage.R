#' @title Extract compound data from a file in SDF format
#'
#' @description
#'
#' `compound_tbl_sdf` extracts basic compound annotations from a file in SDF
#' format (structure-data file). The function currently supports SDF files from:
#' + HMDB (Human Metabolome Database): http://www.hmdb.ca
#' + ChEBI (Chemical Entities of Biological Interest): http://ebi.ac.uk/chebi
#' + LMSD (LIPID MAPS Structure Database).
#'
#' @details
#'
#' Column `"compound_name"` reports for HMDB files the `"GENERIC_NAME"`, for
#' ChEBI the `"ChEBI Name"` and for Lipid Maps the `"COMMON_NAME"`, if that is
#' not available, the first of the compounds synonyms and, if that is also not
#' provided, the `"SYSTEMATIC_NAME"`.
#' 
#' @param file `character(1)` with the name of the SDF file.
#'
#' @param collapse optional `character(1)` to be used to collapse multiple
#'     values in the columns `"synonyms"`. See examples for details.
#' 
#' @return A [tibble::tibble] with general compound information (one row per
#' compound):
#' + `compound_id`: the ID of the compound.
#' + `compound_name`: the compound's name.
#' + `inchi`: the inchi of the compound.
#' + `formula`: the chemical formula.
#' + `mass`: the compound's mass.
#' + `synonyms`: the compound's synonyms (aliases). This type of this column is
#'   by default a `list` to support multiple aliases per compound, unless
#'   argument `collapse` is provided, in which case multiple synonyms are pasted
#'   into a single element separated by the value of `collapse`.
#'
#' @family compound table creation functions
#'
#' @author Johannes Rainer and Jan Stanstrup
#'
#' @export
#' 
#' @md
#' 
#' @examples
#'
#' ## Read compound information from a subset of HMDB
#' fl <- system.file("sdf/HMDB_sub.sdf", package = "CompoundDb")
#' cmps <- compound_tbl_sdf(fl)
#' cmps
#' 
#' ## Column synonyms contains a list
#' cmps$synonyms
#'
#' ## If we provide the optional argument collapse, multiple entries will be
#' ## collapsed.
#' cmps <- compound_tbl_sdf(fl, collapse = "|")
#' cmps
#' cmps$synonyms
compound_tbl_sdf <- function(file, collapse) {
    if (missing(file))
        stop("Please provide the file name using 'file'")
    if (!file.exists(file))
        stop("Can not fine file ", file)
    res <- .simple_import_compounds_sdf(file)
    if (!missing(collapse)) {
        ## collapse elements from lists.
        res$synonyms <- vapply(res$synonyms, paste0, collapse = collapse,
                               FUN.VALUE = "character")
    }
    res
}


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
    sep <- get(paste0(".", source_db, "_separator"))
    ## Fix missing COMMON_NAME entries in LipidMaps (see issue #1)
    nms <- full_mat[, colmap["name"]]
    syns <- strsplit(full_mat[, colmap["synonyms"]], split = sep)
    if (source_db == "lipidmaps") {
        nas <- is.na(nms)
        if (any(nas))
            nms[nas] <- vapply(syns[nas], `[[`, 1, FUN.VALUE = "character")
        nas <- is.na(nms)
        if (any(nas))
            nms[nas] <- full_mat[nas, "SYSTEMATIC_NAME"]
    }
    data_frame(compound_id = full_mat[, colmap["id"]],
               compound_name = nms,
               inchi = full_mat[, colmap["inchi"]],
               formula = full_mat[, colmap["formula"]],
               mass = as.numeric(full_mat[, colmap["mass"]]),
               synonyms = syns
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
                  mass = "EXACT_MASS",
                  synonyms = "SYNONYMS"
                  )
.hmdb_separator <- "; "
.chebi_colmap <- c(id = "ChEBI ID",
                   name = "ChEBI Name",
                   inchi = "InChI",
                   formula = "Formulae",
                   mass = "Monoisotopic Mass",
                   synonyms = "Synonyms"
                   )
.chebi_separator <- " __ "
.lipidmaps_colmap <- c(id = "LM_ID",
                       name = "COMMON_NAME",
                       inchi = "INCHI",
                       formula = "FORMULA",
                       mass = "EXACT_MASS",
                       synonyms = "SYNONYMS"
                       )
.lipidmaps_separator <- "; "

#' @description Import compound information from a LipidBlast file in json
#'     format.
#'
#' @note This is a modified version from Jan's generate_lipidblast_tbl that
#'     extracts the mass also from the json and does not calculate it.
#' 
#' @author Jan Stanstrup and Johannes Rainer
#'
#' @importFrom jsonlite read_json
#' @importFrom dplyr bind_rows
#' 
#' @md
#'
#' @noRd
.import_lipidblast <- function(file) {
    lipidb <- read_json(file)

    parse_element <- function(x) {
        id <- x$id
        cmp <- x$compound[[1]]
        ## get the name(s) -> name + aliases
        nms <- vapply(cmp$names, `[[`, "name", FUN.VALUE = "character")
        mass <- unlist(lapply(cmp$metaData, function(z) {
            if (z$name == "total exact mass")
                z$value
        }))
        if (is.null(mass))
            mass <- NA_character_
        frml <- unlist(lapply(cmp$metaData, function(z) {
            if (z$name == "molecular formula")
                z$value
        }))
        if (is.null(frml))
            mass <- NA_character_
        list(
            id = x$id,
            name = nms[1],
            inchi = cmp$inchi,
            formula = frml,
            mass = mass,
            synonyms = nms[-1]
        )
    }
    
    res <- lapply(lipidb, parse_element)
    bind_rows(res)
}
