#' @description
#'
#' Given m/z values and names of potential adducts they represent, get all
#' rows of the `cmps` data frame with a matching mass (allowing `ppm`
#' difference).
#'
#' @param x `numeric` with m/z values.
#'
#' @param cmps `data.frame` with a column containing monoisotopic masses in a
#'     column named `"exactmass"`.
#'
#' @param adduct Adduct definition. Either a `data.frame` as returned by
#'     [adducts()] or a `character` with the names of the adducts in that
#'     `data.frame`.
#'
#' @param ppm `numeric(1)` defining the allowed mass difference in ppm.
#'
#' @return a `list` of `data.frame`. Each element containing the results for
#'     each input m/z. Column `"ppm"` in the result `data.frame`s contains the
#'     mass difference in ppm. If nothing is found a `data.frame` with `NA`
#'     values is returned.
#'
#' @noRd
#'
#' @importMethodsFrom S4Vectors extractROWS
#'
#' @importFrom MetaboCoreUtils mz2mass
#'
#' @importFrom MsCoreUtils closest
#'
#' @examples
#'
#' x <- c(76.039304409, 127.036558, 197.078426936, 160.0756)
#' names(x) <- c("Glycine", "3-Hydroxybutyric Acid", "Suberic Acid",
#'     "Serotonin")
#' ## The reported masses are for adducts:
#' ## Glycine: [M+H]+
#' ## 3-Hydroxybutyric Acid: [M+Na]+
#' ## Suberic Acid: [M+Na]+
#' ## Serotonin: [M-NH3+H]+
#'
#' x <- c(105.0546, 209.0851, 97.07362, 100000)
#' names(x) <- c("a", "b", "c", "d")
#' cmps <- compounds(cdb)
#' res <- CompoundDb:::.annotate_adduct_mz(x, cmps, ppm = 5, adduct = c("[M+H]+", "[M+Na]+"))
#'
#' res <- CompoundDb:::.annotate_adduct_mz(x, cmps, ppm = 5, adduct = c("[M+H]+", "[M-NH3+H]+"))
#'
#' ## Negative
#' x <- c(266.089709, 219.0433)
#' names(x) <- c("Adenosine", "Phosphorylcholine")
#'
#' ## Adenosine: [M-H]-
#' ## Phosphorylcholine: [M+Cl]-
#' res <- CompoundDb:::.annotate_adduct_mz(x, cmps, ppm = 5, adduct = c("[M-H]-", "[M+Cl]-"))
.annotate_adduct_mz <- function(x, cmps, adduct = "[M+H]+", ppm = 10,
                                massCol = "exactmass", return.NA = FALSE, ...) {
    x_adds <- split.data.frame(mz2mass(x, adduct = adduct),
                               as.factor(seq_len(length(x))))
    if (!is.null(names(x)))
        names(x_adds) <- names(x)
    cmps <- cmps[order(cmps[, massCol]), ]
    exact_mass <- cmps[, massCol, drop = TRUE]

    lapply(x_adds, function(z) {
        z <- z[, order(z)]                  # get a numeric then.
        hits <- closest(exact_mass, z, tolerance = 0, ppm = ppm, .check = FALSE)
        keep <- !is.na(hits)
        res <- extractROWS(cmps, keep)
        if (any(keep)) {
            keep_adduct <- hits[keep]
            res$.adduct <- names(z)[keep_adduct]
            res$.difference <- abs(exact_mass[keep] - z[keep_adduct])
        } else {
            if (return.NA) {
                res <- cmps[1, ]
                res[,] <- NA
                res$.adduct <- NA_character_
                res$.difference <- NA_real_
            } else {
                res$.adduct <- character()
                res$.difference <- numeric()
            }
        }
        res
    })
}

#' @title Annotate m/z values for ion adducts
#'
#' @aliases annotateMz
#'
#' @rdname annotateMz
#'
#' @description
#'
#' `annotateMz` annotates m/z values using either a `data.frame` containing
#' compound definitions or a [CompDb()] and a user provided list of potential
#' (expected) adducts. In detail, the function calculates for each input m/z
#' value the masses of the potential ion adducts and compares these with
#' (monoisotopic) masses provided in the reference annotation. A hit is returned
#' if the mass difference is smaller than `ppm`.
#'
#' For `object` being a `numeric` and `compounds` a `data.frame` or equivalend
#' class, the function returns a `list` of subsets of `compounds` with rows
#' matching the adducts' mass. If no match is found for the mass of any
#' adducts, a 0-rows `data.frame` (same columns than `compounds` is returned).
#'
#' @param object either a `numeric` of m/z values or a `data.frame` (or
#'     equivalent) with a column containing the m/z values.
#'
#' @param compounds either a `data.frame` (or equivalent) or a [CompDb()] with
#'     the reference annotations. A column named `"exactmass"` is mandatory if
#'     `compounds` is a `data.frame`.
#'
#' @param adduct adduct definition. A `character` defining the adduct(s) to
#'     look for or a `data.frame` with adduct definition. See help on
#'     [mz2mass()] from the `MetaboCoreUtils` package.
#'
#' @param ppm `numeric(1)` defining the allowed mass difference in ppm.
#'
#' @param mzcol for `object` being a `data.frame`: `character(1)` defining the
#'     column name in `object` containing the m/z values. Defaults to `"mz"`.
#'
#' @param ... additional parameters for the `compounds` call if `compounds` is
#'     a `CompDb`.
#'
#' @return See description above for details.
#'
#' @author Johannes Rainer
#'
#' @export
#'
#' @examples
#'
#' ## Read compound annotations from a package-internal SDF file containing
#' ## a small subset of HMDB annotations
#' fl <- system.file("sdf/HMDB_sub.sdf.gz", package = "CompoundDb")
#' cmps <- compound_tbl_sdf(fl)
#' cmps
#'
#' ## Annotate provided m/z values assuming they represent mass-to-charge
#' ## values for [M+H]+ adducts.
#' mzs <- c(105.0546, 75.09168, 127.0365, 113.04756, 113.1210)
#'
#' res <- annotateMz(mzs, cmps, adduct = "[M+H]+")
#'
#' ## res is a `list` with the matching hits for each m/z and the specified
#' ## adduct. If no match was found the result contains only `NA` values.
#' res
#'
#' ## Using all adducts supported by the `mz2mass` function.
#' annotateMz(mzs, cmps, ppm = 10, adduct = MetaboCoreUtils::adductNames())
#'
#' ## It is also possible to annotate an input data.frame, in which case each
#' ## input will be duplicated depending on the number of hits for each m/z
#' mzs_df <- data.frame(id = letters[1:5], mzmed = mzs)
#'
#' res <- annotateMz(mzs_df, cmps, mzcol = "mzmed", adduct = "[M+H]+")
#' res
#'
#' ## Note that `compounds` could also be `CompDb` database, in which a search
#' ## against the annotations provided in that database would be performed.
setMethod(
    "annotateMz", signature(object = "numeric",
                            compounds = "DataFrameOrEquivalent"),
    function(object, compounds, adduct = "[M+H]+", ppm = 10, ...) {
        massCol <- "exactmass"
        if (!any(colnames(compounds) == massCol))
            stop("Required column \"", massCol, "\" not found in 'compounds'",
                 call. = FALSE)
        .annotate_adduct_mz(object, compounds, adduct = adduct, ppm = ppm,
                            massCol = massCol, ...)
    })

#' @rdname annotateMz
setMethod(
    "annotateMz", signature(object = "DataFrameOrEquivalent",
                            compounds = "DataFrameOrEquivalent"),
    function(object, compounds, adduct = "[M+H]+", ppm = 10,
             mzcol = "mz", ...) {
        if (!any(colnames(object) == mzcol))
            stop("Column \"", mzcol, "\" with input m/z values not found ",
                 "in 'object'", call. = FALSE)
        res <- annotateMz(object[, mzcol], compounds = compounds,
                          adduct = adduct, ppm = ppm, return.NA = TRUE)
        nrows <- vapply(res, nrow, integer(1))
        object <- object[rep(1:length(res), nrows), , drop = FALSE]
        rownames(object) <- NULL
        cbind(object, do.call(rbind, res))
    }
)

#' @rdname annotateMz
setMethod(
    "annotateMz", signature(object = c("numericOrDataFrameOrEquivalent"),
                            compounds = "CompDb"),
    function(object, compounds, adduct = "[M+H]+", ppm = 10, ...) {
        annotateMz(object,
                   compounds(object = compounds, ...),
                   adduct = adduct,
                   ppm = ppm, ...)
    }
)
