## Utility functions related to mass/adduct calculations and similar.

.ppm <- function(x, ppm = 10) {
    ppm * x / 1e6
}

#' @title Match values allowing a certain difference in ppm
#'
#' @description
#'
#' Match numeric values in `x` against values in `y` and return indices of
#' **all** matches of `x` in `y` allowing differences between the values
#' specified with parameter `ppm` (part per million).
#'
#' @param x `numeric` with input values to match in `y`. Supposed to be a
#'     shorter `numeric` than `y`.
#'
#' @param y `numeric` with values to match `x` against.
#'
#' @param ppm `numeric(1)` defining the allowed difference between values to be
#'     still considered a match. Differences smaller than +/- `ppm` are
#'     accepted.
#'
#' @return `list` with `integer` representing the index in `y` where each
#'     element in `x` matches (given the provided ppm).
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @export
#'
#' @examples
#'
#' yvals <- abs(rnorm(1000))
#' yvals[123] <- yvals[2] - yvals[2] * 10 * 1e-6
#' yvals[124] <- yvals[2] + yvals[2] * 10 * 1e-6
#' yvals[125] <- yvals[2] + yvals[2] * 12 * 1e-6
#' xvals <- yvals[c(2, 3, 3, 20, 21, 20)]
#' xvals[2] <- xvals[2] + (10 * xvals[2] / 1e6)
#' xvals[3] <- xvals[3] - (10 * xvals[3] / 1e6)
#' xvals[6] <- xvals[6] + (12 * xvals[6] / 1e6)
#'
#' ## Perfect matches:
#' matchWithPpm(xvals, yvals)
#'
#' ## Match allowing +/- 10ppm difference
#' matchWithPpm(xvals, yvals, ppm = 10)
#'
#' ## Match allowing +/- 20ppm difference
#' matchWithPpm(xvals, yvals, ppm = 20)
matchWithPpm <- function(x, y, ppm = 0) {
    lapply(x, function(z, ppm) {
        which(abs(z - y) <= (.ppm(z, ppm)) + 1e-9)
    }, ppm = force(ppm))
}

#' @title Conversion functions between mass and m/z
#'
#' @description
#'
#' `mass2mz` and `mz2mass` allow to convert between (monoisotopic) mass and m/z
#' of specified ion adducts an *vice versa*. The functions return a `list` with
#' length equal to `x`, each element being a `numeric` with the m/z or mass for
#' the specified adducts. See below for examples.
#'
#' `adducts` retrieves the definitions for the supported ion adducts. The
#' function returns a `data.frame`.
#'
#' @param x `numeric` with the masses (or m/z values) to be converted.
#'
#' @param adduct either a `data.frame` with required columns `"name"`, `"nmol"`
#'     (number of molecules), `"charge"` (the total charge of the molecule)
#'     and `"massdiff"` (total mass difference), such as returned by the
#'     `adducts` function, or a `character` with the names of the adducts
#'     (see `adducts` for supported adducts).
#'
#' @rdname mass2mz
#'
#' @author Johannes Rainer, Jan Stanstrup
#'
#' @export
#'
#' @examples
#'
#' masses <- c(75.032028409, 105.042595, 162.115698458, 180.063385)
#' names(masses) <- c("a", "b", "c", "d")
#'
#' ## Calculate mz for adducts [M+H]+ and [M+Na]+
#' mzs <- mass2mz(masses, adduct = c("[M+H]+", "[M+Na]+"))
mass2mz <- function(x, adduct = adducts()) {
    if (is.character(adduct))
        adduct <- adducts(name = adduct)
    adduct$charge <- abs(adduct$charge)
    lapply(x, function(z) {
        res <- (adduct$nmol * z + adduct$massdiff) / adduct$charge
        names(res) <- adduct$name
        res
    })
}

#' @export
#'
#' @rdname mass2mz
mz2mass <- function(x, adduct = adducts()) {
    if (is.character(adduct))
        adduct <- adducts(name = adduct)
    adduct$charge <- abs(adduct$charge)
    lapply(x, function(z) {
        res <- (adduct$charge * z - adduct$massdiff) / adduct$nmol
        names(res) <- adduct$name
        res
    })
}

#' @param pattern For `adducts`: optional `character(1)` specifying a pattern to
#'     be used to retrieve selected adducts, e.g. containing a hydrogen.
#'
#' @param polarity For `adducts`: optional `numeric(1)` to retrieve only adducts
#'     with positive (`polarity = 1`) or negative (`polarity = -1`) polarity.
#'
#' @param name For `adducts`: define the names of the adducts for which the
#'     definition should be returned.
#'
#' @param set For `adducts`: `character(1)` defining sets of adducts.
#'
#' @param ... For `adducts`: additional parameters to be passed to the [grep()]
#'     function.
#'
#' @export
#'
#' @rdname mass2mz
adducts <- function(pattern, polarity, name, set, ...) {
    adds <- ADDUCTS
    if (!missing(polarity)) {
        if (polarity < 0)
            adds <- adds[adds$charge < 0, ]
        else adds <- adds[adds$charge > 0, ]
    }
    if (!missing(pattern))
        adds <- adds[grep(pattern, adds$name, ...), ]
    if (!missing(name))
        adds <- adds[match(name, adds$name), ]
    adds
}

#' @description
#'
#' Given m/z values and names of potential adducts they represent, get all
#' rows of the `cmps` data frame with a matching mass (allowing `ppm`
#' difference).
#'
#' @param x `numeric` with m/z values.
#'
#' @param cmps `data.frame` with a column containing monoisotopic masses in a
#'     column named `"mass"`.
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
.annotate_adduct_mz <- function(x, cmps, adduct = adducts(), ppm = 10) {
    if (is.character(adduct))
        adduct <- adducts(name = adduct)
    adduct$charge <- abs(adduct$charge)
    not_found <- cmps[1, , drop = FALSE]
    rownames(not_found) <- NULL
    not_found[1, ] <- lapply(sapply(not_found, class), as, object = NA)
    not_found$adduct <- NA_character_
    not_found$ppm <- NA_real_
    lapply(x, function(z) {
        mss <- (adduct$charge * z - adduct$massdiff) / adduct$nmol
        idx <- matchWithPpm(mss, cmps$mass, ppm = ppm)
        hits <- lengths(idx)
        idx <- unlist(idx, use.names = FALSE)
        if (length(idx)) {
            res <- cmps[idx, , drop = FALSE]
            rownames(res) <- NULL
            res$adduct <- rep(adduct$name, hits)
            rep_mss <- rep(mss, hits)
            res$ppm <- abs(res$mass - rep_mss) * 1e6 / res$mass
            res[order(res$ppm), ]
        } else not_found
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
#' adducts, a row of `compounds` with all values set to `NA` is returned. Note
#' that if `compounds` is ` [CompDb()], a `list` of `data.frame`s is returned.s
#'
#' @param object either a `numeric` of m/z values or a `data.frame` (or
#'     equivalent) with a column containing the m/z values.
#'
#' @param compounds either a `data.frame` (or equivalent) or a [CompDb()] with
#'     the reference annotations. A column named `"mass"` is mandatory if
#'     `compounds` is a `data.frame`.
#'
#' @param adduct adduct definition. Either a `data.frame` as returned by
#'     [adducts()] or a `character` with the names of the adduct(s) in that
#'     `data.frame`.
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
#' ## Annotating using all adducts defined by `adducts()` and accepting a
#' ## 20ppm difference between the masses
#' annotateMz(mzs, cmps, ppm = 20)
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
    function(object, compounds, adduct = adducts(), ppm = 10, ...) {
        if (!any(colnames(compounds) == "mass"))
            stop("Required column \"mass\" not found in 'compounds'",
                 call. = FALSE)
        .annotate_adduct_mz(object, compounds, adduct = adduct, ppm = ppm)
    })

#' @rdname annotateMz
setMethod(
    "annotateMz", signature(object = "DataFrameOrEquivalent",
                            compounds = "DataFrameOrEquivalent"),
    function(object, compounds, adduct = adducts(), ppm = 10,
             mzcol = "mz", ...) {
        if (!any(colnames(object) == mzcol))
            stop("Column \"", mzcol, "\" with input m/z values not found ",
                 "in 'object'", call. = FALSE)
        res <- annotateMz(object[, mzcol], compounds = compounds,
                          adduct = adduct, ppm = ppm)
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
    function(object, compounds, adduct = adducts(), ppm = 10, ...) {
        annotateMz(object, compounds(x = compounds, ...), adduct = adduct,
                   ppm = ppm, ...)
    }
)
