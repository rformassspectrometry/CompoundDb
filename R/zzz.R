#' @importFrom commonMZ MZ_CAMERA
#'
#' @noRd
.onLoad <- function(libname, pkgname) {
    dta <- as.data.frame(rbind(commonMZ::MZ_CAMERA("pos", warn_clash = FALSE),
                               commonMZ::MZ_CAMERA("neg", warn_clash = FALSE)))
    assign("ADDUCTS", dta, envir = asNamespace(pkgname))
}
