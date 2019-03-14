#' @importFrom utils read.table
#'
#' @noRd
.onLoad <- function(libname, pkgname) {
    path <- system.file("extdata", package = pkgname)
    dta <- read.table(paste0(path, "/adducts.txt"), sep = "\t", as.is = TRUE,
                      header = TRUE)
    add_list <- split(dta, dta$name)
    add_env <- list2env(split(dta, dta$name))
    assign("ADDUCTS", add_env, envir = asNamespace(pkgname))
}
