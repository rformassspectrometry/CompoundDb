test_that("matchWithPpm works", {
    set.seed(123)
    yvals <- abs(rnorm(1000))
    yvals[123] <- yvals[2] - yvals[2] * 10 * 1e-6
    yvals[124] <- yvals[2] + yvals[2] * 10 * 1e-6
    yvals[125] <- yvals[2] + yvals[2] * 12 * 1e-6
    xvals <- yvals[c(2, 3, 3, 20, 21, 20)]
    xvals[2] <- xvals[2] + (10 * xvals[2] / 1e6)
    xvals[3] <- xvals[3] - (10 * xvals[3] / 1e6)
    xvals[6] <- xvals[6] + (12 * xvals[6] / 1e6)

    res <- matchWithPpm(xvals, yvals)
    expect_true(is.list(res))
    expect_equal(length(res), length(xvals))
    expect_equal(res[[1]], 2L)
    expect_equal(res[[2]], integer())
    expect_equal(res[[3]], integer())
    expect_equal(res[[4]], 20L)


    res <- matchWithPpm(xvals, yvals, ppm = 10)
    expect_equal(res[[1]], c(2L, 123L, 124L))
    expect_equal(res[[2]], 3L)
    expect_equal(res[[3]], 3L)

    res <- matchWithPpm(xvals, yvals, ppm = 20)
    expect_equal(res[[1]], c(2L, 123L, 124L, 125L))
})

test_that("adducts works", {
    all <- CompoundDb:::ADDUCTS
    res <- adducts()
    expect_equal(all, res)
    res <- adducts(polarity = -4)
    expect_equal(res, all[all$charge < 0, ])
    res <- adducts(pattern = "H+")
    expect_equal(res, all[grep("H+", all$name), ])
    res <- adducts(pattern = "H+", fixed = TRUE)
    expect_equal(res, all[grep("H+", all$name, fixed = TRUE), ])
    res <- adducts(name = "[M+H]+")
    expect_true(nrow(res) == 1)
})

test_that("mass2mz works", {
    library(testthat)
    library(CompoundDb)
    masses <- c(75.032028409, 105.042595)
    names(masses) <- c("Glycine", "Serine")
    ## [M+H]+
    res <- mass2mz(masses, adduct = "[M+H]+")
    expect_true(is.list(res))
    expect_equal(names(res), names(masses))
    expect_true(all(lengths(res) == 1))
    expect_equal(unname(res[[1]]), 76.039304409)
    expect_equal(unname(res[[2]]), 106.049871)

    ## more complicated stuff - cross checked with HMDB
    res <- mass2mz(masses, adduct = "[2M+K]+")
    expect_true(all.equal(unname(res[[1]]), 189.0272, tolerance = 1e-5))
    res <- mass2mz(masses, adduct = "[M+2H]2+")
    expect_true(all.equal(unname(res[[1]]), 38.5233, tolerance = 1e-5))

    ## [M-H]-
    expct <- c(74.024752409, 104.035319)
    res <- mass2mz(masses, adduct = "[M-H]-")
    expect_equal(unname(res[[1]]), 74.024752409)
    expect_equal(unname(res[[2]]), 104.035319)

    ## All of them
    res <- mass2mz(masses)

    ## The adduct m/z for Phosphorylcholine taken from HMDB: commented adducts
    ## are not in Jan's list.
    hmdb_mz <- rbind(`[M+H]+` = c(1, 185.0811),
                     ## `[M-2H2O+H]+` = c(1, 149.0612),
                     ## `[M-H2O+H]+` = c(1, 167.0712),
                     ## `[M-H2O+NH4]+` = c(1, 184.0966),
                     ## `[M+Li]+` = c(1, 191.0899),
                     `[M+NH4]+` = c(1, 202.1077),
                     `[M+Na]+` = c(1, 207.0631),
                     ## `[M+CH3OH+H]+` = c(1, 217.1074),
                     `[M+K]+` = c(1, 223.0370),
                     ## `[M+ACN+H]+` = c(1, 226.1077),
                     `[M+2Na-H]+` = c(1, 229.0450),
                     ## `[M+IsoProp+H]+` = c(1, 245.1392),
                     ## `[M+ACN+Na]+` = c(1, 248.0896),
                     `[M+2K-H]+` = c(1, 260.9929),
                     ## `[M+DMSO+H]+` = c(1, 263.0951),
                     ## `[M+2ACN+H]+` = c(1, 267.1342),
                     ## `[M+IsoProp+Na+H]+` = c(1, 268.1290),
                     `[2M+H]+` = c(1, 369.1550),
                     ## `[2M+NH4]+` = c(1, 386.1816),
                     `[2M+Na]+` = c(1, 391.1370),
                     ## `[2M+3H2O+2H]+` = c(1, 396.1709),
                     `[2M+K]+` = c(1, 407.1109),
                     ## `[2M+ACN+H]+` = c(1, 410.1816),
                     ## `[2M+ACN+Na]+` = c(1, 432.1635),
                     `[M+2H]2+` = c(2, 93.0442),
                     ## `[M+H+NH4]2+` = c(2, 101.5575),
                     `[M+H+Na]2+` = c(2, 104.0352),
                     `[M+H+K]2+` = c(2, 112.0222),
                     ## `[M+ACN+2H]2+`	= c(2, 113.5575),
                     `[M+2Na]2+` = c(2, 115.0262),
                     ## `[M+2ACN+2H]2+` = c(2, 134.0708),
                     ## `[M+3ACN+2H]2+` = c(2, 154.5840),
                     `[M+3H]3+` = c(3, 62.3652),
                     ## `[M+Na+2H]3+` = c(3, 69.6926),
                     ## `[M+2Na+H]3+` = c(3, 77.1242),
                     ## `[M+3Na]3+` = c(3, 84.3472),
                     `[M-H]-` = c(-1, 183.0666),
                     `[M-H-H2O]-` = c(-1, 165.0555),
                     ## `[M+F]-` = c(-1, 203.0723),
                     `[M-2H+Na]-` = c(-1, 205.0485),
                     `[M+Cl]-` = c(-1, 219.0433),
                     `[M-2H+K]-` = c(-1, 221.0225),
                     ## `[M+FA-H]-` = c(-1, 229.0721),
                     ## `[M+HAc-H]-` = c(-1, 243.0877),
                     ## `[M+Br]-` = c(-1, 262.9928),
                     ## `[M+TFA-H]-` = c(-1, 297.0595),
                     `[2M-H]-` = c(-1, 367.1405),
                     ## `[2M+FA-H]-` = c(-1, 413.1459),
                     ## `[2M+HAc-H]-` = c(-1, 427.1616),
                     `[3M-H]-` = c(-1, 551.2143),
                     `[M-2H]2-` = c(-2, 91.0297),
                     `[M-3H]3-` = c(-3, 60.3507))
    mass <- 184.073869485
    mzs <- CompoundDb:::mass2mz(mass, adduct = rownames(hmdb_mz))
    expect_true(all.equal(mzs[[1]], hmdb_mz[, 2], tolerance = 1e-6))
})

test_that("mz2mass works", {
    mass <- 75.032028409

    ## m/z taken from HMDB
    res <- mz2mass(76.0393, adduct = "[M+H]+")
    expect_equal(names(res[[1]]), "[M+H]+")
    expect_true(all.equal(unname(res[[1]]), mass, tolerance = 1e-6))
    res <- mz2mass(120.0032, adduct = "[M+2Na-H]+")
    expect_true(all.equal(unname(res[[1]]), mass, tolerance = 1e-6))
    res <- mz2mass(151.0713, adduct = "[2M+H]+")
    expect_true(all.equal(unname(res[[1]]), mass, tolerance = 1e-6))
    res <- mz2mass(38.5233, adduct = "[M+2H]2+")
    expect_true(all.equal(unname(res[[1]]), mass, tolerance = 1e-6))

    tmp <- mass2mz(mass, adduct = "[2M+2H]2+")
    res <- mz2mass(tmp, adduct = "[2M+2H]2+")
    expect_equal(unname(res[[1]]), mass)
})

test_that(".annotate_adduct works", {
    x <- c(105.0546, 209.0851, 97.07362, 100000)
    names(x) <- c("a", "b", "c", "d")
    ## 1: [M+H]+ for formula C4H8O3 (have 2)
    ## 2: [M+Cl]- for C11H14N2
    ## 3: [M+Na]+ for C3H10N2
    ## 4: does not exist.
    cmps <- compounds(cmp_db)
    res <- CompoundDb:::.annotate_adduct_mz(x, cmps = cmps, adduct = "[M+H]+")
    expect_equal(length(x), length(res))
    expect_equal(names(x), names(res))
    expect_equal(nrow(res[[2]]), 1)
    expect_equal(nrow(res[[3]]), 1)
    expect_equal(nrow(res[[1]]), 2)
    expect_true(all(res[[1]]$formula == "C4H8O3"))
    expect_true(all(is.na(res[[2]]$formula)))
    expect_true(all(is.na(res[[3]]$formula)))
    expect_true(all(is.na(res[[4]]$formula)))

    res <- CompoundDb:::.annotate_adduct_mz(x, cmps = cmps,
                                            adduct = c("[M+H]+", "[M+Cl]-"))
    expect_equal(length(x), length(res))
    expect_equal(nrow(res[[2]]), 1)
    expect_equal(nrow(res[[3]]), 1)
    expect_equal(nrow(res[[1]]), 2)
    expect_true(all(res[[1]]$formula == "C4H8O3"))
    expect_true(all(res[[2]]$formula == "C11H14N2"))
    expect_true(all(is.na(res[[3]]$formula)))
    expect_true(all(is.na(res[[4]]$formula)))

    res <- CompoundDb:::.annotate_adduct_mz(x, cmps = cmps,
                                            adduct = c("[M+H]+", "[M+Cl]-",
                                                       "[M+Na]+"))
    expect_equal(length(x), length(res))
    expect_equal(nrow(res[[2]]), 1)
    expect_equal(nrow(res[[3]]), 1)
    expect_equal(nrow(res[[1]]), 2)
    expect_true(all(res[[1]]$formula == "C4H8O3"))
    expect_true(all(res[[2]]$formula == "C11H14N2"))
    expect_true(all(res[[3]]$formula == "C3H10N2"))
    expect_true(all(is.na(res[[4]]$formula)))

    ## Only empty.
    res <- CompoundDb:::.annotate_adduct_mz(x[4], cmps)
    expect_true(all(is.na(res[[1]]$formula)))
})

test_that("annotateMz,ANY,ANY works", {
    mzs <- c(105.0546, 75.09168, 127.0365, 113.04756, 113.1210)
    ## 1: [M+H+]+ for HMDB0000008 and HMDB0000011
    ## 2: [M+H+]+ for HMDB0000002
    ## 3: [M+Na+]+ for HMDB0000008 and HMDB0000011
    ## 4: [M+K+]+ for HMDB0000002
    ## 5: unknown

    expect_error(annotateMz(), "unable to find")

    ## numeric, data.frame
    expect_error(annotateMz(mzs, data.frame()), "Required column")
    res <- annotateMz(mzs, cmps, ppm = 5)
    expect_true(length(res) == length(mzs))
    expect_true(all(unlist(lapply(res, class)) == class(cmps)))
    expect_equal(res[[1]]$compound_id[1], "HMDB0000008")
    expect_equal(res[[1]]$compound_id[2], "HMDB0000011")
    expect_equal(res[[2]]$compound_id[1], "HMDB0000002")
    expect_equal(res[[3]]$compound_id[1], "HMDB0000008")
    expect_equal(res[[4]]$compound_id[1], "HMDB0000002")
    expect_equal(res[[4]]$adduct[1], "[M+K]+")

    res <- annotateMz(mzs, as.data.frame(cmps), adduct = "[M+H]+")
    expect_true(all(unlist(lapply(res, class)) == "data.frame"))
    expect_equal(res[[1]]$compound_id, c("HMDB0000008", "HMDB0000011"))
    expect_equal(res[[2]]$compound_id, c("HMDB0000002"))
    expect_true(is.na(res[[3]]$compound_id))
    expect_true(is.na(res[[4]]$compound_id))
    expect_true(is.na(res[[5]]$compound_id))

    ## data.frame, data.frame
    mzs_df <- data.frame(id = letters[1:5], mzmed = mzs,
                         stringsAsFactors = FALSE)
    expect_error(annotateMz(mzs_df, cmps), "Column")
    res <- annotateMz(mzs_df, cmps, mzcol = "mzmed", adduct = "[M+K]+")
    expect_true(is.data.frame(res))
    expect_true(all(colnames(res) == c(colnames(mzs_df), colnames(cmps),
                                       "adduct", "ppm")))
    expect_equal(res$compound_id, c(NA, NA, NA, "HMDB0000002", NA))
    res <- annotateMz(mzs_df, cmps, mzcol = "mzmed",
                      adduct = c("[M+H]+", "[M+Na]+"))
    expect_equal(res$id, c("a", "a", "b", "c", "c", "d", "e"))
    expect_equal(res$compound_id, c("HMDB0000008", "HMDB0000011", "HMDB0000002",
                                    "HMDB0000008", "HMDB0000011", NA, NA))

    ## numeric, CompDb
    res <- annotateMz(mzs, cmp_db, adduct = "[M+H]+")
    expect_true(length(res) == length(mzs))
    expect_equal(res[[1]]$compound_id, c("HMDB0000008", "HMDB0000011"))
    expect_equal(res[[2]]$compound_id, c("HMDB0000002"))

    ## data.frame, CompDb
    expect_error(annotateMz(mzs_df, cmp_db, adduct = "[M+H]+"), "Column")
    res <- annotateMz(mzs_df, cmp_db, adduct = "[M+H]+", mzcol = "mzmed")
    expect_true(is.data.frame(res))
    expect_equal(res$id, c("a", "a", "b", "c", "d", "e"))
})
