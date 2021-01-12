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
    expect_equal(nrow(res[[2]]), 0)
    expect_equal(nrow(res[[3]]), 0)
    expect_equal(nrow(res[[1]]), 2)
    expect_true(all(res[[1]]$formula == "C4H8O3"))
    expect_true(all(is.na(res[[2]]$formula)))
    expect_true(all(is.na(res[[3]]$formula)))
    expect_true(all(is.na(res[[4]]$formula)))

    res <- CompoundDb:::.annotate_adduct_mz(x, cmps = cmps,
                                            adduct = c("[M+H]+", "[M+Cl]-"))
    expect_equal(length(x), length(res))
    expect_equal(nrow(res[[2]]), 1)
    expect_equal(nrow(res[[3]]), 0)
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
    expect_true(nrow(res[[1]]) == 0)
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
    res <- annotateMz(mzs, cmps, ppm = 5,
                      adduct = c("[M+H]+", "[M+Na]+", "[M+K]+"))
    expect_true(length(res) == length(mzs))
    expect_true(all(unlist(lapply(res, class)) == class(cmps)))
    expect_equal(unname(res[[1]]$compound_id[1]), "HMDB0000008")
    expect_equal(unname(res[[1]]$compound_id[2]), "HMDB0000011")
    expect_equal(unname(res[[2]]$compound_id[1]), "HMDB0000002")
    expect_equal(unname(res[[3]]$compound_id[1]), "HMDB0000008")
    expect_equal(unname(res[[4]]$compound_id[1]), "HMDB0000002")
    expect_equal(res[[4]]$.adduct[1], "[M+K]+")

    res <- annotateMz(mzs, as.data.frame(cmps), adduct = "[M+H]+")
    expect_true(all(unlist(lapply(res, class)) == "data.frame"))
    expect_equal(unname(res[[1]]$compound_id), c("HMDB0000008", "HMDB0000011"))
    expect_equal(unname(res[[2]]$compound_id), c("HMDB0000002"))
    expect_true(length(res[[3]]$compound_id) == 0)
    expect_true(length(res[[4]]$compound_id) == 0)
    expect_true(length(res[[5]]$compound_id) == 0)

    ## data.frame, data.frame
    mzs_df <- data.frame(id = letters[1:5], mzmed = mzs,
                         stringsAsFactors = FALSE)
    expect_error(annotateMz(mzs_df, cmps), "Column")
    res <- annotateMz(mzs_df, cmps, mzcol = "mzmed", adduct = "[M+K]+")
    expect_true(is.data.frame(res))
    expect_true(all(colnames(res) == c(colnames(mzs_df), colnames(cmps),
                                       ".adduct", ".difference")))
    expect_equal(unname(res$compound_id), c(NA, NA, NA, "HMDB0000002", NA))
    res <- annotateMz(mzs_df, cmps, mzcol = "mzmed",
                      adduct = c("[M+H]+", "[M+Na]+"))
    expect_equal(res$id, c("a", "a", "b", "c", "c", "d", "e"))
    expect_equal(unname(res$compound_id), c("HMDB0000008", "HMDB0000011",
                                            "HMDB0000002", "HMDB0000008",
                                            "HMDB0000011", NA, NA))

    ## numeric, CompDb
    res <- annotateMz(mzs, cmp_db, adduct = "[M+H]+")
    expect_true(length(res) == length(mzs))
    expect_true(all(is.na(res[[3]])))
    expect_true(all(is.na(res[[4]])))
    expect_true(all(is.na(res[[5]])))
    res <- annotateMz(mzs, cmp_db, columns = c("compound_id", "exactmass"),
                      adduct = "[M+H]+")
    expect_true(length(res) == length(mzs))
    expect_true(all(is.na(res[[3]])))
    expect_true(all(is.na(res[[4]])))
    expect_true(all(is.na(res[[5]])))
    expect_equal(res[[1]]$compound_id, c("HMDB0000008", "HMDB0000011"))
    expect_equal(res[[2]]$compound_id, c("HMDB0000002"))

    ## data.frame, CompDb
    expect_error(annotateMz(mzs_df, cmp_db, adduct = "[M+H]+"), "Column")
    res <- annotateMz(mzs_df, cmp_db, adduct = "[M+H]+", mzcol = "mzmed")
    expect_true(is.data.frame(res))
    expect_equal(res$id, c("a", "a", "b", "c", "d", "e"))
})
