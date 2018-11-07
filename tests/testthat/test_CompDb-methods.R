
test_that("show,CompDb works", {
    expect_output(show(cmp_db))
    db <- new("CompDb")
    expect_output(show(cmp_db))
})

test_that("dbconn,CompDb works", {
    expect_true(!is.null(dbconn(cmp_db)))
    expect_true(is(dbconn(cmp_db), "DBIConnection"))
})

test_that("spectra,CompDb works", {
    expect_error(spectra(cmp_db))
    res <- spectra(cmp_spctra_db, return.type = "Spectra")
    expect_true(is(res, "Spectra"))
    expect_true(all(msLevel(res) == 2))
    expect_true(length(res) == 4)
    res <- spectra(cmp_spctra_db, return.type = "data.frame")
    expect_true(is.data.frame(res))
    ## mzs <- split(res$mz, res$spectrum_id)
    mzs <- res$mz
    expect_true(all(!vapply(mzs, is.unsorted, logical(1))))
    res <- spectra(cmp_spctra_db, return.type = "tibble")
    expect_true(is(res, "tbl"))

    ## Specific columns
    res <- spectra(cmp_spctra_db, columns = c("compound_name", "inchi"))
    expect_true(all(c("spectrum_id", "compound_name", "inchi") %in%
                    colnames(mcols(res))))
    expect_true(length(res) == 4)

    ## With filter.
    res <- spectra(cmp_spctra_db, filter = ~ compound_id == "HMDB0000001")
    expect_true(length(res) == 2)
    res_all <- spectra(cmp_spctra_db)
    expect_equal(res[[which(mcols(res)$spectrum_id == "1")]],
                 res_all[[which(mcols(res_all)$spectrum_id == "1")]])
    expect_equal(res[[which(mcols(res)$spectrum_id == "2")]],
                 res_all[[which(mcols(res_all)$spectrum_id == "2")]])
})
