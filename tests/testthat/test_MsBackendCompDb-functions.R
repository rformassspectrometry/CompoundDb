test_that(".valid_spectra_data_required_columns works", {
    df <- DataFrame()
    expect_null(.valid_spectra_data_required_columns(df))
    df <- DataFrame(msLevel = 1L)
    expect_match(.valid_spectra_data_required_columns(df),
                 "Required column")
    df$dataStorage <- "some"
    expect_null(.valid_spectra_data_required_columns(df))
})

test_that(".valid_ms_backend_dbcon works", {
    res <- .valid_ms_backend_dbcon(cmp_db@dbcon)
    expect_true(length(res) == 1)
    res <- .valid_ms_backend_dbcon(cmp_spctra_db@dbcon)
    expect_true(length(res) == 0)
})

test_that("MsBackendCompDb works", {
    res <- MsBackendCompDb()
    expect_true(is(res, "MsBackendCompDb"))
    expect_true(is.null(res@dbcon))
})

test_that(".peaks works", {
    be <- Spectra(cmp_spctra_db)@backend
    pks <- .peaks(be)
    expect_equal(names(pks), as.character(be$spectrum_id))
    expect_true(is.list(pks))
    expect_true(is.matrix(pks[[1]]))
    expect_true(all(colnames(pks[[2]]) == c("mz", "intensity")))

    be <- be[c(3, 4, 2)]
    pks <- .peaks(be)
    expect_equal(names(pks), as.character(be$spectrum_id))

    mzs <- .peaks(be, column = "mz")
    expect_true(is.list(mzs))
    expect_true(is.numeric(mzs[[2]]))
    expect_identical(mzs, lapply(pks, function(z) z[, 1]))

    ints <- .peaks(be, column = "intensity")
    expect_true(is.list(ints))
    expect_true(is.numeric(ints[[2]]))
    expect_identical(ints, lapply(pks, function(z) z[, 2]))
})
