test_that(".has_mz and hasMz work", {
    mzs <- c(23.231, 123.43, 255.231, 432.0952)
    expect_true(.has_mz(mzs, c(123.1, 432.0931)))
    expect_false(.has_mz(mzs, c(123.1, 432.0931), which = "all"))
    expect_false(.has_mz(mzs, c(123.1, 432.0931), ppm = 0))
    expect_true(.has_mz(mzs, c(123.425, 432.0931), which = "all", ppm = 50))

    expect_true(.has_mz(mzs, c(432.0931)))
    expect_true(.has_mz(mzs, c(432.0931), which = "all"))

    ## On a full spectra...
    sps <- Spectra(cmp_spctra_db)
    mzs <- mz(sps)
    res <- vapply(mzs, .has_mz, logical(1), mz = 10)
    expect_true(all(res == FALSE))
    res <- vapply(mzs, .has_mz, logical(1), mz = 71.2)
    expect_true(any(res))
})
