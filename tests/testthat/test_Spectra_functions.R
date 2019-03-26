test_that(".spectrum_has_mz and hasMz work", {
    sp1 <- new("Spectrum1", mz = c(23.231, 123.43, 255.231, 432.0952),
               intensity = c(123, 3432, 45432, 423))
    sp2 <- new("Spectrum1", mz = c(123.099, 344.531, 453.2313),
               intensity = c(231, 431, 413))
    sp3 <- new("Spectrum1", mz = c(123.1001, 343.4321, 432.0921),
               intensity = c(542, 4524, 32))
    spl <- Spectra(sp1, sp2, sp3)

    expect_error(hasMz(sp1))
    expect_error(hasMz(3, mz = 34))

    mzs <- c(123.1, 432.0931)
    expect_true(hasMz(sp1, mzs))
    expect_false(hasMz(sp1, mzs, which = "all"))
    expect_false(hasMz(sp1, mzs, ppm = 0))

    res <- hasMz(spl, mzs)
    expect_equal(unname(res), c(TRUE, FALSE, TRUE))
    names(spl) <- c("a", "b", "c")
    res <- hasMz(spl, mzs)
    expect_equal(res, c(a = TRUE, b = FALSE, c = TRUE))
    res <- hasMz(spl, mzs, which = "all")
    expect_equal(res, c(a = FALSE, b = FALSE, c = TRUE))
})
