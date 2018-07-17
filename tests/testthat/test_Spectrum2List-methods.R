test_that("$,Spectrum2List and $<-,Spectrum2List works", {
    spl <- spl_
    expect_equal(spl$spectrum_id, mcols(spl)$spectrum_id)
    expect_equal(spl$splash, mcols(spl)$splash)

    spl$other <- 4
    expect_true(any(colnames(mcols(spl)) == "other"))
    expect_equal(spl$other, c(4, 4, 4, 4))
})

test_that("mz, intensity, rtime work", {
    spl <- spl_
    expect_true(length(mz(spl)) == 4)
    expect_equal(mz(spl[[2]]), mz(spl)[[2]])

    expect_true(length(intensity(spl)) == 4)
    expect_equal(intensity(spl[[2]]), intensity(spl)[[2]])
    
    expect_true(length(rtime(spl)) == 4)
    expect_equal(rtime(spl[[2]]), rtime(spl)[[2]])

    spl <- c(spl, Spectrum2List(new("Spectrum2")))
    expect_true(lengths(mz(spl))[5] == 0)
    expect_true(lengths(intensity(spl))[5] == 0)
    expect_equal(rtime(spl), rep(NA_real_, 5))

    ## Put names on it.
    names(spl) <- c("a", "b", "c", "d", "e")
    expect_equal(names(rtime(spl)), c("a", "b", "c", "d", "e"))
    expect_equal(names(mz(spl)), c("a", "b", "c", "d", "e"))
    expect_equal(names(intensity(spl)), c("a", "b", "c", "d", "e"))
})
