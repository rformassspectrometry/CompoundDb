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

test_that("precursor* work", {
    sp1 <- new("Spectrum2", precursorMz = 123.3, precursorCharge = 1L,
               precursorIntensity = 1234.4)
    sp2 <- new("Spectrum2", precursorMz = NA_real_, precursorCharge = integer(),
               precursorIntensity = NA_real_, precScanNum = 34L)
    spl <- Spectrum2List(sp1, sp2)
    expect_equal(precursorMz(spl), c(123.3, NA))
    expect_equal(precursorCharge(spl), c(1L, NA))
    expect_true(is.integer(precursorCharge(spl)))
    expect_equal(precursorIntensity(spl), c(1234.4, NA))
    expect_equal(precScanNum(spl), c(NA_integer_, 34L))
    
    expect_equal(precursorMz(spl_), rep(NA_real_, length(spl_)))
    expect_equal(precursorCharge(spl_), rep(NA_integer_, length(spl_)))
    expect_equal(precursorIntensity(spl_), rep(NA_real_, length(spl_)))
    expect_true(is.integer(precScanNum(spl_)))
})

test_that("acquisitionNum and scanIndex work", {
    sp1 <- new("Spectrum2", acquisitionNum = 2L, scanIndex = 1L)
    sp2 <- new("Spectrum2", acquisitionNum = 4L)
    spl <- Spectrum2List(sp1, sp2)
    expect_identical(acquisitionNum(spl), c(2L, 4L))
    expect_identical(scanIndex(spl), c(1L, NA_integer_))
        
    expect_equal(acquisitionNum(spl_), rep(NA_integer_, length(spl_)))
    expect_equal(scanIndex(spl_), rep(NA_integer_, length(spl_)))
    expect_true(is.integer(acquisitionNum(spl_)))
    expect_true(is.integer(scanIndex(spl_)))
})

test_that("peaksCount, msLevel, tic and ionCount work", {
    sp1 <- new("Spectrum2", msLevel = 3L, tic = 5)
    sp2 <- new("Spectrum2")
    spl <- Spectrum2List(sp1, sp2)

    expect_true(is.integer(peaksCount(spl)))
    expect_equal(peaksCount(spl), c(0, 0))
    expect_true(is.integer(msLevel(spl)))
    expect_equal(msLevel(spl), c(3, 2))
    expect_true(is.numeric(tic(spl)))
    expect_equal(tic(spl), c(5, 0))
    expect_true(is.numeric(ionCount(spl)))
    expect_equal(ionCount(spl), c(0, 0))

    expect_equal(peaksCount(spl_), c(5, 7, 33, 10))
    expect_true(all(msLevel(spl_) == 2))
    expect_true(all(is.na(tic(spl_))))
    expect_equal(ionCount(spl_), unlist(lapply(intensity(spl_), sum)))
})

test_that("collisionEnergy works", {
    sp1 <- new("Spectrum2")
    sp2 <- new("Spectrum2", collisionEnergy = 23.3)
    spl <- Spectrum2List(sp1, sp2)

    expect_true(is.numeric(collisionEnergy(spl)))
    expect_equal(collisionEnergy(spl), c(NA, 23.3))

    expect_true(is.numeric(collisionEnergy(spl_)))
    expect_equal(collisionEnergy(spl_), c(10, 25, NA, 20))
})

test_that("fromFile and polarity work", {
    sp1 <- new("Spectrum2", polarity = 1L, fromFile = 5L)
    sp2 <- new("Spectrum2", fromFile = 3L)
    spl <- Spectrum2List(sp1, sp2)

    expect_true(is.integer(fromFile(spl)))
    expect_equal(fromFile(spl), c(5, 3))
    expect_true(is.integer(polarity(spl)))
    expect_equal(polarity(spl), c(1, NA))

    expect_true(is.integer(fromFile(spl_)))
    expect_true(all(is.na(fromFile(spl_))))
    expect_true(is.integer(polarity(spl_)))
    expect_equal(polarity(spl_), c(1, 1, 1, 0))
})

test_that("smoothed, isEmpty, centroided and isCentroided work", {
    sp1 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(4, 2, 4, 5))
    sp2 <- new("Spectrum2", centroided = TRUE, smoothed = TRUE)
    spl <- Spectrum2List(sp1, sp2)

    expect_true(is.logical(smoothed(spl)))
    expect_equal(smoothed(spl), c(NA, TRUE))
    expect_true(is.logical(isEmpty(spl)))
    expect_equal(isEmpty(spl), c(FALSE, TRUE))
    expect_true(is.logical(centroided(spl)))
    expect_equal(centroided(spl), c(NA, TRUE))
    expect_true(is.logical(isCentroided(spl)))
    expect_equal(isCentroided(spl), c(NA, NA))

    expect_true(is.logical(smoothed(spl_)))
    expect_equal(smoothed(spl_), rep(NA, length(spl_)))
    expect_true(is.logical(isEmpty(spl_)))
    expect_true(all(!isEmpty(spl_)))
    expect_true(is.logical(centroided(spl_)))
    expect_true(all(is.na(centroided(spl_))))
    expect_true(is.logical(isCentroided(spl_)))
    expect_equal(isCentroided(spl_), c(NA, NA, TRUE, NA))
})
