test_that("backendInitialize,MsBackendCompDb works", {
    res <- backendInitialize(MsBackendCompDb(), cmp_spctra_db)
    expect_true(is(res, "MsBackendCompDb"))
    expect_true(length(res) == 4)
    expect_true(!is.null(res@dbcon))

    res <- backendInitialize(MsBackendCompDb(), cmp_spctra_db,
                             columns = c("instrument", "compound_id"))
    expect_true(all(colnames(res@spectraData) == c("compound_id", "instrument",
                                                   "spectrum_id", "dataStorage",
                                                   "dataOrigin")))
    res <- backendInitialize(MsBackendCompDb(), cmp_spctra_db,
                             filter = ~ compound_id == "HMDB0000001")
    expect_true(all(res$compound_id == "HMDB0000001"))
    res <- backendInitialize(MsBackendCompDb(), cmp_spctra_db,
                             filter = ~ compound_id == "HMDB0000001",
                             columns = "polarity")
    expect_true(all(res$compound_id == "HMDB0000001"))
    expect_true(all(colnames(res@spectraData) == c("compound_id", "polarity",
                                                   "spectrum_id", "dataStorage",
                                                   "dataOrigin")))
})

test_that("peaksData,MsBackendCompDb works", {
    be <- MsBackendCompDb()
    res <- peaksData(be)
    expect_true(length(res) == 0)
    expect_true(is(res, "list"))

    be <- Spectra(cmp_spctra_db)@backend
    res <- peaksData(be)
    expect_true(is(res, "list"))
    expect_true(length(res) == length(be))
    expect_true(is.matrix(res[[1]]))
    expect_true(all(colnames(res[[2]]) %in% c("mz", "intensity")))
})

test_that("intensity,intensity<-,MsBackendCompDb works", {
    be <- MsBackendCompDb()
    res <- intensity(be)
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == length(be))

    be <- Spectra(cmp_spctra_db)@backend
    res <- intensity(be)
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == length(be))

    expect_error(intensity(be) <- res, "not support replacing")
})

test_that("mz,mz<-,MsBackendCompDb works", {
    be <- MsBackendCompDb()
    res <- mz(be)
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == length(be))

    be <- Spectra(cmp_spctra_db)@backend
    res <- mz(be)
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == length(be))

    expect_error(mz(be) <- res, "not support replacing")
})

test_that("spectraData,spectraData<-,MsBackendCompDb works", {
    be <- MsBackendCompDb()
    res <- spectraData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0)

    be <- Spectra(cmp_spctra_db)@backend
    res <- spectraData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 4)
    expect_equal(res$mz, mz(be))
    expect_equal(res$intensity, intensity(be))

    be <- be[c(3, 4, 2, 1)]
    res <- spectraData(be, c("compound_id", "mz", "polarity"))
    expect_true(all(c("compound_id", "mz", "polarity") == colnames(res)))
    expect_equal(res$mz, mz(be))

    expect_error(spectraData(be, "sorry"), "column 'sorry' not available")

    ## spectraData<-
    res <- spectraData(be)
    res$compound_id <- c("A", "B", "C", "D")
    res$polarity <- rep(0L, 4)
    expect_warning(spectraData(be) <- res, "Ignoring columns")
})

test_that("$<-,MsBackendCompDb works", {
    be <- Spectra(cmp_spctra_db)@backend
    be$polarity <- rep(1L, 4)
    expect_equal(be$polarity, rep(1L, 4))

    expect_error(be$spectrum_id <- rep("a", 4), "does not support")
    expect_error(be$mz <- be$mz, "does not support")
})

test_that("show,MsBackendCompDb doesn't break", {
    be <- MsBackendCompDb()
    show(be)
    be <- Spectra(cmp_spctra_db)@backend
    show(be)
})
