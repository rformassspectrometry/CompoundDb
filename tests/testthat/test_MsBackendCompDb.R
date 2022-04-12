test_that("backendInitialize,MsBackendCompDb works", {
    res <- backendInitialize(MsBackendCompDb(), cmp_spctra_db)
    expect_true(is(res, "MsBackendCompDb"))
    expect_true(length(res) == 4)
    expect_true(!is.null(res@dbcon))

    expect_error(backendInitialize(MsBackendCompDb(), 4), "'CompDb'")
    expect_error(backendInitialize(MsBackendCompDb(), cmp_db), "no MS/MS")

    res <- backendInitialize(MsBackendCompDb(), cmp_spctra_db,
                             filter = ~ compound_id == "HMDB0000008")
    expect_true(length(res) == 0)
    expect_equal(res@spectraIds, character())

    res <- backendInitialize(MsBackendCompDb(), cmp_spctra_db,
                             filter = ~ compound_id == "HMDB0000001")
    expect_true(length(res) == 2)
    expect_true(all(res$compound_id == "HMDB0000001"))

    res <- backendInitialize(MsBackendCompDb(), cmp_spctra_db,
                             filter = ~ compound_id == "bla")
    expect_true(length(res) == 0)
})

test_that("peaksData,MsBackendCompDb works", {
    be <- MsBackendCompDb()
    res <- peaksData(be)
    expect_true(length(res) == 0)
    expect_true(is(res, "list"))

    expect_error(peaksData(be, columns = c("other")))

    be <- Spectra(cmp_spctra_db)@backend
    res <- peaksData(be)
    expect_true(is(res, "list"))
    expect_true(length(res) == length(be))
    expect_true(is.matrix(res[[1]]))
    expect_equal(colnames(res[[2]]), c("mz", "intensity"))
    res_2 <- peaksData(be, c("intensity", "mz"))
    expect_equal(colnames(res_2[[2]]), c("intensity", "mz"))
    expect_equal(res[[2]][, 1], res_2[[2]][, 2])
    res_2 <- peaksData(be, c("intensity"))
    expect_equal(colnames(res_2[[1]]), "intensity")
    expect_equal(res[[2]][, 2], res_2[[2]][, 1])

    be <- be[c(2, 4, 2)]
    res_2 <- peaksData(be)
    expect_equal(res_2, res[c(2, 4, 2)])
})

test_that("peaksVariables,MsBackendCompDb works", {
    expect_equal(peaksVariables(MsBackendCompDb()), character())

    res <- peaksVariables(Spectra(cmp_spctra_db)@backend)
    expect_equal(res, c("mz", "intensity"))
})

test_that("dataStorage,MsBackendCompDb works", {
    be <- MsBackendCompDb()
    res <- dataStorage(be)
    expect_equal(res, character())

    be <- backendInitialize(MsBackendCompDb(), cdb)
    res <- dataStorage(be)
    expect_equal(res, rep("<db>", length(be)))
})

test_that("intensity,intensity<-,MsBackendCompDb works", {
    be <- MsBackendCompDb()
    res <- be$intensity
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == length(be))

    be <- Spectra(cmp_spctra_db)@backend
    res <- be$intensity
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == length(be))

    expect_error(intensity(be) <- res, "not")
})

test_that("mz,mz<-,MsBackendCompDb works", {
    be <- MsBackendCompDb()
    res <- be$mz
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == length(be))

    be <- Spectra(cmp_spctra_db)@backend
    res <- be$mz
    expect_true(is(res, "NumericList"))
    expect_true(length(res) == length(be))

    expect_error(mz(be) <- res, "not replace")
})

test_that("spectraData,spectraData<-,MsBackendCompDb works", {
    be <- MsBackendCompDb()
    res <- spectraData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 0)

    be <- backendInitialize(MsBackendCompDb(), cmp_spctra_db)
    res <- spectraData(be)
    expect_true(is(res, "DataFrame"))
    expect_true(nrow(res) == 4)
    expect_equal(res$mz, be$mz)
    expect_equal(res$intensity, be$intensity)

    be <- be[c(3, 4, 2, 1)]
    res_2 <- spectraData(be, c("compound_id", "mz", "polarity"))
    expect_true(all(c("compound_id", "mz", "polarity") == colnames(res_2)))
    expect_equal(res_2$mz, be$mz)
    expect_equal(res_2$mz, res$mz[c(3, 4, 2, 1)])

    expect_error(spectraData(be, "sorry"), "not available")

    be <- backendInitialize(MsBackendCompDb(), cmp_spctra_db)
    res <- spectraData(be)
    be <- be[c(3, 1, 2, 2, 2)]
    expect_equal(be@spectraIds, as.character(c(3, 1, 2, 2, 2)))
    res_2 <- spectraData(be)
    expect_equal(res_2$spectrum_id, res$spectrum_id[c(3, 1, 2, 2, 2)])
    expect_equal(res_2$intensity, res$intensity[c(3, 1, 2, 2, 2)])
})

test_that("spectraNames,MsBackendCompDb works", {
    be <- backendInitialize(MsBackendCompDb(), cdb)
    res <- spectraNames(be)
    expect_equal(res, be@spectraIds)
})

test_that("$<-,MsBackendCompDb works", {
    be <- backendInitialize(MsBackendCompDb(), cdb)
    be$polarity <- 0L
    expect_true(any(colnames(be@localData) == "polarity"))
    expect_equal(be$polarity, rep(0L, length(be)))

    be$new_col <- "a"
    expect_equal(be$new_col, rep("a", length(be)))

    expect_error(be$spectrum_id <- "a", "not")
    expect_error(be$mz <- be$mz, "not supported")
})

test_that("[,MsBackendCompDb works", {
    be <- backendInitialize(MsBackendCompDb(), cmp_spctra_db)
    res <- be[c(2, 4)]
    expect_true(length(res) == 2)
    expect_equal(res$polarity, be$polarity[c(2, 4)])
    expect_equal(res@spectraIds, be@spectraIds[c(2, 4)])
    expect_equal(res$mz, be$mz[c(2, 4)])

    ## Arbitrary order and duplicates
    be$my_index <- seq_along(be)
    idx <- c(3, 1, 2, 1, 1, 1, 2)
    res <- be[idx]
    expect_equal(res$polarity, be$polarity[idx])
    expect_equal(res@spectraIds, be@spectraIds[idx])
    expect_equal(res$my_index, be$my_index[idx])
    expect_equal(res$intensity, be$intensity[idx])
    expect_equal(res$mz, be$mz[idx])
    expect_equal(res$compound_id, be$compound_id[idx])
})

test_that("show,MsBackendCompDb doesn't break", {
    be <- MsBackendCompDb()
    show(be)
    be <- Spectra(cmp_spctra_db)@backend
    show(be)
})
