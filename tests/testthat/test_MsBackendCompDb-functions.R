test_that(".valid_dbcon MsBackendCompDb works", {
    expect_equal(.valid_dbcon(.dbconn(cdb)), NULL)
    expect_match(.valid_dbcon(.dbconn(cmp_db)), "no MS/MS")
})

test_that("MsBackendCompDb works", {
    res <- MsBackendCompDb()
    expect_true(is(res, "MsBackendCompDb"))
    expect_true(is.null(.dbconn(res)))
})

test_that(".map_spectraVariables_to_sql works", {
    res <- .map_spectraVariables_to_sql("msLevel")
    expect_equal(res, "ms_level")
    res <- .map_spectraVariables_to_sql(c("mz", "collisionEnergy", "ms_level"))
    expect_equal(res, c("mz", "collision_energy", "ms_level"))
})

test_that(".map_sql_to_spectraVariables works", {
    res <- .map_sql_to_spectraVariables("ms_level")
    expect_equal(res, "msLevel")
    res <- .map_sql_to_spectraVariables(c("precursor_mz", "intensity"))
    expect_equal(res, c("precursorMz", "intensity"))
})

test_that(".fetch_peaks works", {
    res <- .fetch_peaks(MsBackendCompDb())
    expect_true(is.data.frame(res))
    expect_true(nrow(res) == 0)
    expect_true(all(colnames(res) %in% c("spectrum_id", "mz", "intensity")))

    be <- backendInitialize(MsBackendCompDb(), cdb)
    res <- .fetch_peaks(be)
    expect_true(is.data.frame(res))
    expect_true(nrow(res) > 0)
    expect_true(all(colnames(res) %in% c("spectrum_id", "mz", "intensity")))
})

test_that(".peaks_data works", {
    res <- .peaks_data(MsBackendCompDb())
    expect_true(length(res) == 0)
    expect_true(is.list(res))

    be <- backendInitialize(MsBackendCompDb(), cdb)
    pks <- .fetch_peaks(be)
    res <- .peaks_data(be)
    expect_true(length(res) > 1)
    expect_true(is.matrix(res[[1L]]))
})

test_that(".spectra_data works", {
    res <- .spectra_data(MsBackendCompDb())
    expect_true(nrow(res) == 0)
    expect_true(is(res, "DataFrame"))
    expect_equal(colnames(res), spectraVariables(MsBackendCompDb()))

    be <- backendInitialize(MsBackendCompDb(), cdb)
    res <- .spectra_data(be)
    expect_true(nrow(res) == length(be))
    expect_equal(colnames(res), spectraVariables(be))
    expect_equal(be@spectraIds, as.character(res$spectrum_id))

    expect_true(is(res$mz, "NumericList"))

    res_2 <- .spectra_data(be, "mz")
    expect_true(is(res_2, "DataFrame"))
    expect_true(colnames(res_2) == "mz")
    expect_equal(res$mz, res_2$mz)

    res_2 <- .spectra_data(be, c("name", "compound_id"))
    expect_true(is(res_2, "DataFrame"))
    expect_true(all(colnames(res_2) == c("name", "compound_id")))
    expect_equal(nrow(res_2), length(be))

    be <- backendInitialize(MsBackendCompDb(), cmp_spctra_db)
    res_2 <- .spectra_data(be, c("name", "compound_id"))
    expect_true(is(res_2, "DataFrame"))
    expect_true(all(colnames(res_2) == c("name", "compound_id")))
    expect_equal(nrow(res_2), length(be))
})

test_that(".available_peaks_variables works", {
    expect_equal(.available_peaks_variables(MsBackendCompDb()), character())
    res <- .available_peaks_variables(Spectra(cmp_spctra_db)@backend)
    expect_equal(res, c("mz", "intensity"))
})
