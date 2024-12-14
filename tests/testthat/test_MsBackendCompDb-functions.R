test_that(".valid_dbcon MsBackendCompDb works", {
    expect_equal(.valid_dbcon(.dbconn(cdb)), NULL)
    expect_match(.valid_dbcon(.dbconn(cmp_db)), "no MS/MS")
    expect_match(.valid_dbcon("a"), "connection to a database")
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
    tmp <- data.frame(
        spectrum_id = c(1, 1, 1,
                        2, 2, 2,
                        3, 3, 3,
                        5, 5),
        mz = c(1.2, 1.3, 1.4,
               2.4, 2.5, 2.6,
               4.2, 5.3, 5.7,
               6.4, 6.7),
        intensity = c(1, 2, 3,
                      4, 5, 6,
                      7, 8, 9,
                      10, 11))
    ## Missing peaks for spectrum4
    x <- MsBackendCompDb()
    x@spectraIds <- as.character(1:5)
    res <- CompoundDb:::.peaks_data(x, p = tmp)
    expect_true(is.list(res))
    expect_equal(length(res), 5)
    expect_true(all(vapply(res, is.matrix, NA)))
    expect_equal(res[[1]][, "intensity"], 1:3)
    expect_equal(res[[2]][, "intensity"], 4:6)
    expect_equal(res[[3]][, "intensity"], 7:9)
    expect_equal(res[[4]][, "intensity"], numeric())
    expect_equal(res[[5]][, "intensity"], 10:11)

    ## Arbitrary order and duplication
    x@spectraIds <- as.character(c(4, 2, 5, 2))
    res <- CompoundDb:::.peaks_data(x, p = tmp)
    expect_true(is.list(res))
    expect_equal(length(res), 4)
    expect_true(all(vapply(res, is.matrix, NA)))
    expect_equal(res[[1]][, "intensity"], numeric())
    expect_equal(res[[2]][, "intensity"], 4:6)
    expect_equal(res[[3]][, "intensity"], 10:11)
    expect_equal(res[[4]][, "intensity"], 4:6)

    ## only intensity
    res <- CompoundDb:::.peaks_data(x, columns = "intensity",
                                    p = tmp[, c("spectrum_id", "intensity")])
    expect_true(is.list(res))
    expect_equal(length(res), 4)
    expect_true(all(vapply(res, is.matrix, NA)))
    expect_equal(res[[1]][, "intensity"], numeric())
    expect_equal(res[[2]][, "intensity"], 4:6)
    expect_equal(res[[3]][, "intensity"], 10:11)
    expect_equal(res[[4]][, "intensity"], 4:6)

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

    expect_error(
        .spectra_data(be, columns = c("name", "compound_id", "not_exists")),
        "not available")
})

test_that(".available_peaks_variables works", {
    expect_equal(.available_peaks_variables(MsBackendCompDb()), character())
    res <- .available_peaks_variables(Spectra(cmp_spctra_db)@backend)
    expect_equal(res, c("mz", "intensity"))
})
