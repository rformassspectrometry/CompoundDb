
test_that("show,CompDb works", {
    expect_output(show(cmp_db))
    db <- new("CompDb")
    expect_output(show(cmp_db))
})

test_that("dbconn,CompDb works", {
    expect_true(!is.null(dbconn(cmp_db)))
    expect_true(is(dbconn(cmp_db), "DBIConnection"))
})

test_that("Spectra,CompDb works", {
    expect_warning(res <- Spectra(cmp_db), "No spectrum data")
    expect_true(is(res, "Spectra"))

    res <- Spectra(cmp_spctra_db)
    expect_true(is(res, "Spectra"))
    expect_true(length(res) == 4)
    expect_true(all(c("instrument", "predicted") %in% spectraVariables(res)))

    ## filter
    res <- Spectra(cmp_spctra_db, filter = ~ compound_id == "HMDB0000001")
    expect_true(length(res) == 2)

    expect_error(Spectra(cmp_spctra_db, filter = "ad"), "'filter' has to")
    expect_error(Spectra(cmp_spctra_db, filter = ~ gene_name == "b"),
                 "not supported")
})

test_that("supportedFilters works", {
    res <- supportedFilters(cmp_db)
    expect_equal(colnames(res), c("filter", "field"))
    res_2 <- supportedFilters(cmp_spctra_db)
    expect_true(nrow(res) < nrow(res_2))
})

test_that("metadata works", {
    res <- metadata(cdb)
    expect_true(is.data.frame(res))
    expect_true(all(colnames(res) == c("name", "value")))
})

test_that("spectraVariables,CompDb works", {
    db <- new("CompDb")
    expect_equal(spectraVariables(db), character())

    res <- spectraVariables(cdb)
    expect_true(is.character(res))
    expect_true(length(res) > 0)
    expect_true(all(c("spectrum_id", "ms_level") %in% res))
})

test_that("compoundVariables,CompDb works", {
    db <- new("CompDb")
    expect_equal(compoundVariables(db), character())

    res <- compoundVariables(cdb)
    expect_true(is.character(res))
    expect_true(length(res) > 0)
    expect_true(all(c("formula", "inchi") %in% res))

    expect_true(any(compoundVariables(cdb, TRUE) == "compound_id"))
})

test_that("compounds works", {
    res <- compounds(cmp_db, columns = character())
    expect_true(is.data.frame(res))
    expect_true(ncol(res) == 0)
    expect_true(nrow(res) == 0)
    cmps <- compounds(cmp_db)
    expect_true(is(cmps, "data.frame"))
    cmps_tbl <- compounds(cmp_db, columns = c("compound_id", "name"),
                          return.type = "tibble")
    expect_true(is(cmps_tbl, "tbl"))
    expect_equal(colnames(cmps_tbl), c("compound_id", "name"))

    expect_error(compounds(cmp_db, filter = "something"))

    expect_true(
        nrow(compounds(cmp_db, filter = ~ compound_id == "HMDB0000005")) == 1)
    res <- compounds(cmp_spctra_db,
                     columns = c("compound_id", "spectrum_id", "splash"))
    cmp_ids <- compounds(cmp_spctra_db, columns = "compound_id")$compound_id
    expect_true(all(cmp_ids %in% res$compound_id))
    expect_true(sum(is.na(res$spectrum_id)) == 6)

    ## compounds with filters
    res <- compounds(cdb, filter = ~ exactmass > 300)
    expect_true(all(res$exactmass > 300))
    res_2 <- compounds(cdb, filter = ~ exactmass > 300 & exactmass < 340)
    expect_true(nrow(res_2) < nrow(res))

    res <- compounds(cdb, filter = FormulaFilter("C17", "startsWith"))
    expect_true(nrow(res) > 0)
})
