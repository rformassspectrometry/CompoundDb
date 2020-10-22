
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

    ## columns
    res <- Spectra(cmp_spctra_db, columns = "compound_id")
    expect_false(all(c("instrument", "predicted") %in% spectraVariables(res)))

    ## filter
    res <- Spectra(cmp_spctra_db, filter = ~ compound_id == "HMDB0000001")
    expect_true(length(res) == 2)

    ## filter and columns
    res <- Spectra(cmp_spctra_db, filter = ~ compound_id == "HMDB0000001",
                   columns = c("inchi", "name"))
    expect_true(all(c("spectrum_id", "name", "inchi") %in%
                    spectraVariables(res)))
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
