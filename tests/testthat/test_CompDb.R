test_that("CompDb constructor and low level functions", {
    expect_error(CompDb())
    expect_error(CompDb(3))
    cmp <- new("CompDb")
    expect_true(is.null(.dbconn(cmp)))
    expect_true(is.null(dbconn(cmp)))

    ## Create a simple database from the internal files.
    metadata <- data.frame(name = c("source", "url", "source_version",
                                   "source_date", "organism"),
                           value = c("HMDB_tmp", "http://www.hmdb.ca", "4",
                                     "2017", "Hsapiens"),
                           stringsAsFactors = FALSE)
    fl <- system.file("sdf/HMDB_sub.sdf.gz", package = "CompoundDb")
    cmps <- compound_tbl_sdf(fl)
    db_f <- createCompDb(cmps, metadata = metadata, path = tempdir())
    cmp <- CompDb(db_f)

    expect_true(!is.null(.dbconn(cmp)))
    expect_true(.validCompDb(dbconn(cmp)))
    res <- .metadata(cmp)
    expect_equal(metadata, res[1:nrow(metadata), ])
    res <- .metadata(dbconn(cmp))
    expect_equal(metadata, res[1:nrow(metadata), ])
    res <- .metadata_value(cmp, "organism")
    expect_equal(res, "Hsapiens")
    res <- .metadata_value(dbconn(cmp), "source")
    expect_equal(res, "HMDB_tmp")

    ## .tables
    tbls <- .tables(cmp)
    expect_equal(length(tbls), 2)
    expect_equal(names(tbls), c("compound", "synonym"))
    tbls <- .tables(cmp, metadata = TRUE)
    expect_equal(length(tbls), 3)
    expect_equal(names(tbls), c("compound", "metadata", "synonym"))
    tbls <- .tables(cmp, name = "not_there")
    expect_equal(length(tbls), 1)
    tbls <- tables(cmp)
    expect_equal(length(tbls), 2)

    ## tables with spectra
    tbls <- tables(cmp_spctra_db)
    expect_equal(length(tbls), 4)
    expect_equal(names(tbls), c("compound",
                                "msms_spectrum",
                                "msms_spectrum_peak",
                                "synonym"))

    ## .get_property
    prps <- .get_property(cmp, "tables")
    expect_equal(prps, .tables(cmp, metadata = TRUE))
    prps <- .get_property(cmp, "not_there")
    expect_equal(prps, NULL)
})

test_that("compounds works", {
    cmps <- compounds(cmp_db)
    expect_true(is(cmps, "data.frame"))
    cmps_tbl <- compounds(cmp_db, columns = c("compound_id", "compound_name"),
                          return.type = "tibble")
    expect_true(is(cmps_tbl, "tbl"))
    expect_equal(colnames(cmps_tbl), c("compound_id", "compound_name"))

    expect_error(compounds(cmp_db, filter = "something"))

    expect_true(
        nrow(compounds(cmp_db, filter = ~ compound_id == "HMDB0000005")) == 1)
    ## compounds should always return all entries from the compound table, even
    ## if entries from another table are queried.
    res <- compounds(cmp_spctra_db, columns = c("spectrum_id", "splash"))
    cmp_ids <- compounds(cmp_spctra_db, columns = "compound_id")$compound_id
    expect_true(all(cmp_ids %in% res$compound_id))
    expect_true(sum(is.na(res$spectrum_id)) == 6)
})

test_that("src_compound works", {
    src_cmp <- src_compdb(cmp_db)
    expect_true(is(src_cmp, "src_dbi"))
    expect_error(src_compdb(5))
})

test_that(".has_msms_spectra/hasMsMsSpectra works", {
    expect_false(.has_msms_spectra(cmp_db))
    expect_false(hasMsMsSpectra(cmp_db))
    expect_true(.has_msms_spectra(cmp_spctra_db))
    expect_true(hasMsMsSpectra(cmp_spctra_db))
})