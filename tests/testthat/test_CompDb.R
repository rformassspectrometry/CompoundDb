test_that("CompDb works with SQLite database name and connection", {
    a <- CompDb(db_file)
    expect_true(validObject(a))
    expect_identical(a@dbcon, NULL)
    expect_identical(a@dbname, db_file)

    db_con <- dbConnect(SQLite(), db_file)
    b <- CompDb(db_con)
    expect_true(validObject(b))
    expect_identical(b@dbname, character())
    expect_true(length(b@dbcon) > 0)

    expect_equal(compounds(a), compounds(b))
    expect_equal(metadata(a), metadata(b))
    expect_equal(a@.properties, b@.properties)
    dbDisconnect(db_con)
})

test_that(".validCompDb works", {
    tf <- tempfile()
    tmp_db <- dbConnect(SQLite(), tf)
    expect_match(.validCompDb(tmp_db), "found in the database")
    dbWriteTable(tmp_db, name = "metadata", data.frame(a = 1, b = 3))
    dbWriteTable(tmp_db, name = "ms_compound", data.frame(a = 1, b = 2))
    expect_match(.validCompDb(tmp_db), "Miss required columns")
    dbExecute(tmp_db, "drop table metadata")
    dbExecute(tmp_db, "drop table ms_compound")
    .copy_compdb(.dbconn(cmp_spctra_db), tmp_db)
    dbExecute(tmp_db, "drop table msms_spectrum")
    expect_match(.validCompDb(tmp_db), "msms_spectrum and")
    dbDisconnect(tmp_db)
    rm(tf)
    tf <- tempfile()
    tmp_db <- dbConnect(SQLite(), tf)
    .copy_compdb(.dbconn(cmp_spctra_db), tmp_db)
    with_mocked_bindings(
        ".valid_metadata" = function(x, error = FALSE) return("ERROR"),
        code = expect_match(.validCompDb(tmp_db), "ERROR")
    )
    with_mocked_bindings(
        ".valid_msms_spectrum" = function(x, error = FALSE) return("ERROR"),
        code = expect_match(.validCompDb(tmp_db), "ERROR")
    )
    dbExecute(tmp_db, "drop index msms_mid_idx")
    dbExecute(tmp_db, "alter table msms_spectrum drop column 'spectrum_id'")
    expect_match(.validCompDb(tmp_db), "Required column")
    dbDisconnect(tmp_db)
    rm(tf)
    tf <- tempfile()
    tmp_db <- dbConnect(SQLite(), tf)
    .copy_compdb(.dbconn(cmp_spctra_db), tmp_db)
    dbExecute(tmp_db, "alter table msms_spectrum_peak drop column 'mz'")
    expect_match(.validCompDb(tmp_db), "Required columns")
    dbDisconnect(tmp_db)
    rm(tf)
    tf <- tempfile()
    tmp_db <- dbConnect(SQLite(), tf)
    .copy_compdb(.dbconn(cmp_spctra_db), tmp_db)
    dbExecute(tmp_db, "delete from ms_compound where compound_id='HMDB0000001'")
    expect_match(.validCompDb(tmp_db), "Not all compound ids")
    dbDisconnect(tmp_db)
    expect_match(.validCompDb(tmp_db), "not available or closed")
    rm(tf)
})

test_that("CompDb constructor and low level functions", {
    expect_error(CompDb(), "database")
    expect_error(CompDb(3), "database")
    expect_error(CompDb(NA), "database")
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
    expect_equal(names(tbls), c("ms_compound", "synonym"))
    tbls <- .tables(cmp, metadata = TRUE)
    expect_equal(length(tbls), 3)
    expect_equal(names(tbls), c("metadata", "ms_compound", "synonym"))
    tbls <- .tables(cmp, name = "not_there")
    expect_equal(length(tbls), 1)
    tbls <- tables(cmp)
    expect_equal(length(tbls), 2)

    ## tables with spectra
    tbls <- tables(cmp_spctra_db)
    expect_equal(length(tbls), 4)
    expect_equal(names(tbls), c("ms_compound",
                                "msms_spectrum",
                                "msms_spectrum_peak",
                                "synonym"))

    ## .get_property
    prps <- .get_property(cmp, "tables")
    expect_equal(prps, .tables(cmp, metadata = TRUE))
    prps <- .get_property(cmp, "not_there")
    expect_equal(prps, NULL)
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

test_that(".dbconn works", {
    tmp <- new("CompDb")
    expect_identical(.dbconn(tmp), NULL)
    expect_true(is(.dbconn(cmp_db), "SQLiteConnection"))
    tmp <- cmp_db
    tmp@dbcon <- dbConnect(SQLite(), tmp@dbname)
    tmp@dbname <- character()
    expect_true(is(.dbconn(cmp_db), "SQLiteConnection"))
    dbDisconnect(tmp@dbcon)
})

test_that(".dbname works", {
    expect_identical(.dbname(3), character())
    expect_identical(.dbname(cmp_db), db_file)
})

test_that(".require_spectra works", {
    expect_true(.require_spectra())
})

test_that(".initialize_compdb works", {
    tf <- tempfile()
    tmp_db <- dbConnect(SQLite(), tf)
    .copy_compdb(.dbconn(cmp_spctra_db), tmp_db)
    tmp <- .CompDb(dbcon = tmp_db)
    with_mocked_bindings(
        ".validCompDb" = function(x) return("ERROR"),
        code = expect_error(.initialize_compdb(tmp), "ERROR")
    )
    dbDisconnect(tmp_db)
    rm(tf)
})

test_that("copyCompDb works", {
    tf <- tempfile()
    tmp_db <- dbConnect(SQLite(), tf)
    res <- copyCompDb(cmp_spctra_db, tmp_db)
    expect_equal(
        dbGetQuery(.dbconn(cmp_spctra_db), "select * from ms_compound"),
        dbGetQuery(tmp_db, "select * from ms_compound"))
    dbDisconnect(tmp_db)
    rm(tf)
})
