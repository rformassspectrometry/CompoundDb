test_that(".create_ion_table works", {
    hmdb <- system.file("sdf/HMDB_sub.sdf.gz", package = "CompoundDb")
    cmps <- compound_tbl_sdf(hmdb)
    metad <- data.frame(name = c("source", "url", "source_version",
                                 "source_date", "organism"),
                        value = c("HMDB_test", "http://www.hmdb.ca",
                                  "v4", "2017-08-27", "Hsapiens"))
    dr <- paste0(tempdir(), "/rw")
    dir.create(dr, showWarnings = FALSE)
    db_file <- createCompDb(cmps, metadata = metad, path = dr)
    conn <- dbConnect(SQLite(), db_file)
    .create_ion_table(conn)
    tbls <- dbListTables(conn)
    expect_true(all(c("ms_compound", "ms_ion", "metadata", "synonym")
                    %in% tbls))
    tbl <- dbGetQuery(conn, "select * from ms_ion")
    expect_true(
        all(c("ion_id", "compound_id", "ion_adduct", "ion_mz", "ion_rt") %in%
            colnames(tbl)))

    dbExecute(conn, "drop table ms_ion;")

    ions <- data.frame(ion_id = character(), compound_id = character(),
                       mz = numeric(), rt = numeric(), adduct = character(),
                       other_col = integer())
    .create_ion_table(conn, ions = ions)
    tbl <- dbGetQuery(conn, "select * from ms_ion")
    expect_true(
        all(c("ion_id", "compound_id", "adduct", "mz", "rt", "other_col") %in%
            colnames(tbl)))
})

test_that(".valid_ion works", {
    ions <- data.frame(compound_id = "a", ion_adduct = "f", ion_mz = 4.2,
                       ion_rt = 4.2)
    expect_true(.valid_ion(ions))
    ions$compound_id <- 4
    expect_error(.valid_ion(ions), "should be character")
    ions$compound_id <- "a"
    ions$ion_mz <- NULL
    expect_error(.valid_ion(ions), "ion_mz")
    ions$ion_mz <- 45.3
    ions$add_col <- "d"
    expect_true(.valid_ion(ions))
    ions$compound_id <- NA
    expect_error(.valid_ion(ions), "missing")
})

test_that(".copy_compdb works", {
    a <- dbConnect(SQLite(), tempfile())
    .copy_compdb(dbconn(cmp_db), a)
    expect_equal(dbListTables(a), dbListTables(dbconn(cmp_db)))

    idx <- "SELECT name FROM sqlite_master WHERE type = 'index';"
    expect_equal(sort(dbGetQuery(dbconn(cmp_db), idx)$name),
                 sort(dbGetQuery(a, idx)$name))

    a <- dbConnect(SQLite(), tempfile())
    .copy_compdb(dbconn(cmp_spctra_db), a)
    expect_equal(dbListTables(a), dbListTables(dbconn(cmp_spctra_db)))
    expect_equal(sort(dbGetQuery(dbconn(cmp_spctra_db), idx)$name),
                 sort(dbGetQuery(a, idx)$name))
})
