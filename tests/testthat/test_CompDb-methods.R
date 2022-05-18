
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

test_that("insertSpectra,CompDb works", {
    spd <- DataFrame(
        msLevel = c(2L, 2L),
        polarity = c(1L, 1L),
        other_column = "b")
    spd$mz <- list(
        c(109.2, 124.2, 124.5, 170.16, 170.52),
        c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16))
    spd$intensity <- list(
        c(3.407, 47.494, 3.094, 100.0, 13.240),
        c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643))
    sps <- Spectra(spd)
    expect_error(insertSpectra(
        cmp_spctra_db, sps, c("msLevel", "polarity", "other_column")),
                 "Column 'compound_id'")

    sps$compound_id <- c("HMDB0000008", "b")
    expect_error(insertSpectra(
        cmp_spctra_db, sps, c("msLevel", "polarity", "other_column",
                              "compound_id")),
        "variable 'compound_id'")
    sps$compound_id <- c("HMDB0000008", "HMDB0000008")
    expect_error(insertSpectra(
        cmp_spctra_db, sps, c("msLevel", "polarity", "other_column",
                              "compound_id")), "readonly")

    ## Insert to database without spectra data.
    tmp_con <- dbConnect(SQLite(), tempfile())
    CompoundDb:::.copy_compdb(cmp_db@dbcon, tmp_con)

    tmp_db <- CompDb(tmp_con)
    expect_false(CompoundDb:::.has_msms_spectra(tmp_db))
    tmp_db <- insertSpectra(
        tmp_db, sps, c("msLevel", "polarity", "other_column", "compound_id"))
    expect_true(CompoundDb:::.has_msms_spectra(tmp_db))
    res <- dbGetQuery(tmp_con, "select * from msms_spectrum")
    expect_true(all(c("ms_level", "polarity", "other_column", "compound_id")
                    %in% colnames(res)))
    expect_equal(tmp_db@.properties$tables$msms_spectrum, colnames(res))
    expect_true(sum(res$compound_id == "HMDB0000008") == 2)
    expect_true(all(res$other_column[res$compound_id == "HMDB0000008"] == "b"))
    expect_true(length(unique(res$spectrum_id)) == nrow(res))

    res <- dbGetQuery(tmp_con, "select * from msms_spectrum_peak")
    expect_true(sum(res$spectrum_id %in% 1:2) == 12)
    expect_true(length(unique(res$peak_id)) == nrow(res))

    ## Append to existing database.
    tmp_con <- dbConnect(SQLite(), tempfile())
    CompoundDb:::.copy_compdb(cmp_spctra_db@dbcon, tmp_con)

    tmp_db <- CompDb(tmp_con)
    tmp_db <- insertSpectra(
        tmp_db, sps, c("msLevel", "polarity", "other_column", "compound_id"))
    expect_true(CompoundDb:::.has_msms_spectra(tmp_db))
    res <- dbGetQuery(tmp_con, "select * from msms_spectrum")
    expect_true(all(c("ms_level", "polarity", "other_column", "compound_id")
                    %in% colnames(res)))
    expect_equal(tmp_db@.properties$tables$msms_spectrum, colnames(res))
    expect_true(sum(res$compound_id == "HMDB0000008") == 2)
    expect_true(all(res$other_column[res$compound_id == "HMDB0000008"] == "b"))
    expect_true(length(unique(res$spectrum_id)) == nrow(res))

    res <- dbGetQuery(tmp_con, "select * from msms_spectrum_peak")
    expect_true(sum(res$spectrum_id %in% 5:6) == 12)
    expect_true(length(unique(res$peak_id)) == nrow(res))
})

test_that("deleteSpectra,CompDb works", {
    expect_error(deleteSpectra(cmp_spctra_db, ids = c("1", "2")), "readonly")

    tmp_con <- dbConnect(SQLite(), tempfile())
    CompoundDb:::.copy_compdb(cmp_db@dbcon, tmp_con)
    tmp_db <- CompDb(tmp_con)
    expect_false(CompoundDb:::.has_msms_spectra(tmp_db))
    expect_error(deleteSpectra(tmp_db, ids = c("1", "2")), "not contain msms")



    tmp_con <- dbConnect(SQLite(), tempfile())
    CompoundDb:::.copy_compdb(cmp_spctra_db@dbcon, tmp_con)
    tmp_db <- CompDb(tmp_con)
    tmp_db <- deleteSpectra(tmp_db) #should instead the default be delete evrything?
    expect_equal(dbReadTable(dbconn(tmp_db), "msms_spectrum"),
                 dbReadTable(dbconn(cmp_spctra_db), "msms_spectrum"))
    expect_equal(dbReadTable(dbconn(tmp_db), "msms_spectrum_peak"),
                 dbReadTable(dbconn(cmp_spctra_db), "msms_spectrum_peak"))


    tmp_db <- deleteSpectra(tmp_db, ids = c("1", "2"))
    tmp_msms_sp <- dbReadTable(dbconn(cmp_spctra_db), "msms_spectrum")
    exp_msms_sp <- tmp_msms_sp[!tmp_msms_sp$spectrum_id %in% c("1", "2"), ]
    rownames(exp_msms_sp) <- NULL
    expect_equal(dbReadTable(dbconn(tmp_db), "msms_spectrum"), exp_msms_sp)
    tmp_msms_p <- dbReadTable(dbconn(cmp_spctra_db), "msms_spectrum_peak")
    exp_msms_p <- tmp_msms_p[!tmp_msms_p$spectrum_id %in% c("1", "2"), ]
    rownames(exp_msms_p) <- NULL
    expect_equal(dbReadTable(dbconn(tmp_db), "msms_spectrum_peak"), exp_msms_p)
})

test_that("mass2mz,CompDb works",{
    ads <- c("[M+H]+", "[M+Na]+", "[M+K]+")

    #Default adduct as [M+H]+
    expect_identical(mass2mz(cmp_db), mass2mz(cmp_db, "[M+H]+"))

    output <- mass2mz(cmp_db, ads, "compound_id")
    expect_equal(nrow(output), nrow(compounds(cmp_db, "compound_id")))
    expect_equal(ncol(output), length(ads))
    expect_equal(rownames(output), compounds(cmp_db, "compound_id")$compound_id)

    output <- mass2mz(cmp_db, ads, "formula")
    expect_equal(nrow(output), nrow(compounds(cmp_db, "formula")))
    expect_equal(ncol(output), length(ads))
    expect_equal(rownames(output), compounds(cmp_db, "formula")$formula)
})

test_that("mass2mz,ANY works", {
    cmps <- compounds(cdb, c("formula", "exactmass"))
    res <- mass2mz(cmps$exactmass, adduct = c("[M+H]+", "[M+Na]+"))
    res_2 <- mass2mz(cdb, adduct = c("[M+H]+", "[M+Na]+"))
    rownames(res_2) <- NULL
    expect_equal(res, res_2)
})

test_that("insertCompound,CompDb works", {
    db <- emptyCompDb(tempfile())
    res <- insertCompound(db, compounds = data.frame())
    expect_equal(compounds(db), compounds(res))

    cmp <- data.frame(compound_id = 1:3, name = c("a", "b", "c"))
    res <- insertCompound(db, compounds = cmp)
    res_c <- compounds(res)
    expect_equal(colnames(res_c), compoundVariables(res))
    expect_equal(res_c$name, cmp$name)

    ## additional columns.
    cmp$add_col <- 5
    res <- insertCompound(db, compounds = cmp)
    res_c <- compounds(res)
    expect_equal(colnames(res_c), compoundVariables(res))
    expect_equal(res_c$name, cmp$name)
    library(RSQLite)
    all <- dbGetQuery(dbconn(res), "select * from ms_compound")
    expect_equal(all$name, c(cmp$name, cmp$name))
    expect_true(!any(colnames(all) == "add_col"))

    res <- insertCompound(res, compounds = cmp, addColumns = TRUE)
    expect_true(any(compoundVariables(res) == "add_col"))
    all <- dbGetQuery(dbconn(res), "select * from ms_compound")
    expect_equal(all$name, c(cmp$name, cmp$name, cmp$name))
    expect_true(any(colnames(all) == "add_col"))
    expect_true(all(all$add_col[7:9] == 5))

    ## synonyms.
    cmp <- data.frame(compound_id = c("8", "9"), name = c("first", "second"),
                      synonyms = c("primo", NA))
    res <- insertCompound(res, compounds = cmp)
    syns <- dbGetQuery(dbconn(res), "select * from synonym")
    expect_equal(syns$compound_id, "8")
    expect_equal(syns$synonym, "primo")

    cmp$synonyms <- list(c(), c("secondo", "segundo", "zweiter"))
    res <- insertCompound(res, compounds = cmp)
    syns <- dbGetQuery(dbconn(res), "select * from synonym")
    expect_equal(syns$compound_id, c("8", "9", "9", "9"))
    expect_equal(syns$synonym, c("primo", "secondo", "segundo", "zweiter"))

    ## errors
    expect_error(insertCompound(db, compounds = "d"), "data.frame")
    expect_error(insertCompound(new("CompDb"), cmp), "not initialized")
})
