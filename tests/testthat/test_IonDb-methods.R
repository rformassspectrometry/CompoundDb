test_that("show,IonDb works", {
    expect_output(show(ion_db))
    db <- new("IonDb")
    expect_output(show(ion_db))
})

test_that("IonDb works", {
    ## empty.
    res <- IonDb()
    expect_true(is(res, "IonDb"))

    ## errors
    tf <- tempfile()
    tc <- dbConnect(SQLite(), tf)
    expect_error(IonDb(tc), "not found in")
    dbDisconnect(tc)
    rm(tf)

    tc <- dbConnect(SQLite(), db_file)
    expect_error(IonDb(tc), "not found in")
    dbDisconnect(tc)

    ## from existing CompDb using character
    dbf <- tempfile()
    res <- IonDb(dbf, cmp_db)
    expect_true(is(res, "IonDb"))
    expect_true(nrow(ions(res)) == 0)
    expect_equal(compounds(cmp_db), compounds(res))
    rm(dbf)

    ## with ions
    ins <- data.frame(compound_id = "HMDB0000002", ion_adduct = "some",
                      ion_mz = 13.2, ion_rt = 123)
    dbf <- tempfile()
    res <- IonDb(dbf, cmp_db, ions = ins)
    expect_true(is(res, "IonDb"))
    expect_true(nrow(ions(res)) == 1)
    expect_equal(ions(res), ins)
    expect_equal(compounds(cmp_db), compounds(res))

    ## from existing CompDb using dbconnection
    rm(dbf)
    dbf <- tempfile()
    con <- dbConnect(dbDriver("SQLite"), dbf)
    res <- IonDb(con, cmp_spctra_db)
    expect_true(is(res, "IonDb"))
    expect_true(nrow(ions(res)) == 0)
    expect_equal(compounds(cmp_db), compounds(res))
    rm(dbf)

    ## with ions

    ## from existing (ro) CompDb
    expect_error(res <- IonDb(cmp_db), "readonly")

    ## from existing (rw) CompDb
    dbf <- tempfile()
    con <- dbConnect(dbDriver("SQLite"), dbf)
    .copy_compdb(dbconn(cmp_db), con)
    cdb_dbf <- CompDb(con)
    res <- IonDb(cdb_dbf)
    expect_true(is(res, "IonDb"))
    expect_true(nrow(ions(res)) == 0)
    expect_equal(compounds(cmp_db), compounds(res))

    ## load previous IonDb from character
    res <- IonDb(dbf)
    expect_true(is(res, "IonDb"))
    expect_true(nrow(ions(res)) == 0)
    expect_equal(compounds(cmp_db), compounds(res))

    ## load previous IonDb from dbconnection
    res <- IonDb(con)
    expect_true(is(res, "IonDb"))
    expect_true(nrow(ions(res)) == 0)
    expect_equal(compounds(cmp_db), compounds(res))
    dbDisconnect(con)
    rm(dbf)
})

# test_that("Spectra,IonDb works", {
#
#     res <- Spectra(ion_spctra_db)
#     expect_true(is(res, "Spectra"))
#     expect_true(length(res) == 4)
#     expect_true(all(c("instrument", "predicted") %in% spectraVariables(res)))
#
#     ## columns
#     res <- Spectra(ion_spctra_db, columns = "compound_id")
#     expect_false(all(c("instrument", "predicted") %in% spectraVariables(res)))
#
#     ## filter
#     res <- Spectra(ion_spctra_db, filter = ~ compound_id == "HMDB0000001")
#     expect_true(length(res) == 2)
#
#     ## filter and columns
#     res <- Spectra(ion_spctra_db, filter = ~ compound_id == "HMDB0000001",
#                    columns = c("inchi", "name"))
#     expect_true(all(c("spectrum_id", "name", "inchi") %in%
#                     spectraVariables(res)))
#     expect_true(length(res) == 2)
#
#     expect_error(Spectra(ion_spctra_db, filter = "ad"), "'filter' has to")
#     expect_error(Spectra(ion_spctra_db, filter = ~ gene_name == "b"),
#                  "not supported")
# })
#
# test_that("supportedFilters works", {
#     res <- supportedFilters(ion_db)
#     expect_equal(colnames(res), c("filter", "field"))
#     res_2 <- supportedFilters(ion_spctra_db)
#     expect_true(nrow(res) < nrow(res_2))
# })
#
# test_that("metadata works", {
#     res <- metadata(ion_db)
#     expect_true(is.data.frame(res))
#     expect_true(all(colnames(res) == c("name", "value")))
# })
#
# test_that("spectraVariables,IonDb works", {
#     db <- new("IonDb")
#     expect_equal(spectraVariables(db), character())
#
#     res <- spectraVariables(ion_spctra_db)
#     expect_true(is.character(res))
#     expect_true(length(res) > 0)
# })
#
# test_that("compoundVariables,IonDb works", {
#     db <- new("CompDb")
#     expect_equal(compoundVariables(db), character())
#
#     res <- compoundVariables(ion_db)
#     expect_true(is.character(res))
#     expect_true(length(res) > 0)
#     expect_true(all(c("formula", "inchi") %in% res))
#
#     expect_true(any(compoundVariables(ion_db, TRUE) == "compound_id"))
# })
#
# test_that("compounds works", {
#     res <- compounds(ion_db, columns = character())
#     expect_true(is.data.frame(res))
#     expect_true(ncol(res) == 0)
#     expect_true(nrow(res) == 0)
#     cmps <- compounds(ion_db)
#     expect_true(is(cmps, "data.frame"))
#     cmps_tbl <- compounds(ion_db, columns = c("compound_id", "name"),
#                           return.type = "tibble")
#     expect_true(is(cmps_tbl, "tbl"))
#     expect_equal(colnames(cmps_tbl), c("compound_id", "name"))
#
#     expect_error(compounds(ion_db, filter = "something"))
#
#     expect_true(
#         nrow(compounds(ion_db, filter = ~ compound_id == "HMDB0000005")) == 1)
#     res <- compounds(cmp_spctra_db,
#                      columns = c("compound_id", "spectrum_id", "splash"))
#     cmp_ids <- compounds(cmp_spctra_db, columns = "compound_id")$compound_id
#     expect_true(all(cmp_ids %in% res$compound_id))
#     expect_true(sum(is.na(res$spectrum_id)) == 6)
#
#     ## compounds with filters
#     res <- compounds(cdb, filter = ~ exactmass > 300)
#     expect_true(all(res$exactmass > 300))
#     res_2 <- compounds(cdb, filter = ~ exactmass > 300 & exactmass < 340)
#     expect_true(nrow(res_2) < nrow(res))
#
#     res <- compounds(cdb, filter = FormulaFilter("C17", "startsWith"))
#     expect_true(nrow(res) > 0)
# })


test_that("ionVariables,ionDb works", {
    db <- new("IonDb")
    expect_equal(ionVariables(db), character())

    res <- ionVariables(ion_db)
    expect_true(is.character(res))
    expect_true(length(res) > 0)
    expect_true(all(c("compound_id", "ion_adduct", "ion_mz", "ion_rt")
                    %in% res))
    expect_true(any(ionVariables(ion_db, TRUE) == "ion_id"))
})

test_that("ions,ionDb works", {
    res <- ions(ion_db, columns = character())
    expect_true(is.data.frame(res))
    expect_true(ncol(res) == 0)
    expect_true(nrow(res) == 0)
    expect_true(is(ions(ion_db), "data.frame"))
    ions_tbl <- ions(ion_db, columns = c("ion_id", "compound_id", "ion_adduct"),
                          return.type = "tibble")
    expect_true(is(ions_tbl, "tbl"))
    expect_equal(colnames(ions_tbl), c("ion_id", "compound_id", "ion_adduct"))

    expect_error(compounds(ion_db, filter = "something"))

    ## ions with filters
    res <- ions(ion_db, filter = ~ ion_rt < 100)
    expect_true(all(res$ion_rt < 100))
    res_2 <- ions(ion_db, filter = ~ ion_rt < 100 & ion_mz < 200)
    expect_true(nrow(res_2) < nrow(res))

    res <- ions(ion_db, filter = CompoundIdFilter("01", "endsWith"))
    expect_true(nrow(res) > 0)
})

test_that("insertIon,ionDb works", {
    dbf <- tempfile()
    con <- dbConnect(dbDriver("SQLite"), dbf)
    .copy_compdb(dbconn(ion_spctra_db), con)
    idb <- IonDb(con)

    more_ions <- data.frame(compound_id = c("HMDB0000005", "HMDB0000008"),
                            ion_adduct = c("C", "D"),
                            ion_mz = c(220, 300),
                            ion_rt = c(90, 140),
                            ion_id = c(1, 2))
    expect_warning(idb <- insertIon(idb, more_ions), "replaced with")
    expect_true(nrow(ions(idb)) == 7)
    expect_equal(ions(idb), rbind(ions(ion_spctra_db), more_ions[, -5L]))

    ## Different ordering of columns
    more_ions <- more_ions[, c(3, 1, 2, 4)]
    idb <- insertIon(idb, more_ions)
    expect_true(nrow(ions(idb, c("ion_id"))) == 9)
    res <- ions(idb, c("ion_id", "compound_id", "ion_adduct"))
    expect_equal(res$ion_adduct, c("A", "B", "B", "C", "D", "C", "D", "C", "D"))

    ## Additional columns
    more_ions$add_column <- "Z"
    expect_error(insertIon(idb, more_ions), "no column")

    idb <- insertIon(idb, more_ions, addColumns = TRUE)
    expect_true(any(ionVariables(idb) == "add_column"))
    res <- ions(idb)
    expect_identical(res$add_column, c(rep(NA_character_, 7), "Z", "Z"))
    expect_true(any(idb@.properties$tables$ms_ion == "add_column"))

    ## Errors
    expect_error(insertIon(idb, more_ions[, 1:3]), "required")
    more_ions$compound_id <- c("a", "b")
    expect_error(insertIon(idb, more_ions), "compound_id")
    ## dbDisconnect(con)

    tmp <- IonDb()
    expect_error(insertIon(tmp, more_ions, addColumns = TRUE),
                 "not initialized")
})

test_that("insertSpectra,IonDb works", {
    spd <- DataFrame(
        msLevel = c(2L, 2L, 2L),
        polarity = c(1L, 1L, 1L),
        collisionEnergy = c(20, 30, 30))
    spd$mz <- list(
        c(109.2, 124.2, 124.5, 170.16, 170.52),
        c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16), c(50, 100))
    spd$intensity <- list(
        c(3.407, 47.494, 3.094, 100.0, 13.240),
        c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643), c(50, 100))
    sps <- Spectra(spd)
    dbf <- tempfile()
    con <- dbConnect(dbDriver("SQLite"), dbf)
    .copy_compdb(dbconn(ion_spctra_db), con)
    idb <- IonDb(con)
    msms_sp <- dbReadTable(dbconn(idb), "msms_spectrum")
    msms_sp_peak <- dbReadTable(dbconn(idb), "msms_spectrum_peak")
    expect_error(insertSpectra(idb, sps), "'compound_id' needs")
    sps$compound_id <- c("HMDB0000001", "HMDB0000001", "HMDB1000000")
    expect_error(insertSpectra(idb, sps), paste("All values of spectra variable",
                                                "'compound_id' must be in"))
    sps <- sps[1:2]
    insertSpectra(idb, sps, columns = c("msLevel", "polarity",
                                        "collisionEnergy", "compound_id"))
    msms_sp2 <- dbReadTable(dbconn(idb), "msms_spectrum")
    msms_sp_peak2 <- dbReadTable(dbconn(idb), "msms_spectrum_peak")
    ns <- nrow(msms_sp)
    ns2 <- nrow(msms_sp2)
    expect_equal(ns2, ns + length(sps))
    expect_equal(msms_sp, msms_sp2[1:ns, colnames(msms_sp)])
    np <- nrow(msms_sp_peak)
    np2 <- nrow(msms_sp_peak2)
    expect_equal(np2, np + sum(lengths(peaksData(sps))))
    expect_equal(msms_sp_peak2$peak_id, 1:np2)
    expect_equal(msms_sp_peak2[, c("mz", "intensity")],
                 rbind(msms_sp_peak[, c("mz", "intensity")],
                       do.call(rbind, Spectra:::.peaksapply(sps))))
    expect_equal(msms_sp_peak2$spectrum_id,
                 c(msms_sp_peak$spectrum_id,
                   ns + rep(1:length(sps), lengths(peaksData(sps)))))

    # Add additional column in spectra object to the database object
    sps$newvariable <- c("value1", "value2")
    dbf <- tempfile()
    con <- dbConnect(dbDriver("SQLite"), dbf)
    .copy_compdb(dbconn(ion_spctra_db), con)
    idb <- IonDb(con)
    msms_sp <- dbReadTable(dbconn(idb), "msms_spectrum")
    expect_error(insertSpectra(idb, sps, columns = c("newvariable2",
                                                     "compound_id")),
                 "not available")
    insertSpectra(idb, sps, columns = c("newvariable", "compound_id"))
    msms_sp2 <- dbReadTable(dbconn(idb), "msms_spectrum")
    ns <- nrow(msms_sp)
    ns2 <- nrow(msms_sp2)
    expect_equal(ns2, ns + length(sps))
    expect_true("newvariable" %in% colnames(msms_sp2))
    expect_equal(c(rep(NA, ns), sps$newvariable), msms_sp2$newvariable)
})

test_that("deleteIon,IonDb works", {
    dbf <- tempfile()
    con <- dbConnect(dbDriver("SQLite"), dbf)
    .copy_compdb(dbconn(ion_spctra_db), con)
    idb <- IonDb(con)
    ids <- c(1, 2)
    deleteIon(idb, ids)
    tmp <- ions(ion_spctra_db, ionVariables(ion_spctra_db, includeId = TRUE))
    expected_ions <- tmp[-match(ids, tmp$ion_id), ]
    rownames(expected_ions) <- NULL
    expect_equal(ions(idb, ionVariables(idb, includeId = TRUE)), expected_ions)
    dbDisconnect(con)
    rm(dbf)

    dbf <- tempfile()
    con <- dbConnect(dbDriver("SQLite"), dbf)
    CompoundDb:::.copy_compdb(dbconn(ion_spctra_db), con)
    idb <- IonDb(dbf)
    ids <- c(1, 2, 30)
    expect_warning(deleteIon(idb, ids), "Some IDs")
    tmp <- ions(ion_spctra_db, ionVariables(ion_spctra_db, includeId = TRUE))
    expected_ions <- tmp[-match(ids[-3], tmp$ion_id), ]
    rownames(expected_ions) <- NULL
    expect_equal(ions(idb, ionVariables(idb, includeId = TRUE)), expected_ions)
    rm(fbf)

    expect_error(deleteIon(IonDb(), ids = 1:3), "not initialized")
})

test_that("deleteCompound,IonDb works", {
    db <- IonDb(emptyCompDb(tempfile()))
    res <- deleteCompound(db, 1:4)
    expect_equal(compounds(res), compounds(db))
    expect_equal(ions(res), ions(db))

    res <- deleteCompound(db)
    expect_equal(compounds(res), compounds(db))

    cmp <- data.frame(compound_id = 1:4, name = letters[1:4], synonyms = 1:4)
    db <- insertCompound(db, cmp)
    ins <- data.frame(compound_id = c("2", "3"), ion_adduct = "[M+H]+",
                      ion_mz = c(1.2, 4.3), ion_rt = c(4.3, 4.5))
    db <- insertIon(db, ins)

    library(RSQLite)
    db <- deleteCompound(db, ids = 1)
    expect_true(!any(compounds(db, "compound_id")$compound_id == 1))

    expect_error(deleteCompound(db, ids = c(3, 6, 8, 10)), "compounds present")

    ## With Spectra data
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
    sps$compound_id <- as.character(c(2, 4))

    db <- insertSpectra(db, sps)

    expect_error(deleteCompound(db, ids = c(1, 2)), "Ions")
    expect_error(deleteCompound(db, ids = c(1, 4)), "MS2")

    db <- deleteCompound(db, ids = c(1, 2), recursive = TRUE)
    expect_true(!any(ions(db)$compound_id %in% c(1, 2)))
    expect_equal(ions(db)$compound_id, "3")

    expect_true(!any(compounds(db, "compound_id")$compound_id %in% c(1, 2)))
    expect_equal(
        dbGetQuery(dbconn(db), "select * from msms_spectrum")$compound_id, "4")

    ## errors
    expect_error(deleteCompound(new("IonDb"), 2), "not initialized")
})
