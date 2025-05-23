test_that(".simple_extract_compounds_sdf works", {
    hmdb <- system.file("sdf/HMDB_sub.sdf.gz", package = "CompoundDb")

    with_mocked_bindings(
        ".guess_sdf_source" = function(x) return(NULL),
        code = expect_error(.simple_extract_compounds_sdf(hmdb), "The SDF file")
    )
    cmps <- .simple_extract_compounds_sdf(
        datablock2ma(datablock(read.SDFset(hmdb))))
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "name", "inchi",
                                   "inchikey", "formula", "exactmass",
                                   "synonyms", "smiles"))
    expect_true(nrow(cmps) == 9)

    chebi <- system.file("sdf/ChEBI_sub.sdf.gz", package = "CompoundDb")
    cmps <- .simple_extract_compounds_sdf(
        datablock2ma(datablock(read.SDFset(chebi))))
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "name", "inchi",
                                   "inchikey", "formula", "exactmass",
                                   "synonyms", "smiles"))
    expect_true(nrow(cmps) == 6)

    lm <- system.file("sdf/LipidMaps_sub.sdf.gz", package = "CompoundDb")
    cmps <- CompoundDb:::.simple_extract_compounds_sdf(
        datablock2ma(datablock(read.SDFset(lm))))
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "name", "inchi",
                                   "inchikey", "formula", "exactmass",
                                   "synonyms", "smiles"))
    expect_true(nrow(cmps) == 7)

    lm <- system.file("sdf/LipidMaps_mock.sdf", package = "CompoundDb")
    cmps2 <- CompoundDb:::.simple_extract_compounds_sdf(
                              datablock2ma(datablock(read.SDFset(lm))))
    expect_true(any(cmps2$name == "syns"))
    expect_true(any(cmps2$name == "sysname"))

    pubchem <- system.file("sdf/PubChem_sub.sdf.gz", package = "CompoundDb")
    cmps <- .simple_extract_compounds_sdf(
        datablock2ma(datablock(read.SDFset(pubchem))))
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "name", "inchi",
                                   "inchikey", "formula", "exactmass",
                                   "synonyms", "smiles"))
    expect_true(nrow(cmps) == 12)

    mona <- system.file("sdf/MoNa_export-All_Spectra_sub.sdf.gz",
                        package = "CompoundDb")
    cmps <- .simple_extract_compounds_sdf(
        datablock2ma(datablock(read.SDFset(mona))))
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "name", "inchi",
                                   "inchikey", "formula", "exactmass",
                                   "synonyms", "smiles"))
    expect_true(nrow(cmps) == 7)
})

test_that("compound_tbl_sdf works", {
    expect_error(compound_tbl_sdf())
    expect_error(compound_tbl_sdf("somefile"))

    hmdb <- system.file("sdf/HMDB_sub.sdf.gz", package = "CompoundDb")
    cmps <- compound_tbl_sdf(hmdb)
    cmps2 <- compound_tbl_sdf(hmdb, onlyValid = FALSE)
    expect_equal(cmps, cmps2)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "name", "inchi",
                                   "inchikey", "formula", "exactmass",
                                   "synonyms", "smiles"))
    expect_true(nrow(cmps) == 9)
    expect_true(is(cmps$synonyms, "list"))
    cmps <- compound_tbl_sdf(hmdb, collapse = "|")
    expect_true(is(cmps$synonyms, "character"))

    with_mocked_bindings(
        "validSDF" = function(x) c(1, 3),
        code = expect_message(compound_tbl_sdf(hmdb), "Skipped import")
    )

    chebi <- system.file("sdf/ChEBI_sub.sdf.gz", package = "CompoundDb")
    cmps <- compound_tbl_sdf(chebi)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "name", "inchi",
                                   "inchikey", "formula", "exactmass",
                                   "synonyms", "smiles"))
    expect_true(nrow(cmps) == 6)
    expect_true(is(cmps$synonyms, "list"))
    cmps <- compound_tbl_sdf(chebi, collapse = "|")
    expect_true(is(cmps$synonyms, "character"))

    lm <- system.file("sdf/LipidMaps_sub.sdf.gz", package = "CompoundDb")
    cmps <- compound_tbl_sdf(lm)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "name", "inchi",
                                   "inchikey", "formula", "exactmass",
                                   "synonyms", "smiles"))
    expect_true(nrow(cmps) == 7)
    expect_true(is(cmps$synonyms, "list"))
    cmps <- compound_tbl_sdf(lm, collapse = "|")
    expect_true(is(cmps$synonyms, "character"))

    pc <- system.file("sdf/PubChem_sub.sdf.gz", package = "CompoundDb")
    cmps <- compound_tbl_sdf(pc)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "name", "inchi",
                                   "inchikey", "formula", "exactmass",
                                   "synonyms", "smiles"))
    expect_true(nrow(cmps) == 12)
    expect_true(is(cmps$synonyms, "list"))
    cmps <- compound_tbl_sdf(pc, collapse = "|")
    expect_true(is(cmps$synonyms, "character"))
})

test_that(".guess_sdf_source works", {
    hmdb <- system.file("sdf/HMDB_sub.sdf.gz", package = "CompoundDb")
    library(ChemmineR)
    cn <- colnames(datablock2ma(datablock(read.SDFset(hmdb))))
    expect_true(.guess_sdf_source(cn) == "hmdb")

    chebi <- system.file("sdf/ChEBI_sub.sdf.gz", package = "CompoundDb")
    cn <- colnames(datablock2ma(datablock(read.SDFset(chebi))))
    expect_true(.guess_sdf_source(cn) == "chebi")

    lm <- system.file("sdf/LipidMaps_sub.sdf.gz", package = "CompoundDb")
    cn <- colnames(datablock2ma(datablock(read.SDFset(lm))))
    expect_true(.guess_sdf_source(cn) == "lipidmaps")

    pc <- system.file("sdf/PubChem_sub.sdf.gz", package = "CompoundDb")
    cn <- colnames(datablock2ma(datablock(read.SDFset(pc))))
    expect_true(.guess_sdf_source(cn) == "pubchem")

    expect_equal(.guess_sdf_source(c("some", "thing")), NULL)
})

test_that(".import_lipidblast works", {
    fl <- system.file("json/MoNa-LipidBlast_sub.json", package = "CompoundDb")
    cmps <- .import_lipidblast(fl)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_true(all(c("compound_id", "name", "inchi", "inchikey",
                      "formula", "exactmass", "synonyms") %in% colnames(cmps)))
    expect_true(nrow(cmps) == 8)
})

test_that("compound_tbl_lipidblast works", {
    expect_error(compound_tbl_lipidblast())
    expect_error(compound_tbl_lipidblast("sddfd"))

    lb <- system.file("json/MoNa-LipidBlast_sub.json", package = "CompoundDb")
    cmps <- compound_tbl_lipidblast(lb)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_true(all(c("compound_id", "name", "inchi", "inchikey", "formula",
                       "exactmass", "synonyms") %in% colnames(cmps)))
    expect_true(nrow(cmps) == 8)
    expect_true(is.list(cmps$synonyms))
    cmps <- compound_tbl_lipidblast(lb, collapse = ";")
    expect_true(is.character(cmps$synonyms))
})

test_that(".valid_metadata works", {
    ## Check errors.
    expect_error(.valid_metadata("something"))
    expect_true(is.character(
        .valid_metadata("something", error = FALSE)))
    metadata <- data.frame(name = c("source", "url", "source_version",
                                    "source_date", "organism"),
                           value = c("HMDB1", "http://www.hmdb.ca", "4", "2017",
                                     "Hsapiens"),
                           stringsAsFactors = FALSE)
    expect_error(.valid_metadata(metadata[, 1, drop = FALSE]))
    expect_error(.valid_metadata(metadata[1:2, ]))

    ## Valid one.
    expect_true(.valid_metadata(metadata))
})

test_that(".db_file_from_metadata works", {
    metadata <- data.frame(name = c("source", "url", "source_version",
                                    "source_date", "organism"),
                           value = c("HMDB", "http://www.hmdb.ca", "v4", "2017",
                                     "Hsapiens"),
                           stringsAsFactors = FALSE)
    db_file <- .db_file_from_metadata(metadata)
    expect_equal(db_file, "CompDb.Hsapiens.HMDB.v4")
})

test_that(".valid_compound works", {
    cmps <- data.frame(compound_id = c("01", "02"), name = c("a", "b"),
                       inchi = c("i1", "i2"), inchikey = c("k1", "k2"),
                       formula = c("some", "thing"),
                       exactmass = c(1, 3), synonyms = c("a", "b"))
    expect_true(.valid_compound(cmps, db = FALSE))
    expect_true(.valid_compound(cmps[, 1:6]))
    expect_error(.valid_compound(cmps[, 1:6], db = FALSE))
    expect_true(.valid_compound(cmps, db = TRUE))

    ## Errors
    expect_error(.valid_compound("b"))
    expect_true(is.character(.valid_compound("b", error = FALSE)))
    expect_error(.valid_compound(data.frame()))
    expect_error(.valid_compound(cmps[, 1:3]))
    cmps$exactmass <- c("1", "2")
    expect_error(.valid_compound(cmps))

    cmps$compound_id <- ""
    cmps$exactmass <- 1.3
    expect_error(.valid_compound(cmps), "not allowed")
    cmps$compound_id <- NA_character_
    expect_error(.valid_compound(cmps), "not allowed")
})

test_that("createCompDb and createCompDbPackage works", {
    fl <- system.file("sdf/ChEBI_sub.sdf.gz", package = "CompoundDb")
    cmps <- compound_tbl_sdf(fl)

    metad <- data.frame(name = c("source", "url", "source_version",
                                 "source_date", "organism"),
                        value = c("ChEBI", "http://www.hmdb.ca",
                                  "unknown", "2017-08-27", "Hsapiens"),
                        stringsAsFactors = FALSE)
    db_f <- createCompDb(cmps, metadata = metad, path = tempdir())
    expect_true(is(db_f, "character"))

    cdb <- CompDb(db_f)

    ## createCompDbPackage:
    res <- createCompDbPackage(db_f, version = "0.0.1",
                               maintainer = "John Doe <john.doe@mail.com>",
                               author = "J Doe", path = tempdir())
    expect_true(is.character(res))
    expect_equal(basename(res), "CompDb.Hsapiens.ChEBI.unknown")
    expect_error(createCompDbPackage(5, version = "0.0.1",
                                     maintainer = "John Doe <john.doe@mail.com>",
                                     author = "J Doe", path = tempdir()))
    expect_error(createCompDbPackage(dbf, version = "0.0.1",
                                     author = "J Doe", path = tempdir()))


    ## Provide a single file name.
    fl <- system.file("sdf/LipidMaps_sub.sdf.gz", package = "CompoundDb")
    metad <- make_metadata(source = "LipidMaps", source_date = "2016",
                           source_version = "xx", organism = "Hsapiens",
                           url = NA)
    res <- createCompDb(fl, metadata = metad, path = tempdir())
    db <- CompDb(res)
    md <- .metadata(db)
    expect_true(is.na(md$value[md$name == "url"]))
    ## Multiple files.
    fls <- c(system.file("sdf/LipidMaps_sub.sdf.gz", package = "CompoundDb"),
             system.file("sdf/HMDB_sub.sdf.gz", package = "CompoundDb"),
             system.file("sdf/ChEBI_sub.sdf.gz", package = "CompoundDb")
             )
    metad <- make_metadata(source = "Multiple", source_date = "2016",
                           source_version = "xx", organism = "Hsapiens",
                           url = NA)
    res <- createCompDb(fls, metadata = metad, path = tempdir())
    db <- CompDb(res)
    expect_true(nrow(compounds(db)) == 22)

    ## Multiple files including json.
    fls <- c(fls, system.file("json/MoNa-LipidBlast_sub.json",
                              package = "CompoundDb"))
    metad <- make_metadata(source = "EvenMore", source_date = "2016",
                           source_version = "xx", organism = "Hsapiens",
                           url = NA)
    res <- createCompDb(fls, metadata = metad, path = tempdir())
    db <- CompDb(res)
    expect_true(nrow(compounds(db)) == 30)

    ## Error with one unsupported file.
    fls <- c(fls, system.file("NEWS", package = "CompoundDb"))
    metad <- make_metadata(source = "Fails", source_date = "2016",
                           source_version = "xx", organism = "Hsapiens",
                           url = NA)
    expect_error(createCompDb(fls, metadata = metad, path = tempdir()))

    a <- system.file("sdf/HMDB_sub.sdf.gz", package = "CompoundDb")
    cmps <- compound_tbl_sdf(a)
    cmps <- cmps[-1L, ]
    metad <- data.frame(name = c("source", "url", "source_version",
                                 "source_date", "organism"),
                        value = c("HMDB_testx", "http://www.hmdb.ca",
                                  "v4", "2017-08-27", "Hsapiens"))
    dr <- system.file("xml/", package = "CompoundDb")
    msms_spctra <- msms_spectra_hmdb(dr)
    expect_error(createCompDb(cmps, metadata = metad, path = tempdir(),
                              msms_spectra = msms_spctra),
                 "All compound")
})

test_that(".is_sdf_filename works", {
    expect_true(.is_sdf_filename("somesdf.file.sdf"))
    expect_true(.is_sdf_filename("somesdf.file.SDF"))
    expect_true(.is_sdf_filename("somesdf.file.SDF.gz"))
    expect_true(.is_sdf_filename("somesdf.file.SDF.GZ"))
    expect_false(.is_sdf_filename("somesdf.file.SDF.GZIP"))
    expect_false(.is_sdf_filename("somesdf.file.SDF.zup"))

    expect_true(all(.is_sdf_filename(c("bla.sdf", "blu.sdf.gz"))))
})

test_that("make_metadata works", {
    res <- make_metadata(source = "some", source_version = 1,
                         source_date = "now", url = "some url",
                         organism = "Hsapiens")
    expect_true(is(res, "data.frame"))
    expect_true(.valid_metadata(res))
    expect_equal(res[, "value"], c("some", "some url", "1", "now",
                                   "Hsapiens"))
    ## Errors...
    expect_error(make_metadata())
    expect_error(make_metadata(source = "a", source_version = "2",
                               source_date = "now", organism = "MM"))
    expect_error(make_metadata(source = "a", source_version = "2", url = NULL,
                               source_date = "now", organism = "MM"))
})

test_that(".valid_msms_spectrum works", {
    msms_spctra_local <- msms_spctra
    msms_spctra_local$collision_energy <-
        as.character(msms_spctra_local$collision_energy)
    expect_true(.valid_msms_spectrum(msms_spctra_local,
                                     blob = is.list(msms_spctra_local$mz)))
    tmp <- msms_spctra
    expect_error(.valid_msms_spectrum(tmp[, 1:4]))
    tmp$polarity <- as.character(tmp$polarity)
    expect_error(.valid_msms_spectrum(tmp))
    res <- .valid_msms_spectrum(tmp, error = FALSE)
    expect_true(is.character(res))
    dr <- system.file("xml/", package = "CompoundDb")
    tmp <- msms_spectra_hmdb(dr, collapsed = FALSE)
    tmp$collision_energy <- as.character(tmp$collision_energy)
    expect_true(.valid_msms_spectrum(tmp, blob = FALSE))
})

test_that(".insert_msms_spectra works", {
    library(RSQLite)
    con <- dbConnect(dbDriver("SQLite"), dbname = tempfile())
    .insert_msms_spectra(con, msms_spctra)
    res <- dbGetQuery(con, "select * from msms_spectrum")
    expect_equal(nrow(res), nrow(msms_spctra))
    expect_equal(res$splash, msms_spctra$splash)
    res_2 <- dbGetQuery(con, "select * from msms_spectrum_peak")
    expect_equal(colnames(res_2), c("spectrum_id", "mz", "intensity", "peak_id"))
    expect_equal(res_2$mz, unlist(msms_spctra$mz))
    expect_equal(res_2$intensity, unlist(msms_spctra$intensity))
})

test_that(".add_mz_range_column works", {
    z <- msms_spctra
    res <- .add_mz_range_column(z)
    expect_true(ncol(res) == ncol(z) + 2)
    expect_true(all(c("msms_mz_range_min", "msms_mz_range_max") %in% colnames(res)))

    colnames(z)[colnames(z) == "mz"] <- "other"
    expect_error(.add_mz_range_column(z))

    ## Expanded
    z <- .expand_spectrum_df(msms_spctra)
    res <- .add_mz_range_column(z)
    mzs <- split(res$mz, res$spectrum_id)
    mzmin <- split(res$msms_mz_range_min, res$spectrum_id)
    mzmax <- split(res$msms_mz_range_max, res$spectrum_id)
    expect_equal(lapply(mzs, min), lapply(mzmin, min))
    expect_equal(lapply(mzs, max), lapply(mzmax, max))
})

test_that("import_mona_sdf works", {
    mona <- system.file("sdf/MoNa_export-All_Spectra_sub.sdf.gz",
                        package = "CompoundDb")
    chebi <- system.file("sdf/ChEBI_sub.sdf.gz", package = "CompoundDb")

    res <- .import_mona_sdf(mona)
    expect_equal(names(res), c("compound", "msms_spectrum"))
    expect_equal(nrow(res$compound), 7)
    expect_equal(nrow(res$msms_spectrum), 7)

    expect_equal(import_mona_sdf(mona), res)

    res <- .import_mona_sdf(mona, TRUE, FALSE)
    expect_true(nrow(res$compound) == 0)

    expect_error(import_mona_sdf(chebi))
})

test_that(".msms_spectrum_add_missing_columns works", {
    df <- data.frame(spectrum_id = 1:3, polarity = 1L)
    res <- .msms_spectrum_add_missing_columns(df)
    expect_true(all(c("collision_energy", "predicted", "splash") %in%
                    colnames(res)))
    expect_true(is.character(res$instrument))
    expect_equal(res$polarity, df$polarity)
    expect_equal(res$spectrum_id, df$spectrum_id)
})

test_that(".append_msms_spectra works", {
    tmp_con <- dbConnect(SQLite(), tempfile())
    tmp_db <- .copy_compdb(.dbconn(cmp_spctra_db), tmp_con)
    tmp_db <- CompDb(tmp_con)

    spd <- DataFrame(
        msLevel = c(2L, 2L),
        polarity = c(1L, 1L),
        compound_id = c("HMDB0000008", "HMDB0000008"))
    spd$mz <- list(
        c(109.2, 124.2, 124.5, 170.16, 170.52),
        c(83.1, 96.12, 97.14, 109.14, 124.08, 125.1, 170.16))
    spd$intensity <- list(
        c(3.407, 47.494, 3.094, 100.0, 13.240),
        c(6.685, 4.381, 3.022, 16.708, 100.0, 4.565, 40.643))
    sps <- Spectra(spd)

    x <- as.data.frame(spectraData(sps, c("msLevel", "polarity", "precursorMz",
                                          "rtime", "compound_id")))
    x$mz <- as.list(sps$mz)
    x$intensity <- as.list(sps$intensity)
    x$spectrum_id <- seq_len(nrow(x))

    expect_warning(.append_msms_spectra(tmp_con, x), "replaced with internal")
    res <- dbGetQuery(tmp_con, "select * from msms_spectrum")
    expect_true(all(c("msLevel", "polarity", "precursorMz", "rtime") %in%
                    colnames(res)))
    expect_true(all(res$compound_id == c("HMDB0000001", "HMDB0000001",
                                         "HMDB0004370", "HMDB0006719",
                                         "HMDB0000008", "HMDB0000008")))
    res <- dbGetQuery(tmp_con, "select * from msms_spectrum_peak")
    expect_true(sum(res$spectrum_id %in% c(5, 6)) == 12)
})

test_that(".prepare_msms_spectra_table works", {
    tbl <- data.frame(spectrum_id = 1:6, collision_energy = 5, rtime = 1:6,
                      compound_id = "a", polarity = 1L, predicted = FALSE,
                      precursor_mz = as.numeric(1:6))
    tbl$mz <- list(1:2, 4, 1:6, 4, 40:50, 12:15)
    tbl$intensity <- list(1:2, 4, 1:6, 4, 40:50, 12:15)
    res <- .prepare_msms_spectra_table(tbl)
    expect_equal(names(res), c("x", "msms_spectrum_peak"))
    expect_true(is.character(res$x$splash))

    tbl$mz[[2]] <- c("a")
    expect_error(.prepare_msms_spectra_table(tbl), "contain numeric values")
})

test_that(".valid_data_frame_columns works", {
    tbl <- data.frame(a = 1:4, b = "B")
    expect_match(.valid_data_frame_columns("b"), "data.frame")

    expect_equal(.valid_data_frame_columns(tbl, req_cols = c("a", "b")),
                 character())
    expect_match(.valid_data_frame_columns(tbl, req_cols = c("a", "b", "c")),
                 "required")
})

test_that(".throw_error works", {
    expect_true(.throw_error())
    expect_error(.throw_error("abc"), "abc")
    expect_match(.throw_error("abc", error = FALSE), "abc")
})

test_that("emptyCompDb works", {
    fl <- tempfile()
    res <- emptyCompDb(fl)
    expect_s4_class(res, "CompDb")
    expect_true(nrow(compounds(res)) == 0)

    expect_error(emptyCompDb(fl), "exist")
})

test_that(".parse_lipidblast_json_element works", {
    library(jsonlite)
    f <- system.file("json", "MoNa-LipidBlast_sub.json", package = "CompoundDb")
    js <- read_json(f)
    res <- .parse_lipidblast_json_element(js[[1L]])
    expect_true(is.list(res))
    expect_true(all(c("compound_id", "name", "inchi", "inchikey",
                      "formula", "exactmass", "synonyms") %in% names(res)))
    expect_equal(res$name, "CerP 24:0")

    x <- js[[1L]]
    x$spectrum <- NULL

    res <- CompoundDb:::.parse_lipidblast_json_element(x)
    expect_identical(res$spectrum, NA_character_)
})

test_that(".import_lipidblast_json_chunk works", {
    f <- system.file("json", "MoNa-LipidBlast_sub.json", package = "CompoundDb")
    res <- .import_lipidblast_json_chunk(f, n = 3)
    expect_true(is.list(res))
    expect_true(length(res) == 8L)

    res <- .import_lipidblast_json_chunk(f, n = 9, verbose = TRUE)
    expect_true(is.list(res))
    expect_true(length(res) == 8L)

    ref <- .import_lipidblast(f, verbose = TRUE)
    expect_equal(nrow(ref), length(res))
    res <- bind_rows(res)
    expect_equal(ref, res)
})

test_that(".lipidblast_parallel_parse works", {
    f <- system.file("json", "MoNa-LipidBlast_sub.json", package = "CompoundDb")
    l <- readLines(f)
    l <- sub(",$", "", l)
    res <- .lipidblast_parallel_parse(l[2:5], BPPARAM = SerialParam())
    expect_true(length(res) == 4)
})

test_that("compound_tbl_lipidblast works with n > 0", {
    f <- system.file("json", "MoNa-LipidBlast_sub.json", package = "CompoundDb")
    ref <- compound_tbl_lipidblast(f)
    expect_true(is.numeric(ref$exactmass))
    expect_true(is.character(ref$compound_id))
    expect_true(is.character(ref$name))
    expect_true(is.character(ref$inchi))
    expect_true(is.character(ref$inchikey))
    expect_true(is.character(ref$formula))
    expect_true(is.list(ref$synonyms))
    res <- compound_tbl_lipidblast(f, n = 4)
    expect_true(is.numeric(res$exactmass))
    expect_equal(ref, res)
})
