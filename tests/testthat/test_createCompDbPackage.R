test_that(".simple_extract_compounds_sdf works", {

    hmdb <- system.file("sdf/HMDB_sub.sdf.gz", package = "CompoundDb")
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
    cmps <- .simple_extract_compounds_sdf(
        datablock2ma(datablock(read.SDFset(lm))))
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "name", "inchi",
                                   "inchikey", "formula", "exactmass",
                                   "synonyms", "smiles"))
    expect_true(nrow(cmps) == 7)

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
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "name", "inchi",
                                   "inchikey", "formula", "exactmass",
                                   "synonyms", "smiles"))
    expect_true(nrow(cmps) == 9)
    expect_true(is(cmps$synonyms, "list"))
    cmps <- compound_tbl_sdf(hmdb, collapse = "|")
    expect_true(is(cmps$synonyms, "character"))

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
    expect_equal(colnames(cmps), c("compound_id", "name", "inchi",
                                   "inchikey", "formula", "exactmass",
                                   "synonyms"))
    expect_true(nrow(cmps) == 8)
})

test_that("compound_tbl_lipidblast works", {
    expect_error(compound_tbl_lipidblast())
    expect_error(compound_tbl_lipidblast("sddfd"))

    lb <- system.file("json/MoNa-LipidBlast_sub.json", package = "CompoundDb")
    cmps <- compound_tbl_lipidblast(lb)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "name", "inchi",
                                   "inchikey", "formula", "exactmass",
                                   "synonyms"))
    expect_true(nrow(cmps) == 8)
    expect_true(is(cmps$synonyms, "character"))
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
    expect_error(.valid_metadata(metadata[1:4, ]))
    metadata_fail <- metadata

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

    res <- import_mona_sdf(mona)
    expect_equal(names(res), c("compound", "msms_spectrum"))
    expect_equal(nrow(res$compound), 7)
    expect_equal(nrow(res$msms_spectrum), 7)

    expect_error(import_mona_sdf(chebi))
})
