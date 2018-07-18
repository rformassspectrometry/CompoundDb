test_that(".simple_import_compounds_sdf works", {
    hmdb <- system.file("sdf/HMDB_sub.sdf", package = "CompoundDb")
    cmps <- .simple_import_compounds_sdf(hmdb)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass", "synonyms"))
    expect_true(nrow(cmps) == 9)

    chebi <- system.file("sdf/ChEBI_sub.sdf.gz", package = "CompoundDb")
    cmps <- CompoundDb:::.simple_import_compounds_sdf(chebi)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass", "synonyms"))
    expect_true(nrow(cmps) == 6)

    lm <- system.file("sdf/LipidMaps_sub.sdf.gz", package = "CompoundDb")
    cmps <- .simple_import_compounds_sdf(lm)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass", "synonyms"))
    expect_true(nrow(cmps) == 7)
    
    pubchem <- system.file("sdf/PubChem_sub.sdf.gz", package = "CompoundDb")
    cmps <- .simple_import_compounds_sdf(pubchem)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass", "synonyms"))
    expect_true(nrow(cmps) == 12)
})

test_that("compound_tbl_sdf works", {
    hmdb <- system.file("sdf/HMDB_sub.sdf", package = "CompoundDb")
    cmps <- compound_tbl_sdf(hmdb)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass", "synonyms"))
    expect_true(nrow(cmps) == 9)
    expect_true(is(cmps$synonyms, "list"))
    cmps <- compound_tbl_sdf(hmdb, collapse = "|")
    expect_true(is(cmps$synonyms, "character"))

    chebi <- system.file("sdf/ChEBI_sub.sdf.gz", package = "CompoundDb")
    cmps <- compound_tbl_sdf(chebi)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass", "synonyms"))
    expect_true(nrow(cmps) == 6)
    expect_true(is(cmps$synonyms, "list"))
    cmps <- compound_tbl_sdf(chebi, collapse = "|")
    expect_true(is(cmps$synonyms, "character"))

    lm <- system.file("sdf/LipidMaps_sub.sdf.gz", package = "CompoundDb")
    cmps <- compound_tbl_sdf(lm)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass", "synonyms"))
    expect_true(nrow(cmps) == 7)
    expect_true(is(cmps$synonyms, "list"))
    cmps <- compound_tbl_sdf(lm, collapse = "|")
    expect_true(is(cmps$synonyms, "character"))

    pc <- system.file("sdf/PubChem_sub.sdf.gz", package = "CompoundDb")
    cmps <- compound_tbl_sdf(pc)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass", "synonyms"))
    expect_true(nrow(cmps) == 12)
    expect_true(is(cmps$synonyms, "list"))
    cmps <- compound_tbl_sdf(pc, collapse = "|")
    expect_true(is(cmps$synonyms, "character"))
})

test_that(".guess_sdf_source works", {
    hmdb <- system.file("sdf/HMDB_sub.sdf", package = "CompoundDb")
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
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass", "synonyms"))
    expect_true(nrow(cmps) == 8)
})

test_that("compound_tbl_lipidblast works", {
    lb <- system.file("json/MoNa-LipidBlast_sub.json", package = "CompoundDb")
    cmps <- compound_tbl_lipidblast(lb)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass", "synonyms"))
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
    cmps <- data.frame(compound_id = c("01", "02"), compound_name = c("a", "b"),
                       inchi = c("i1", "i2"), formula = c("some", "thing"),
                       mass = c(1, 3), synonyms = c("a", "b"))
    expect_true(.valid_compound(cmps, db = FALSE))
    expect_true(.valid_compound(cmps[, 1:5]))
    expect_error(.valid_compound(cmps[, 1:5], db = FALSE))
    expect_error(.valid_compound(cmps, db = TRUE))

    ## Errors
    expect_error(.valid_compound("b"))
    expect_true(is.character(.valid_compound("b", error = FALSE)))
    expect_error(.valid_compound(data.frame()))
    expect_error(.valid_compound(cmps[, 1:3]))
    cmps$mass <- c("1", "2")
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

    ## Provide a single file name.
    fl <- system.file("sdf/LipidMaps_sub.sdf.gz", package = "CompoundDb")
    metad <- make_metadata(source = "LipidMaps", source_date = "2016",
                           source_version = "xx", organism = "Hsapiens",
                           url = NA)
    res <- createCompDb(fl, metadata = metad, path = tempdir())
    db <- CompDb(res)
    md <- CompoundDb:::.metadata(db)
    expect_true(is.na(md$value[md$name == "url"]))
    ## Multiple files.
    fls <- c(system.file("sdf/LipidMaps_sub.sdf.gz", package = "CompoundDb"),
             system.file("sdf/HMDB_sub.sdf", package = "CompoundDb"),
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

test_that(".check_msms_spectra_columns works", {
    expect_true(.check_msms_spectra_columns(msms_spctra))
    tmp <- msms_spctra
    expect_error(.check_msms_spectra_columns(tmp[, 1:4]))
    tmp$polarity <- as.character(tmp$polarity)
    expect_error(.check_msms_spectra_columns(tmp))
})
