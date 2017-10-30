test_that(".simple_import_compounds_sdf works", {
    hmdb <- system.file("sdf/HMDB_sub.sdf", package = "CompoundDb")
    cmps <- .simple_import_compounds_sdf(hmdb)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass", "synonyms"))
    expect_true(nrow(cmps) == 7)

    chebi <- system.file("sdf/ChEBI_sub.sdf", package = "CompoundDb")
    cmps <- .simple_import_compounds_sdf(chebi)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass", "synonyms"))
    expect_true(nrow(cmps) == 6)

    lm <- system.file("sdf/LipidMaps_sub.sdf", package = "CompoundDb")
    cmps <- .simple_import_compounds_sdf(lm)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass", "synonyms"))
    expect_true(nrow(cmps) == 7)
    
    pubchem <- system.file("sdf/PubChem_sub.sdf", package = "CompoundDb")
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
    expect_true(nrow(cmps) == 7)
    expect_true(is(cmps$synonyms, "list"))
    cmps <- compound_tbl_sdf(hmdb, collapse = "|")
    expect_true(is(cmps$synonyms, "character"))

    chebi <- system.file("sdf/ChEBI_sub.sdf", package = "CompoundDb")
    cmps <- compound_tbl_sdf(chebi)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass", "synonyms"))
    expect_true(nrow(cmps) == 6)
    expect_true(is(cmps$synonyms, "list"))
    cmps <- compound_tbl_sdf(chebi, collapse = "|")
    expect_true(is(cmps$synonyms, "character"))

    lm <- system.file("sdf/LipidMaps_sub.sdf", package = "CompoundDb")
    cmps <- compound_tbl_sdf(lm)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass", "synonyms"))
    expect_true(nrow(cmps) == 7)
    expect_true(is(cmps$synonyms, "list"))
    cmps <- compound_tbl_sdf(lm, collapse = "|")
    expect_true(is(cmps$synonyms, "character"))

    pc <- system.file("sdf/PubChem_sub.sdf", package = "CompoundDb")
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
    
    chebi <- system.file("sdf/ChEBI_sub.sdf", package = "CompoundDb")
    cn <- colnames(datablock2ma(datablock(read.SDFset(chebi))))
    expect_true(.guess_sdf_source(cn) == "chebi")

    lm <- system.file("sdf/LipidMaps_sub.sdf", package = "CompoundDb")
    cn <- colnames(datablock2ma(datablock(read.SDFset(lm))))
    expect_true(.guess_sdf_source(cn) == "lipidmaps")

    pc <- system.file("sdf/PubChem_sub.sdf", package = "CompoundDb")
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
})
