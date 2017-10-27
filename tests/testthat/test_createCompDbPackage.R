test_that(".simple_import_compounds_sdf works", {
    hmdb <- system.file("sdf/HMDB_sub.sdf", package = "CompoundDb")
    cmps <- .simple_import_compounds_sdf(hmdb)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass"))
    expect_true(nrow(cmps) == 7)

    chebi <- system.file("sdf/ChEBI_sub.sdf", package = "CompoundDb")
    cmps <- .simple_import_compounds_sdf(chebi)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass"))
    expect_true(nrow(cmps) == 6)

    lm <- system.file("sdf/LipidMaps_sub.sdf", package = "CompoundDb")
    cmps <- .simple_import_compounds_sdf(lm)
    expect_true(is(cmps, "data.frame"))
    expect_true(is(cmps, "tbl"))
    expect_equal(colnames(cmps), c("compound_id", "compound_name", "inchi",
                                   "formula", "mass"))
    expect_true(nrow(cmps) == 7)
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

    expect_equal(.guess_sdf_source(c("some", "thing")), NULL)
})

