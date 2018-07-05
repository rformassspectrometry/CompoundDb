test_that(".import_hmdb_ms_ms_spectrum works", {
    
    fl <- system.file("xml/HMDB0000001_ms_ms_spectrum_2_experimental.xml",
                      package = "CompoundDb")

    expect_error(.import_hmdb_ms_ms_spectrum())
    expect_error(.import_hmdb_ms_ms_spectrum(4))

    x <- read_xml(fl)
    res <- .import_hmdb_ms_ms_spectrum(fl)
    expect_equal(colnames(res), c("spectrum_id", "compound_id", "polarity",
                                  "collision_energy", "predicted", "splash",
                                  "mz", "intensity"))
    expect_equal(nrow(res), 7)

    expect_equal(res$compound_id[1], xml_text(xml_find_first(x, "database-id")))
    expect_equal(res$spectrum_id[1], xml_text(xml_find_first(x, "id")))    
})
