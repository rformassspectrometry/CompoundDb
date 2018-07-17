test_that(".import_hmdb_ms_ms_spectrum works", {
    
    fl <- system.file("xml/HMDB0000001_ms_ms_spectrum_2_experimental.xml",
                      package = "CompoundDb")

    expect_error(.import_hmdb_ms_ms_spectrum())
    expect_error(.import_hmdb_ms_ms_spectrum(4))

    library(xml2)
    x <- read_xml(fl)
    res <- CompoundDb:::.import_hmdb_ms_ms_spectrum(fl)
    expect_equal(colnames(res), c("spectrum_id", "compound_id", "polarity",
                                  "collision_energy", "predicted", "splash",
                                  "instrument_type", "mz", "intensity"))
    expect_equal(nrow(res), 7)

    expect_equal(res$compound_id[1], xml_text(xml_find_first(x, "database-id")))
    expect_equal(res$spectrum_id[1], xml_text(xml_find_first(x, "id")))

    ## One that should fail.
    fl <- system.file("xml/fail/HMDB0001875_ms_ms_spectrum_1768_experimental.xml",
                      package = "CompoundDb")
    expect_error(CompoundDb:::.import_hmdb_ms_ms_spectrum(fl))
    expect_warning(CompoundDb:::.import_hmdb_ms_ms_spectrum(fl, nonStop = TRUE))
})

test_that("msms_spectra_hmdb works", {
    dr <- system.file("xml", package = "CompoundDb")
    res <- msms_spectra_hmdb(dr)
    expect_true(length(unique(res$spectrum_id)) == 4)
    expect_equal(colnames(res), c("spectrum_id", "compound_id", "polarity",
                                  "collision_energy", "predicted", "splash",
                                  "instrument_type", "mz", "intensity"))

    ## Construct the Spectrum2List of these
    spl <- Spectrum2List(res)
    expect_equal(sapply(spl, polarity), c(1, 1, 1, 0))
    
    dr <- system.file("sdf", package = "CompoundDb")
    expect_error(msms_spectra_hmdb(dr))
})

