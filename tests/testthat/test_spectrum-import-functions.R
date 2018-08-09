test_that(".import_hmdb_ms_ms_spectrum works", {
    
    fl <- system.file("xml/HMDB0000001_ms_ms_spectrum_2_experimental.xml",
                      package = "CompoundDb")

    expect_error(.import_hmdb_ms_ms_spectrum())
    expect_error(.import_hmdb_ms_ms_spectrum(4))

    library(xml2)
    x <- read_xml(fl)
    res <- CompoundDb:::.import_hmdb_ms_ms_spectrum(fl, collapsed = FALSE)
    expect_equal(colnames(res), c("spectrum_id", "compound_id", "polarity",
                                  "collision_energy", "predicted", "splash",
                                  "instrument_type", "mz", "intensity"))
    expect_equal(nrow(res), 7)

    expect_equal(res$compound_id[1], xml_text(xml_find_first(x, "database-id")))
    expect_equal(res$spectrum_id[1], xml_text(xml_find_first(x, "id")))

    ## Checking collapsed input.
    res_clpsd <- CompoundDb:::.import_hmdb_ms_ms_spectrum(fl, collapsed = TRUE)
    expect_equal(nrow(res_clpsd), 1)
    expect_equal(res_clpsd$compound_id, res$compound_id[1])
    expect_equal(res$mz, unlist(res_clpsd$mz))
    expect_equal(res$intensity, unlist(res_clpsd$intensity))

    ## One that should fail.
    fl <- system.file("xml/fail/HMDB0001875_ms_ms_spectrum_1768_experimental.xml",
                      package = "CompoundDb")
    expect_error(CompoundDb:::.import_hmdb_ms_ms_spectrum(fl))
    expect_warning(CompoundDb:::.import_hmdb_ms_ms_spectrum(fl, nonStop = TRUE))

})

test_that("msms_spectra_hmdb works", {
    dr <- system.file("xml", package = "CompoundDb")
    res <- msms_spectra_hmdb(dr, collapsed = FALSE)
    expect_true(length(unique(res$spectrum_id)) == 4)
    expect_equal(colnames(res), c("original_spectrum_id", "compound_id",
                                  "polarity", "collision_energy", "predicted",
                                  "splash", "instrument_type", "mz",
                                  "intensity", "spectrum_id"))
    ## Get it in collapsed form.
    res_clpsd <- msms_spectra_hmdb(dr, collapsed = TRUE)
    expect_equal(colnames(res_clpsd),
                 c("original_spectrum_id", "compound_id", "polarity",
                   "collision_energy", "predicted", "splash",
                   "instrument_type", "mz", "intensity", "spectrum_id"))
    expect_true(nrow(res_clpsd) == 4)
    expect_equal(unlist(res_clpsd$mz), res$mz)
    expect_equal(unlist(res_clpsd$intensity), res$intensity)
    
    ## Construct the Spectrum2List of these
    spl <- Spectrum2List(res)
    expect_equal(sapply(spl, polarity), c(1, 1, 1, 0))
    
    dr <- system.file("sdf", package = "CompoundDb")
    expect_error(msms_spectra_hmdb(dr))
})

