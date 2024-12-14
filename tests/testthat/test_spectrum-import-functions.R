test_that(".import_hmdb_ms_ms_spectrum works", {
    fl <- system.file("xml/fail/no_id.xml", package = "CompoundDb")
    expect_error(CompoundDb:::.import_hmdb_ms_ms_spectrum(fl), "Could not extract")
    expect_warning(res <- .import_hmdb_ms_ms_spectrum(fl, nonStop = TRUE), "Could not extract")
    expect_equal(res, data.frame())

    fl <- system.file("xml/fail/missing_values.xml", package = "CompoundDb")
    res <- CompoundDb:::.import_hmdb_ms_ms_spectrum(fl, nonStop = TRUE)
    expect_equal(res$spectrum_id, as.character(7))
    expect_true(is.na(res$polarity))
    expect_true(is.na(res$predicted))
    expect_true(is.na(res$instrument))
    expect_true(is.na(res$instrument_type))

    fl <- system.file("xml/HMDB0000001_ms_ms_spectrum_2_experimental.xml",
                      package = "CompoundDb")

    expect_error(.import_hmdb_ms_ms_spectrum())
    expect_error(.import_hmdb_ms_ms_spectrum(4))

    library(xml2)
    x <- read_xml(fl)
    res <- .import_hmdb_ms_ms_spectrum(fl, collapsed = FALSE)
    expect_equal(colnames(res), c("spectrum_id", "compound_id", "polarity",
                                  "collision_energy", "predicted", "splash",
                                  "instrument_type", "instrument",
                                  "precursor_mz", "mz", "intensity"))
    expect_equal(nrow(res), 7)

    expect_equal(res$compound_id[1], xml_text(xml_find_first(x, "database-id")))
    expect_equal(res$spectrum_id[1], xml_text(xml_find_first(x, "id")))

    ## Checking collapsed input.
    res_clpsd <- .import_hmdb_ms_ms_spectrum(fl, collapsed = TRUE)
    expect_equal(nrow(res_clpsd), 1)
    expect_equal(res_clpsd$compound_id, res$compound_id[1])
    expect_equal(res$mz, unlist(res_clpsd$mz))
    expect_equal(res$intensity, unlist(res_clpsd$intensity))

    ## One that should fail.
    fl <- system.file("xml/fail/HMDB0001875_ms_ms_spectrum_1768_experimental.xml",
                      package = "CompoundDb")
    expect_error(.import_hmdb_ms_ms_spectrum(fl))
    expect_warning(.import_hmdb_ms_ms_spectrum(fl, nonStop = TRUE))
})

test_that("msms_spectra_hmdb works", {
    dr <- system.file("xml", package = "CompoundDb")
    res <- msms_spectra_hmdb(dr, collapsed = FALSE)
    expect_true(length(unique(res$spectrum_id)) == 4)
    expect_equal(colnames(res), c("original_spectrum_id", "compound_id",
                                  "polarity", "collision_energy", "predicted",
                                  "splash", "instrument_type", "instrument",
                                  "precursor_mz", "mz", "intensity",
                                  "spectrum_id"))
    ## Get it in collapsed form.
    res_clpsd <- msms_spectra_hmdb(dr, collapsed = TRUE)
    expect_equal(colnames(res_clpsd),
                 c("original_spectrum_id", "compound_id", "polarity",
                   "collision_energy", "predicted", "splash",
                   "instrument_type", "instrument", "precursor_mz",
                   "mz", "intensity", "spectrum_id"))
    expect_true(nrow(res_clpsd) == 4)
    expect_equal(unlist(res_clpsd$mz), res$mz)
    expect_equal(unlist(res_clpsd$intensity), res$intensity)

    dr <- system.file("sdf", package = "CompoundDb")
    expect_error(msms_spectra_hmdb(dr))
})

## test_that(".spectra2_from_df works", {
##     sp1 <- new("Spectrum2", mz = c(1, 2, 4), intensity = c(4, 5, 2))
##     sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5))
##     ## Build from data.frame
##     df <- data.frame(spectrum_id = c("b", "b", "b", "b", "a", "a", "a"),
##                      mz = c(1, 2, 3, 4, 1, 2, 4),
##                      intensity = c(5, 3, 2, 5, 4, 5, 2),
##                      polarity = c(1, 1, 1, 1, 0, 0, 0),
##                      compound_id = rep("cp_1", 7), stringsAsFactors = FALSE)
##     ## Test internal function:
##     res <- .spectra2_from_df(df)
##     expect_equal(mz(res$spectra[[1]]), mz(sp2))
##     expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
##     expect_equal(mz(res$spectra[[2]]), mz(sp1))
##     expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
##     expect_equal(polarity(res$spectra[[1]]), 1)
##     expect_equal(polarity(res$spectra[[2]]), 0)
##     expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
##     expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

##     ## provide rt
##     df$rt <- c(1, 1, 1, 1, 3, 3, 3)
##     res <- .spectra2_from_df(df)
##     expect_equal(rtime(res$spectra[[1]]), 1)
##     expect_equal(rtime(res$spectra[[2]]), 3)
##     expect_equal(mz(res$spectra[[1]]), mz(sp2))
##     expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
##     expect_equal(mz(res$spectra[[2]]), mz(sp1))
##     expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
##     expect_equal(polarity(res$spectra[[1]]), 1)
##     expect_equal(polarity(res$spectra[[2]]), 0)
##     expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
##     expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

##     ## provide ms_level
##     df$ms_level <- c(2, 2, 2, 2, 3, 3, 3)
##     res <- .spectra2_from_df(df)
##     expect_equal(msLevel(res$spectra[[1]]), 2)
##     expect_equal(msLevel(res$spectra[[2]]), 3)
##     expect_equal(mz(res$spectra[[1]]), mz(sp2))
##     expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
##     expect_equal(mz(res$spectra[[2]]), mz(sp1))
##     expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
##     expect_equal(polarity(res$spectra[[1]]), 1)
##     expect_equal(polarity(res$spectra[[2]]), 0)
##     expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
##     expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

##     ## provide precursor_mz
##     df$precursor_mz <- c(2.1, 2.1, 2.1, 2.1, 3.4, 3.4, 3.4)
##     res <- .spectra2_from_df(df)
##     expect_equal(precursorMz(res$spectra[[1]]), 2.1)
##     expect_equal(precursorMz(res$spectra[[2]]), 3.4)
##     expect_equal(mz(res$spectra[[1]]), mz(sp2))
##     expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
##     expect_equal(mz(res$spectra[[2]]), mz(sp1))
##     expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
##     expect_equal(polarity(res$spectra[[1]]), 1)
##     expect_equal(polarity(res$spectra[[2]]), 0)
##     expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
##     expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

##     ## provide precursor_charge
##     df$precursor_charge <- c(5, 5, 5, 5, 2, 2, 2)
##     res <- .spectra2_from_df(df)
##     expect_equal(precursorCharge(res$spectra[[1]]), 5)
##     expect_equal(precursorCharge(res$spectra[[2]]), 2)
##     expect_equal(mz(res$spectra[[1]]), mz(sp2))
##     expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
##     expect_equal(mz(res$spectra[[2]]), mz(sp1))
##     expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
##     expect_equal(polarity(res$spectra[[1]]), 1)
##     expect_equal(polarity(res$spectra[[2]]), 0)
##     expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
##     expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

##     ## provide precuror_intensity
##     df$precursor_intensity <- c(3, 3, 3, 3, 1, 1, 1)
##     res <- .spectra2_from_df(df)
##     expect_equal(precursorIntensity(res$spectra[[1]]), 3)
##     expect_equal(precursorIntensity(res$spectra[[2]]), 1)
##     expect_equal(mz(res$spectra[[1]]), mz(sp2))
##     expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
##     expect_equal(mz(res$spectra[[2]]), mz(sp1))
##     expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
##     expect_equal(polarity(res$spectra[[1]]), 1)
##     expect_equal(polarity(res$spectra[[2]]), 0)
##     expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
##     expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

##     ## provide collision_energy
##     df$collision_energy <- c(6, 6, 6, 6, 3, 3, 3)
##     res <- .spectra2_from_df(df)
##     expect_equal(collisionEnergy(res$spectra[[1]]), 6)
##     expect_equal(collisionEnergy(res$spectra[[2]]), 3)
##     expect_equal(mz(res$spectra[[1]]), mz(sp2))
##     expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
##     expect_equal(mz(res$spectra[[2]]), mz(sp1))
##     expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
##     expect_equal(polarity(res$spectra[[1]]), 1)
##     expect_equal(polarity(res$spectra[[2]]), 0)
##     expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
##     expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

##     ## provide acquisition_num
##     df$acquisition_num <- c(1, 1, 1, 1, 2, 2, 2)
##     res <- .spectra2_from_df(df)
##     expect_equal(acquisitionNum(res$spectra[[1]]), 1)
##     expect_equal(acquisitionNum(res$spectra[[2]]), 2)
##     expect_equal(mz(res$spectra[[1]]), mz(sp2))
##     expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
##     expect_equal(mz(res$spectra[[2]]), mz(sp1))
##     expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
##     expect_equal(polarity(res$spectra[[1]]), 1)
##     expect_equal(polarity(res$spectra[[2]]), 0)
##     expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
##     expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

##     ## provide scan_index
##     df$scan_index <- c(1, 1, 1, 1, 2, 2, 2)
##     res <- .spectra2_from_df(df)
##     expect_equal(scanIndex(res$spectra[[1]]), 1)
##     expect_equal(scanIndex(res$spectra[[2]]), 2)
##     expect_equal(mz(res$spectra[[1]]), mz(sp2))
##     expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
##     expect_equal(mz(res$spectra[[2]]), mz(sp1))
##     expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
##     expect_equal(polarity(res$spectra[[1]]), 1)
##     expect_equal(polarity(res$spectra[[2]]), 0)
##     expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
##     expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

##     ## provide from_file
##     df$from_file <- c(1, 1, 1, 1, 2, 2, 2)
##     res <- .spectra2_from_df(df)
##     expect_equal(fromFile(res$spectra[[1]]), 1)
##     expect_equal(fromFile(res$spectra[[2]]), 2)
##     expect_equal(mz(res$spectra[[1]]), mz(sp2))
##     expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
##     expect_equal(mz(res$spectra[[2]]), mz(sp1))
##     expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
##     expect_equal(polarity(res$spectra[[1]]), 1)
##     expect_equal(polarity(res$spectra[[2]]), 0)
##     expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
##     expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

##     ## provide precursor_scan_num
##     df$precursor_scan_num <- c(3, 3, 3, 3, 4, 4, 4)
##     res <- .spectra2_from_df(df)
##     expect_equal(precScanNum(res$spectra[[1]]), 3)
##     expect_equal(precScanNum(res$spectra[[2]]), 4)
##     expect_equal(mz(res$spectra[[1]]), mz(sp2))
##     expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
##     expect_equal(mz(res$spectra[[2]]), mz(sp1))
##     expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
##     expect_equal(polarity(res$spectra[[1]]), 1)
##     expect_equal(polarity(res$spectra[[2]]), 0)
##     expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
##     expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))
## })

test_that("msms_spectra_mona works", {
    fl <- system.file("sdf/MoNa_export-All_Spectra_sub.sdf.gz",
                      package = "CompoundDb")
    res <- msms_spectra_mona(fl)
    expect_equal(nrow(res), 7)
    res <- msms_spectra_mona(fl, collapsed = FALSE)
    expect_equal(nrow(res), 179)
})

test_that(".extract_spectra_mona_sdf works", {
    fl <- system.file("sdf/MoNa_export-All_Spectra_sub.sdf.gz",
                      package = "CompoundDb")
    x <- datablock2ma(datablock(read.SDFset(fl)))
    res <- .extract_spectra_mona_sdf(x)
    expect_true(all(res$polarity == 1L))
    expect_true(is.character(res$collision_energy))
    expect_true(is.numeric(res$precursor_mz))
})

test_that(".compound_id_from_mona_sdf works", {
    res <- .compound_id_from_mona_sdf(matrix(nrow = 4))
    expect_equal(length(res), 4)
    expect_equal(res, c("CMP1", "CMP2", "CMP3", "CMP4"))
})
