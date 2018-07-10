test_that("Spectrum2List construction works as expected", {
    sp1 <- new("Spectrum2", mz = c(1, 2, 4), intensity = c(4, 5, 2))
    sp2 <- new("Spectrum2", mz = c(1, 2, 3, 4), intensity = c(5, 3, 2, 5))

    ## Errors.
    expect_error(new("Spectrum2List", 4))
    expect_error(new("Spectrum2List", list(4)))
    expect_error(new("Spectrum2List", list(sp1, 4)))
    expect_error(Spectrum2List(4))

    spl <- new("Spectrum2List", list(sp1, sp2))
    expect_true(validObject(spl))
    expect_equal(spl[[1]], sp1)
    expect_equal(spl[[2]], sp2)

    spl <- Spectrum2List(sp1)
    expect_true(validObject(spl))
    expect_equal(spl[[1]], sp1)
    spl <- Spectrum2List(sp1, sp2)
    expect_equal(spl[[1]], sp1)
    expect_equal(spl[[2]], sp2)
    spl <- Spectrum2List(list(sp1, sp2))
    expect_equal(spl[[1]], sp1)
    expect_equal(spl[[2]], sp2)
    
    ## Concatenating.
    spl <- c(Spectrum2List(sp1), Spectrum2List(sp2))
    expect_true(is(spl, "Spectrum2List"))
    expect_equal(spl[[1]], sp1)
    expect_equal(spl[[2]], sp2)

    ## Build from data.frame
    df <- data.frame(spectrum_id = c("b", "b", "b", "b", "a", "a", "a"),
                     mz = c(1, 2, 3, 4, 1, 2, 4),
                     intensity = c(5, 3, 2, 5, 4, 5, 2),
                     polarity = c(1, 1, 1, 1, 0, 0, 0),
                     compound_id = rep("cp_1", 7), stringsAsFactors = FALSE)
    ## Test internal function:
    res <- CompoundDb:::.spectra2_from_df(df)
    expect_equal(mz(res$spectra[[1]]), mz(sp2))
    expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
    expect_equal(mz(res$spectra[[2]]), mz(sp1))
    expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
    expect_equal(polarity(res$spectra[[1]]), 1)
    expect_equal(polarity(res$spectra[[2]]), 0)
    expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
    expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))
    
    ## provide rt
    df$rt <- c(1, 1, 1, 1, 3, 3, 3)
    res <- CompoundDb:::.spectra2_from_df(df)
    expect_equal(rtime(res$spectra[[1]]), 1)
    expect_equal(rtime(res$spectra[[2]]), 3)
    expect_equal(mz(res$spectra[[1]]), mz(sp2))
    expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
    expect_equal(mz(res$spectra[[2]]), mz(sp1))
    expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
    expect_equal(polarity(res$spectra[[1]]), 1)
    expect_equal(polarity(res$spectra[[2]]), 0)
    expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
    expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

    ## provide ms_level
    df$ms_level <- c(2, 2, 2, 2, 3, 3, 3)
    res <- CompoundDb:::.spectra2_from_df(df)
    expect_equal(msLevel(res$spectra[[1]]), 2)
    expect_equal(msLevel(res$spectra[[2]]), 3)
    expect_equal(mz(res$spectra[[1]]), mz(sp2))
    expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
    expect_equal(mz(res$spectra[[2]]), mz(sp1))
    expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
    expect_equal(polarity(res$spectra[[1]]), 1)
    expect_equal(polarity(res$spectra[[2]]), 0)
    expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
    expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

    ## provide precursor_mz
    df$precursor_mz <- c(2.1, 2.1, 2.1, 2.1, 3.4, 3.4, 3.4)
    res <- CompoundDb:::.spectra2_from_df(df)
    expect_equal(precursorMz(res$spectra[[1]]), 2.1)
    expect_equal(precursorMz(res$spectra[[2]]), 3.4)
    expect_equal(mz(res$spectra[[1]]), mz(sp2))
    expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
    expect_equal(mz(res$spectra[[2]]), mz(sp1))
    expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
    expect_equal(polarity(res$spectra[[1]]), 1)
    expect_equal(polarity(res$spectra[[2]]), 0)
    expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
    expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

    ## provide precursor_charge
    df$precursor_charge <- c(5, 5, 5, 5, 2, 2, 2)
    res <- CompoundDb:::.spectra2_from_df(df)
    expect_equal(precursorCharge(res$spectra[[1]]), 5)
    expect_equal(precursorCharge(res$spectra[[2]]), 2)
    expect_equal(mz(res$spectra[[1]]), mz(sp2))
    expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
    expect_equal(mz(res$spectra[[2]]), mz(sp1))
    expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
    expect_equal(polarity(res$spectra[[1]]), 1)
    expect_equal(polarity(res$spectra[[2]]), 0)
    expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
    expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

    ## provide precuror_intensity
    df$precursor_intensity <- c(3, 3, 3, 3, 1, 1, 1)
    res <- CompoundDb:::.spectra2_from_df(df)
    expect_equal(precursorIntensity(res$spectra[[1]]), 3)
    expect_equal(precursorIntensity(res$spectra[[2]]), 1)
    expect_equal(mz(res$spectra[[1]]), mz(sp2))
    expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
    expect_equal(mz(res$spectra[[2]]), mz(sp1))
    expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
    expect_equal(polarity(res$spectra[[1]]), 1)
    expect_equal(polarity(res$spectra[[2]]), 0)
    expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
    expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

    ## provide collision_energy
    df$collision_energy <- c(6, 6, 6, 6, 3, 3, 3)
    res <- CompoundDb:::.spectra2_from_df(df)
    expect_equal(collisionEnergy(res$spectra[[1]]), 6)
    expect_equal(collisionEnergy(res$spectra[[2]]), 3)
    expect_equal(mz(res$spectra[[1]]), mz(sp2))
    expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
    expect_equal(mz(res$spectra[[2]]), mz(sp1))
    expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
    expect_equal(polarity(res$spectra[[1]]), 1)
    expect_equal(polarity(res$spectra[[2]]), 0)
    expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
    expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

    ## provide acquisition_num
    df$acquisition_num <- c(1, 1, 1, 1, 2, 2, 2)
    res <- CompoundDb:::.spectra2_from_df(df)
    expect_equal(acquisitionNum(res$spectra[[1]]), 1)
    expect_equal(acquisitionNum(res$spectra[[2]]), 2)
    expect_equal(mz(res$spectra[[1]]), mz(sp2))
    expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
    expect_equal(mz(res$spectra[[2]]), mz(sp1))
    expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
    expect_equal(polarity(res$spectra[[1]]), 1)
    expect_equal(polarity(res$spectra[[2]]), 0)
    expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
    expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

    ## provide scan_index
    df$scan_index <- c(1, 1, 1, 1, 2, 2, 2)
    res <- CompoundDb:::.spectra2_from_df(df)
    expect_equal(scanIndex(res$spectra[[1]]), 1)
    expect_equal(scanIndex(res$spectra[[2]]), 2)
    expect_equal(mz(res$spectra[[1]]), mz(sp2))
    expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
    expect_equal(mz(res$spectra[[2]]), mz(sp1))
    expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
    expect_equal(polarity(res$spectra[[1]]), 1)
    expect_equal(polarity(res$spectra[[2]]), 0)
    expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
    expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

    ## provide from_file
    df$from_file <- c(1, 1, 1, 1, 2, 2, 2)
    res <- CompoundDb:::.spectra2_from_df(df)
    expect_equal(fromFile(res$spectra[[1]]), 1)
    expect_equal(fromFile(res$spectra[[2]]), 2)
    expect_equal(mz(res$spectra[[1]]), mz(sp2))
    expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
    expect_equal(mz(res$spectra[[2]]), mz(sp1))
    expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
    expect_equal(polarity(res$spectra[[1]]), 1)
    expect_equal(polarity(res$spectra[[2]]), 0)
    expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
    expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

    ## provide precursor_scan_num
    df$precursor_scan_num <- c(3, 3, 3, 3, 4, 4, 4)
    res <- CompoundDb:::.spectra2_from_df(df)
    expect_equal(precScanNum(res$spectra[[1]]), 3)
    expect_equal(precScanNum(res$spectra[[2]]), 4)
    expect_equal(mz(res$spectra[[1]]), mz(sp2))
    expect_equal(intensity(res$spectra[[1]]), intensity(sp2))
    expect_equal(mz(res$spectra[[2]]), mz(sp1))
    expect_equal(intensity(res$spectra[[2]]), intensity(sp1))
    expect_equal(polarity(res$spectra[[1]]), 1)
    expect_equal(polarity(res$spectra[[2]]), 0)
    expect_equal(colnames(res$mcols), c("spectrum_id", "compound_id"))
    expect_equal(res$mcols$spectrum_id, unique(df$spectrum_id))

    ## Now the "official" constructor:
    spl <- Spectrum2List(df)
    expect_equal(mz(spl[[1]]), mz(sp2))
    expect_equal(intensity(spl[[1]]), intensity(sp2))
    expect_equal(mz(spl[[2]]), mz(sp1))
    expect_equal(intensity(spl[[2]]), intensity(sp1))

    expect_equal(mcols(spl)$compound_id, c("cp_1", "cp_1"))
    expect_equal(mcols(spl)$spectrum_id, c("b", "a"))
})

test_that(".make_naked_matrix_from_Spectrum2List works", {
    df <- data.frame(spectrum_id = c("b", "b", "b", "b", "a", "a", "a"),
                     mz = c(1, 2, 3, 4, 1, 2, 4),
                     intensity = c(5, 3, 2, 5, 4, 5, 2),
                     polarity = c(1, 1, 1, 1, 0, 0, 0),
                     compound_id = rep("cp_1", 7), stringsAsFactors = FALSE)
    spl <- Spectrum2List(df)
    res <- .make_naked_matrix_from_Spectrum2List(spl)
    expect_equal(ncol(res), 7)
    expect_equal(nrow(res), 2)
    
    spl <- Spectrum2List(new("Spectrum2", mz = c(1, 2, 3), intensity = 1:3))
    res <- .make_naked_matrix_from_Spectrum2List(spl)
    expect_equal(ncol(res), 4)
    expect_equal(nrow(res), 1)
})

test_that("show,Spectrum2List works", {
    df <- data.frame(spectrum_id = c("b", "b", "b", "b", "a", "a", "a"),
                     mz = c(1, 2, 3, 4, 1, 2, 4),
                     intensity = c(5, 3, 2, 5, 4, 5, 2),
                     polarity = c(1, 1, 1, 1, 0, 0, 0),
                     compound_id = rep("cp_1", 7), stringsAsFactors = FALSE)
    spl <- Spectrum2List(df)
    .show_Spectrum2List(spl)
    .show_Spectrum2List(spl, print.classinfo = TRUE)
    show(spl)
})

