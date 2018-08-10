test_that(".expand and .collapse_spectrum_df work", {
    tmp <- msms_spctra                  # collapsed one.
    res <- CompoundDb:::.expand_spectrum_df(tmp)
    expect_equal(nrow(res), length(unlist(tmp$mz)))
    expect_equal(unique(res$spectrum_id), tmp$spectrum_id)
    expect_equal(res$mz, unlist(tmp$mz))
    expect_equal(res$intensity, unlist(tmp$intensity))

    ## Collapsing
    res_2 <- CompoundDb:::.collapse_spectrum_df(res)
    expect_equal(res_2, tmp)
})
