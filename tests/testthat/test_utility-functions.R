test_that("expandMzIntensity, .expand and .collapse_spectrum_df work", {
    tmp <- msms_spctra                  # collapsed one.
    res <- .expand_spectrum_df(tmp)
    expect_equal(nrow(res), length(unlist(tmp$mz)))
    expect_equal(unique(res$spectrum_id), tmp$spectrum_id)
    expect_equal(res$mz, unlist(tmp$mz))
    expect_equal(res$intensity, unlist(tmp$intensity))

    ## Collapsing
    res_2 <- .collapse_spectrum_df(res)
    expect_equal(res_2, tmp)

    res_3 <- expandMzIntensity(tmp)
    expect_equal(res_3, res)
    expect_error(expandMzIntensity(4))
})

test_that(".extract_field_from_string works", {
    strng <- "some=bl df;other=some nice thing;last=the last entry"
    res <- .extract_field_from_string(strng, "some=", ";")
    expect_equal(res, "bl df")
    res <- .extract_field_from_string(strng, "notthere=", ";")
    expect_equal(res, NA_character_)
    res <- .extract_field_from_string(strng, "last=", ";")
    expect_equal(res, "the last entry")

    strngs <- c(strng, "some=second")
    res <- CompoundDb:::.extract_field_from_string(strngs, "some=", ";")
    expect_equal(res, c("bl df", "second"))
    res <- CompoundDb:::.extract_field_from_string(strngs, "last=", ";")
    expect_equal(res, c("the last entry", NA))
})
