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

test_that(".column_indices works", {
    df <- data.frame(A = 1:4, B = 1:4 , C = 1:4, D = 1:4)
    expect_equal(.column_indices(df, c("C", "B")), c(3, 2))
    expect_error(.column_indices(df, c("C", "B", "Z")))
    expect_error(.column_indices(df, 1:5))
    expect_error(.column_indices(df, 6L))
    expect_equal(.column_indices(df, c(2, 3)), c(2, 3))
    expect_true(length(.column_indices(df, integer())) == 0)
    expect_error(.column_indices(df, c(TRUE, FALSE)))
    expect_equal(.column_indices(df, c(TRUE, FALSE, TRUE, FALSE)), c(1L, 3L))
    expect_equal(.column_indices(df, rep(FALSE, 4)), integer())
})

test_that(".collapse_table works", {
    df <- data.frame(A = c("a", "a", "a", "c", "c", "d", "d", "d"),
                     B = c("e", "e", "f", "f", "g", "g", "h", "h"),
                     C = c("m", "n", "o", "p", "q", "r", "s", "t"),
                     D = c(1, 2, 1, 2, 1, 1, 2, 2),
                     stringsAsFactors = FALSE)
    expect_error(.collapse_table(df, by = "G"))
    res <- .collapse_table(df, by = 1)
    expect_equal(nrow(res), length(unique(df$A)))
    expect_equal(res$A, unique(df$A))
    expect_true(is.list(res$B))
    expect_true(is.list(res$C))
    expect_true(is.list(res$D))
    expect_equal(res$B, list(c("e", "f"), c("f", "g"), c("g", "h")))
    expect_equal(res$D, list(c(1, 2), c(2, 1), c(1, 2)))

    res <- .collapse_table(df, by = "D")
    expect_equal(res$D, c(1, 2))
    expect_equal(res$A, list(c("a", "c", "d"), c("a", "c", "d")))

    res <- .collapse_table(df, by = 1:2)
    res_unique <- unique(df[, 1:2])
    rownames(res_unique) <- NULL
    expect_equal(res[, 1:2], res_unique)
    expect_equal(res$C, list(c("m", "n"), c("o"), c("p"), c("q"), c("r"),
                             c("s", "t")))
    expect_equal(res$D, list(c(1, 2), c(1), c(2), c(1), c(1), c(2)))
})


## test_that(".inchikey2id works", {
##     ids <- c("b", "a", "a", "c", "d", "a", "f", "g", "b", "d", "a")
##     ids_new <- .inchikey2id(ids)
##     expect_equal(length(ids_new), length(ids))
##     expect_equal(length(unique(ids_new)), length(unique(ids)))
##     expect_equal(as.numeric(factor(ids_new)), as.numeric(factor(ids)))
## })
