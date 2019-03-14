test_that(".massbank_extract_field works", {
    str <- c("CH$NAME: A", "CH$OTHER: B", "CH$NAME: C")
    res <- .massbank_extract_field(field = "CH$NAME: ", str)
    expect_equal(res, c("A", "C"))
    res <- .massbank_extract_field(field = "CH$OTHER: ", str)
    expect_equal(res, "B")
    res <- .massbank_extract_field(field = "nothing", str)
    expect_equal(res, NA_character_)
    res <- .massbank_extract_field(x = character())
    expect_equal(res, NA_character_)
})

test_that(".import_massbank_file works", {
    fls <- dir(system.file("txt", package = "CompoundDb"), pattern = ".txt$",
               full.names = TRUE)
    res <- CompoundDb:::.import_massbank_file(fls[1])
    expect_equal(nrow(res), 1)
    expect_true(is.character(res$inchi))
    expect_true(is.numeric(res$mass))

    fl <- system.file("sdf/ChEBI_sub.sdf.gz")
    res2 <- CompoundDb:::.import_massbank_file(fl)
    expect_true(nrow(res2) == 0)
    expect_equal(colnames(res), colnames(res2))
})
