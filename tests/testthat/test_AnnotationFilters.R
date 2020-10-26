test_that("CompoundIdFilter, .field, .sql_condition, sql_value work", {
    fl <- CompoundIdFilter("samid")
    expect_true(is(fl, "CompoundIdFilter"))
    expect_true(is(fl, "CharacterFilter"))
    expect_true(is(fl, "AnnotationFilter"))

    expect_error(CompoundIdFilter())

    expect_equal(.field(fl), "compound_id")
    expect_equal(.sql_condition(fl), "=")
    expect_equal(.sql_value(fl), "'samid'")
})

test_that("NameFilter works", {
    fl <- NameFilter("a")
    expect_true(is(fl, "NameFilter"))
    expect_true(is(fl, "CharacterFilter"))
    expect_true(is(fl, "AnnotationFilter"))

    expect_error(NameFilter())

    expect_equal(.field(fl), "name")
    expect_equal(.sql_condition(fl), "=")
    expect_equal(.sql_value(fl), "'a'")
})

test_that("SpectrumIdFilter works", {
    fl <- SpectrumIdFilter("a")
    expect_true(is(fl, "SpectrumIdFilter"))
    expect_true(is(fl, "CharacterFilter"))
    expect_true(is(fl, "AnnotationFilter"))

    expect_error(SpectrumIdFilter())

    expect_equal(.field(fl), "spectrum_id")
    expect_equal(.sql_condition(fl), "=")
    expect_equal(.sql_value(fl), "'a'")
})

test_that("FormulaFilter works", {
    fl <- FormulaFilter("a")
    expect_true(is(fl, "FormulaFilter"))
    expect_true(is(fl, "CharacterFilter"))
    expect_true(is(fl, "AnnotationFilter"))

    expect_error(FormulaFilter())

    expect_equal(.field(fl), "formula")
    expect_equal(.sql_condition(fl), "=")
    expect_equal(.sql_value(fl), "'a'")
})

test_that("InchiFilter works", {
    fl <- InchiFilter("a")
    expect_true(is(fl, "InchiFilter"))
    expect_true(is(fl, "CharacterFilter"))
    expect_true(is(fl, "AnnotationFilter"))

    expect_error(InchiFilter())

    expect_equal(.field(fl), "inchi")
    expect_equal(.sql_condition(fl), "=")
    expect_equal(.sql_value(fl), "'a'")
})

test_that("InchikeyFilter works", {
    fl <- InchikeyFilter("a")
    expect_true(is(fl, "InchikeyFilter"))
    expect_true(is(fl, "CharacterFilter"))
    expect_true(is(fl, "AnnotationFilter"))

    expect_error(InchikeyFilter())

    expect_equal(.field(fl), "inchikey")
    expect_equal(.sql_condition(fl), "=")
    expect_equal(.sql_value(fl), "'a'")
})

test_that("ExactmassFilter works", {
    fl <- ExactmassFilter(12.2)
    expect_true(is(fl, "ExactmassFilter"))
    expect_true(is(fl, "DoubleFilter"))
    expect_true(is(fl, "AnnotationFilter"))

    expect_error(ExactmassFilter())

    expect_equal(.field(fl), "exactmass")
    expect_equal(.sql_condition(fl), "=")
    expect_equal(.sql_value(fl), 12.2)
})

test_that(".field works", {
    library(AnnotationFilter)
    gif <- GeneIdFilter("a")
    sf <- SymbolFilter("b")
    tif <- TxIdFilter("c")
    expect_equal(.field(gif), "gene_id")
    expect_equal(.field(AnnotationFilterList(gif)), "gene_id")
    expect_equal(.field(AnnotationFilterList(gif, tif)), c("gene_id", "tx_id"))
    expect_equal(.field(
        AnnotationFilterList(gif, AnnotationFilterList(tif, sf))),
        c("gene_id", "tx_id", "symbol"))
})

test_that(".process_filter works", {
    library(AnnotationFilter)
    gif <- GeneIdFilter("a")
    fl <- CompoundIdFilter("d")

    expect_error(.process_filter("3"))
    expect_error(.process_filter(gif))
    expect_error(.process_filter(AnnotationFilterList(gif, fl)))

    expect_equal(.process_filter(fl), AnnotationFilterList(fl))
    expect_equal(.process_filter(~compound_id == "d"), AnnotationFilterList(fl))

    expect_error(.process_filter(~ compound_id == "d" & msms_mz_range_min > 12))
    expect_error(.process_filter(~ compound_id == "d" & msms_mz_range_min > 12,
                                 cmp_db))
    res <- .process_filter(~ compound_id == "d" & msms_mz_range_min > 12,
                           cmp_spctra_db)
    expect_true(is(res, "AnnotationFilterList"))
})

test_that(".sql_condition works", {
    fl <- CompoundIdFilter("a")
    expect_equal(.sql_condition(fl), "=")
    fl <- CompoundIdFilter("a", "!=")
    expect_equal(.sql_condition(fl), "!=")
    fl <- CompoundIdFilter(c("a", "b"), "!=")
    expect_equal(.sql_condition(fl), "not in")
    fl <- CompoundIdFilter(c("a", "b"), "==")
    expect_equal(.sql_condition(fl), "in")
    fl <- CompoundIdFilter("a", "startsWith")
    expect_equal(.sql_condition(fl), "like")
})

test_that(".sql_value works", {
    fl <- CompoundIdFilter("a")
    expect_equal(.sql_value(fl), "'a'")
    fl <- CompoundIdFilter(c("a", "b"))
    expect_equal(.sql_value(fl), "('a','b')")
    fl <- CompoundIdFilter("a", condition = "startsWith")
    expect_equal(.sql_value(fl), "'a%'")
    fl <- CompoundIdFilter("a", condition = "endsWith")
    expect_equal(.sql_value(fl), "'%a'")
    fl <- CompoundIdFilter("a", condition = "contains")
    expect_equal(.sql_value(fl), "'%a%'")
})

test_that(".sql_logicOp works", {
    afl <- AnnotationFilter(~ compound_id == "a" & name == "2323434")
    expect_equal(.sql_logicOp(afl), "and")
    afl <- AnnotationFilter(~ compound_id == "a" | name == "2323434")
    expect_equal(.sql_logicOp(afl), "or")
    afl <- AnnotationFilter(~ compound_id == "a" & name == "2323434" |
                            gene_id == "123")
    expect_equal(.sql_logicOp(afl), c("and", "or"))
})

test_that(".where_filter works", {
    fl <- CompoundIdFilter("5")
    afl <- AnnotationFilter(~ compound_id == "a" & name == "1")
    expect_equal(.where_filter(fl), "compound_id = '5'")
    expect_equal(.where_filter(afl),
                 "(compound_id = 'a' and name = '1')")
    afl_2 <- AnnotationFilterList(fl, afl, logicOp = "|")
    expect_equal(.where_filter(afl_2),
                 paste0("(compound_id = '5' or (compound_id =",
                        " 'a' and name = '1'))"))
    afl_2 <- AnnotationFilterList(afl_2, afl, logicOp = "&")
    expect_equal(.where_filter(afl_2),
                 paste0("((compound_id = '5' or (compound_id",
                        " = 'a' and name = '1')",
                        ") and (compound_id = 'a' and ",
                        "name = '1'))"))
    res <- .where_filter(fl, c(compound_id = "test.compound_id"))
    expect_equal(res, "test.compound_id = '5'")
})

test_that("MsmsMzRangeMinFilter works", {
    fl <- MsmsMzRangeMinFilter(12.3)
    expect_true(is(fl, "MsmsMzRangeMinFilter"))
    expect_true(is(fl, "DoubleFilter"))

    expect_equal(value(fl), 12.3)

    expect_error(MsmsMzRangeMinFilter())
    expect_error(MsmsMzRangeMinFilter("d"))
    expect_error(MsmsMzRangeMinFilter(c(2.4, 3.5)))

    expect_equal(.field(fl), "msms_mz_range_min")
    expect_equal(.sql_condition(fl), ">=")
    expect_equal(.sql_value(fl), 12.3)
})

test_that("MsmsMzRangeMaxFilter works", {
    fl <- MsmsMzRangeMaxFilter(12.3)
    expect_true(is(fl, "MsmsMzRangeMaxFilter"))
    expect_true(is(fl, "DoubleFilter"))

    expect_equal(value(fl), 12.3)

    expect_error(MsmsMzRangeMaxFilter())
    expect_error(MsmsMzRangeMaxFilter("d"))
    expect_error(MsmsMzRangeMaxFilter(c(2.4, 3.5)))

    expect_equal(.field(fl), "msms_mz_range_max")
    expect_equal(.sql_condition(fl), "<=")
    expect_equal(.sql_value(fl), 12.3)
})

test_that(".supported_filters works", {
    res <- .supported_filters()
    expect_equal(colnames(res), c("filter", "field"))
    res_2 <- .supported_filters(cmp_db)
    expect_equal(res, res_2)
    res_2 <- .supported_filters(cmp_spctra_db)
    expect_true(nrow(res_2) > nrow(res))
})

test_that(".filter_class works", {
    res <- .filter_class(MsmsMzRangeMinFilter(3))
    expect_equal(res, "MsmsMzRangeMinFilter")
    res <- .filter_class(AnnotationFilterList(CompoundIdFilter("a"),
                                              MsmsMzRangeMinFilter(3)))
    expect_equal(res, c("CompoundIdFilter", "MsmsMzRangeMinFilter"))
})
