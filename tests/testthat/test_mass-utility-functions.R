test_that("matchWithPpm works", {
    set.seed(123)
    yvals <- abs(rnorm(1000))
    yvals[123] <- yvals[2] - yvals[2] * 10 * 1e-6
    yvals[124] <- yvals[2] + yvals[2] * 10 * 1e-6
    yvals[125] <- yvals[2] + yvals[2] * 12 * 1e-6
    xvals <- yvals[c(2, 3, 3, 20, 21, 20)]
    xvals[2] <- xvals[2] + (10 * xvals[2] / 1e6)
    xvals[3] <- xvals[3] - (10 * xvals[3] / 1e6)
    xvals[6] <- xvals[6] + (12 * xvals[6] / 1e6)

    res <- matchWithPpm(xvals, yvals)
    expect_true(is.list(res))
    expect_equal(length(res), length(xvals))
    expect_equal(res[[1]], 2L)
    expect_equal(res[[2]], integer())
    expect_equal(res[[3]], integer())
    expect_equal(res[[4]], 20L)


    res <- matchWithPpm(xvals, yvals, ppm = 10)
    expect_equal(res[[1]], c(2L, 123L, 124L))
    expect_equal(res[[2]], 3L)
    expect_equal(res[[3]], 3L)

    res <- matchWithPpm(xvals, yvals, ppm = 20)
    expect_equal(res[[1]], c(2L, 123L, 124L, 125L))
})

test_that("adducts works", {
    all <- CompoundDb:::ADDUCTS
    res <- adducts()
    expect_equal(all, res)
    res <- adducts(polarity = -4)
    expect_equal(res, all[all$charge < 0, ])
    res <- adducts(pattern = "H+")
    expect_equal(res, all[grep("H+", all$name), ])
    res <- adducts(pattern = "H+", fixed = TRUE)
    expect_equal(res, all[grep("H+", all$name, fixed = TRUE), ])
    res <- adducts(name = "[M+H]+")
    expect_true(nrow(res) == 1)
})

test_that("mass2mz works", {
    library(testthat)
    library(CompoundDb)
    masses <- c(75.032028409, 105.042595)
    names(masses) <- c("Glycine", "Serine")
    ## [M+H]+
    res <- mass2mz(masses, adduct = "[M+H]+")
    expect_true(is.list(res))
    expect_equal(names(res), names(masses))
    expect_true(all(lengths(res) == 1))
    expect_equal(unname(res[[1]]), 76.039304409)
    expect_equal(unname(res[[2]]), 106.049871)

    ## more complicated stuff - cross checked with HMDB
    res <- mass2mz(masses, adduct = "[2M+K]+")
    expect_true(all.equal(unname(res[[1]]), 189.0272, tolerance = 1e-5))
    res <- mass2mz(masses, adduct = "[M+2H]2+")
    expect_true(all.equal(unname(res[[1]]), 38.5233, tolerance = 1e-5))

    ## [M-H]-
    expct <- c(74.024752409, 104.035319)
    res <- mass2mz(masses, adduct = "[M-H]-")
    expect_equal(unname(res[[1]]), 74.024752409)
    expect_equal(unname(res[[2]]), 104.035319)

    ## All of them
    res <- mass2mz(masses)
})

test_that("mz2mass works", {
    mass <- 75.032028409

    ## m/z taken from HMDB
    res <- mz2mass(76.0393, adduct = "[M+H]+")
    expect_equal(names(res[[1]]), "[M+H]+")
    expect_true(all.equal(unname(res[[1]]), mass, tolerance = 1e-6))
    res <- mz2mass(120.0032, adduct = "[M+2Na-H]+")
    expect_true(all.equal(unname(res[[1]]), mass, tolerance = 1e-6))
    res <- mz2mass(151.0713, adduct = "[2M+H]+")
    expect_true(all.equal(unname(res[[1]]), mass, tolerance = 1e-6))
    res <- mz2mass(38.5233, adduct = "[M+2H]2+")
    expect_true(all.equal(unname(res[[1]]), mass, tolerance = 1e-6))

    tmp <- mass2mz(mass, adduct = "[2M+2H]2+")
    res <- mz2mass(tmp, adduct = "[2M+2H]2+")
    expect_equal(unname(res[[1]]), mass)
})
