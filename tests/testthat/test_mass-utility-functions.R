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
