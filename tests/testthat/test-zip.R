test_that("Density of ZIP", {
	expect_equal(dzip(x = 0, lambda = 10, p = 0, log = TRUE), dpois(x = 0, lambda = 10, log = TRUE))
	expect_equal(dzip(x = 0, lambda = 10, p = 0), dpois(x = 0, lambda = 10))
	expect_equal(dzip(x = 1, lambda = 5, p = 0), dpois(x = 1, lambda = 5))
	expect_equal(dzip(x = 5, lambda = 5, p = 0.25), 0.75*dpois(x = 5, lambda = 5))
	expect_equal(dzip(x = 0, lambda = 5, p = 0.25), 0.75*dpois(x = 0, lambda = 5) + 0.25)
})
