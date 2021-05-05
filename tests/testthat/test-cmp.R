test_that("Expectation of CMP", {
	expect_equal(ecmp(lambda = 10, nu = 1), 10, tolerance = 1e-4)
	expect_equal(ecmp(lambda = 0.5, nu = 1), 0.5, tolerance = 1e-4)
	expect_equal(ecmp(lambda = 0.5, nu = 100), 0.5 / (1 + 0.5), tolerance = 1e-2)
	expect_equal(ecmp(lambda = 0.25, nu = 0.0001), 0.25 / (1 - 0.25), tolerance = 1e-2)
})
