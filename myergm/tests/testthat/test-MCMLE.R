library(ergm)
data(florentine)

test_that("Stopping the process when a wrong formula is provided works", {
  model <- flomarriage ~ edges + triangle + kstar(2)
  expect_error(myergm_MCMLE(model = model))

  model <- flomarriage ~ edges + kstar(2) + triangle
  expect_error(myergm_MCMLE(model = model))
})

test_that("Yielding identical estimates from the same seed value works", {
  theta1 <- myergm_MCMLE(model = flomarriage ~ edges + triangle,
                         seed = 334,
                         MCMC_interval = 100,
                         MCMC_samplesize = 100,
                         MCMC_burnin = 100)
  theta2 <- myergm_MCMLE(model = flomarriage ~ edges + triangle,
                         seed = 334,
                         MCMC_interval = 100,
                         MCMC_samplesize = 100,
                         MCMC_burnin = 100)
  expect_equal(theta1, theta2, check.attributes = FALSE, tolerance = 1e-16)
})
