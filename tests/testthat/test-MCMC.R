library(ergm)
data(florentine)
g <- intergraph::asIgraph(flomarriage)
g <- igraph::as_adjacency_matrix(g)
true_stats <- summary(flomarriage ~ edges + triangle)


test_that("counting edges of a network works", {
  expect_equal(count_edges_cpp(g), true_stats[[1]])
})


test_that("counting triangles of a network works", {
  expect_equal(count_triangle_cpp(g), true_stats[[2]])
})


test_that("changing one link of a network works", {
  g_next <- change_one_link(g, count_edges_cpp(g), count_triangle_cpp(g), i = 0, j = 2, verbose = 0)
  expect_false(g[1, 3] == g_next[1, 3])
  expect_false(g[3, 1] == g_next[3, 1])
})

test_that("changing one link of a network works", {
  g_next <- change_one_link(g, count_edges_cpp(g), count_triangle_cpp(g), i = 0, j = 2, verbose = 0)
  expect_false(g[1, 3] == g_next[1, 3])
  expect_false(g[3, 1] == g_next[3, 1])
})


test_that("changing all link of one vertex works", {
  i <- 5
  g_next <- change_all_links_of_one_node(g, i-1)
  index <- 1:nrow(g)
  index <- index[index != i]
  # Check if all the entries except the diagonal are different.
  expect_true(all(g_next[i, index] != g[i, index]))
  expect_true(all(g_next[index, i] != g[index, i]))
  # Check if the diagonal is still zero.
  expect_equal(g[i, i], 0)
})


test_that("changing multiple links in one step works", {
  g_next <- change_multiple_links(g, lambda = 0.5)
  # Check if all the entries except the diagonal are different.
  expect_true(any(g != g_next))
  # Check if the diagonals are zero.
  expect_true(all(Matrix::diag(g_next) == 0))
})


test_that("Creating identical MCMC networks from the same seed value works", {
  set.seed(334)
  mcmc1 <- create_MCMC(adjmat = g,
                       coefEdges = -1.8,
                       coefTriangle = 0.2,
                       MCMC_interval = 100,
                       MCMC_samplesize = 100,
                       MCMC_burnin = 100,
                       full_sample = TRUE)

  set.seed(334)
  mcmc2 <- create_MCMC(adjmat = g,
                       coefEdges = -1.8,
                       coefTriangle = 0.2,
                       MCMC_interval = 100,
                       MCMC_samplesize = 100,
                       MCMC_burnin = 100,
                       full_sample = TRUE)
  expect_true(all(mcmc1 == mcmc2))
})

