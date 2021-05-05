#' Implement an MCMLE to estimate parameters.
#' @param model a model to be estimated. It only supports `network ~ edges + triangle`.
#' @param seed Seed value (integer) for the random number generator.
#' @param MCMC_interval Number of proposals between sampled statistics.
#' @param MCMC_samplesize Number of network statistics, randomly drawn from a given distribution on the set of all networks,
#' returned by the Metropolis-Hastings algorithm.
#' @param MCMC_burnin Number of proposals before any MCMC sampling is done.
#' @param p_one_node_swap Probabiluty of changing all links of aone node.
#' @param p_large_step Probability of swapping multiple links in one step.
#' The number of links to be swapped is determined by lambda * N.
#' @param p_invert Probability of inverting the adjacency matrix.
#' When the number of vertices is large, it is recommended to set this probaility to be zero.
#' @param lambda Fraction of links updated for large step.
#' @param verbose A logical or an integer: if this is TRUE/1, the program will print out additional information about the progress of estimation and simulation.
#' Higher values produce more verbosity.
#' @param second_round A logical: If TRUE, MCMLE is implemented twice. Otherwise, only once.
#' @param ... Additional arguments, passed to other functions.
#' @import ergm
#' @importFrom ergm ergm.getnetwork
#' @importFrom intergraph asIgraph
#' @importFrom igraph as_adjacency_matrix
#' @importFrom glue glue
#' @importFrom stats cov as.formula
#' @importFrom stringr str_split
#' @return Estimates of `edges` and `triangle` parameters
#' @examples
#' # Load a network object
#' library(ergm)
#' data(florentine)
#'
#' # Estimate the parameters
#' theta <- myergm_MCMLE(model = flomarriage ~ edges + triangle, seed = 334)
#' @export
myergm_MCMLE <- function(model,
                         seed = NULL,
                         MCMC_interval = 1024,
                         MCMC_samplesize = 1024,
                         MCMC_burnin = 1024 * 16,
                         p_one_node_swap = 0.01,
                         p_large_step = 0.01,
                         p_invert = 0.001,
                         lambda = 0.5,
                         verbose = 0,
                         second_round = TRUE,
                         ...) {
  # Check if the formula is valid.
  rhs <- as.character(as.formula(model))[[3]]
  terms_rhs <- unlist(stringr::str_split(string = rhs, pattern = " \\+ "))
  if (length(terms_rhs) != 2 | terms_rhs[[1]] != "edges" | terms_rhs[[2]] != "triangle") {
    stop("myergm_MCMLE() only accepts '~ edges + triangle' model.")
  }

  # Set seed
  set.seed(seed)

  # Estimate theta0 by MPLE. I am a lazy man, so here I use the `ergm` library to get theta0.
  if (verbose >= 1) {
    message("Starting MPLE...")
  }
  MPLE <- ergm(model, estimate = "MPLE")
  theta0 <- MPLE$coef

  # Get the adjacency matrix
  if (verbose >= 1) {
    message("Converting the network object into an adjacency matrix...")
  }
  network <- ergm::ergm.getnetwork(model)
  g <- intergraph::asIgraph(network)
  g <- igraph::as_adjacency_matrix(g)

  # Get observed network statistics
  if (verbose >= 1) {
    message("Getting observed network statistics...")
  }
  g_obs <- matrix(c(count_edges_cpp(g), count_triangle_cpp(g)), nrow = 1, ncol = 2)

  # Get network statistics, randomly drawn from a given distribution on set of all networks, returned by the Metropolis-Hastings algorithm.
  if (verbose >= 1) {
    message("Starting the first round of MCMC...")
  }
  mcmc <- create_MCMC(adjmat = g,
                      coefEdges = theta0[[1]],
                      coefTriangle = theta0[[2]],
                      MCMC_samplesize = MCMC_samplesize,
                      MCMC_burnin = MCMC_burnin,
                      MCMC_interval = MCMC_interval,
                      p_one_node_swap = p_one_node_swap,
                      p_large_step = p_large_step,
                      p_invert = p_invert,
                      lambda = lambda,
                      full_sample = FALSE,
                      verbose = verbose)

  # Get the estimates
  if (verbose >= 1) {
    message("Using log-normal approximation (no optim)")
  }
  m0 <- matrix(c(mean(mcmc[,1]), mean(mcmc[,2])), nrow = 1, ncol = 2)
  theta <- theta0 + solve(cov(mcmc)) %*% t(g_obs - m0)

  # Print the estimates
  message(glue::glue("Parameter estimates in the first round\n edges: {theta[[1]]}, triangle: {theta[[2]]}"))

  # Start the second round of MCMC if second_round = TRUE.
  if (second_round == TRUE) {
    if (verbose >= 1) {
      message(glue::glue("Starting the second round of MCMC...\n Increasing MCMC sample size from {MCMC_samplesize} to {4 * MCMC_samplesize}."))
    }
    mcmc <- create_MCMC(adjmat = g,
                        coefEdges = theta[[1]],
                        coefTriangle = theta[[2]],
                        MCMC_samplesize = MCMC_samplesize * 4,
                        MCMC_burnin = MCMC_burnin,
                        MCMC_interval = MCMC_interval,
                        p_one_node_swap = p_one_node_swap,
                        p_large_step = p_large_step,
                        p_invert = p_invert,
                        lambda = lambda,
                        full_sample = FALSE,
                        verbose = verbose)

    # Get the estimates
    if (verbose >= 1) {
      message("Using log-normal approximation (no optim)")
    }
    m0 <- matrix(c(mean(mcmc[,1]), mean(mcmc[,2])), nrow = 1, ncol = 2)
    theta <- theta + solve(cov(mcmc)) %*% t(g_obs - m0)

    # Print the estimates if verbose >= 1
    message(glue::glue("Parameter estimates in the second round\n edges: {theta[[1]]}, triangle: {theta[[2]]}"))
  }

  # Return the output
  rownames(theta) <- c("edges", "triangle")
  return(theta)
}
