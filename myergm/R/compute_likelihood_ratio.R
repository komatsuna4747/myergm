#' Compute the loglikelihood ratio
compute_loglikelihood_ratio <- function(theta, theta0, g_obs, g_mcmc) {
  # Compute the first term
  theta <- matrix(theta, nrow = 2, ncol = 1)
  first_term <- sum(g_obs %*% (theta - theta0))
  # Compute the second term
  M <- nrow(g_mcmc)
  second_term <- exp(g_mcmc %*% (theta - theta0))
  second_term <- sum(second_term) / M
  second_term <- log(second_term)
  # Compute the likelihood ratio
  loglikelihood_ratio <- first_term - second_term
  # Return the output
  return(loglikelihood_ratio)
}
