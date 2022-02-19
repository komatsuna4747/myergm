#include <RcppArmadillo.h>
#include "MCMC.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat run_network_sampler(
    const arma::mat& adjacency_matrix,
    double coefEdges,
    double coefTriangle,
    int MCMC_interval = 1024,
    int MCMC_samplesize = 1024,
    int MCMC_burnin = 1024 * 16,
    double p_one_node_swap = 0.01,
    double p_large_step = 0.01,
    double p_invert = 0.01,
    double lambda = 0.5,
    int verbose = 0
){
  // Necessary for R random number generator
  GetRNGstate();
  arma::mat adjmat = adjacency_matrix;

  // Instantiation
  NetworkSampler network_sampler(
      adjmat,
      coefEdges,
      coefTriangle,
      MCMC_interval,
      MCMC_samplesize,
      MCMC_burnin,
      p_one_node_swap,
      p_large_step,
      p_invert,
      lambda,
      verbose
  );

  // Sample networks and store their sufficient statistics
  network_sampler.run_simulation();

  // This must be called after GetRNGstate before returning to R.

  PutRNGstate();
  return network_sampler.netstats;
}
