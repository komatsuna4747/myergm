#include <RcppArmadillo.h>
#include <RcppSpdlog>
#include <spdlog/stopwatch.h>   		// also support stopwatch feature
#include "MCMC.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat run_network_sampler(
    const arma::mat& adjacency_matrix,
    double coefEdges,
    double coefTwostars,
    double coefTriangle,
    int MCMC_interval = 1024,
    int MCMC_samplesize = 1024,
    int MCMC_burnin = 1024 * 16,
    double p_one_node_swap = 0.001,
    double p_large_step = 0.001,
    double p_invert = 0.001,
    double lambda = 0.5,
    std::string level = "info"
){
  // Necessary for R random number generator
  GetRNGstate();
  arma::mat adjmat = adjacency_matrix;

  std::string logname = "myergm"; 							// fix a name for this logger
  auto sp = spdlog::get(logname);       					// retrieve existing one
  if (sp == nullptr) sp = spdlog::r_sink_mt(logname);   	// or create new one if needed
  spdlog::set_default_logger(sp);                         // and set as default
  spdlog::stopwatch sw;       							// instantiate a stop watch
  spdlog::set_level(spdlog::level::from_str(level));

  // Instantiation
  NetworkSampler network_sampler(
      adjmat,
      coefEdges,
      coefTwostars,
      coefTriangle,
      MCMC_interval,
      MCMC_samplesize,
      MCMC_burnin,
      p_one_node_swap,
      p_large_step,
      p_invert,
      lambda
  );

  // Sample networks and store their sufficient statistics
  network_sampler.run_simulation();

  // This must be called after GetRNGstate before returning to R.

  PutRNGstate();
  return network_sampler.netstats;
}


// [[Rcpp::export]]
Rcpp::List simulate_network(
    const arma::mat& adjacency_matrix,
    double coefEdges,
    double coefTwostars,
    double coefTriangle,
    int MCMC_interval = 1024,
    int MCMC_samplesize = 1024,
    int MCMC_burnin = 1024 * 16,
    double p_one_node_swap = 0.001,
    double p_large_step = 0.001,
    double p_invert = 0.001,
    double lambda = 0.5,
    std::string level = "info"
){
  // Necessary for R random number generator
  GetRNGstate();
  arma::mat adjmat = adjacency_matrix;

  std::string logname = "myergm"; 							// fix a name for this logger
  auto sp = spdlog::get(logname);       					// retrieve existing one
  if (sp == nullptr) sp = spdlog::r_sink_mt(logname);   	// or create new one if needed
  spdlog::set_default_logger(sp);                         // and set as default
  spdlog::stopwatch sw;       							// instantiate a stop watch
  spdlog::set_level(spdlog::level::from_str(level));

  // Instantiation
  NetworkSampler network_sampler(
      adjmat,
      coefEdges,
      coefTwostars,
      coefTriangle,
      MCMC_interval,
      MCMC_samplesize,
      MCMC_burnin,
      p_one_node_swap,
      p_large_step,
      p_invert,
      lambda
  );

  // Sample networks and store their sufficient statistics
  network_sampler.run_simulation();

  // This must be called after GetRNGstate before returning to R.

  PutRNGstate();
  Rcpp::List output(2);
  output[0] = network_sampler.adjmat;
  output[1] = network_sampler.netstats;

  return output;
}
