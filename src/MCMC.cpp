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

  std::string logname = "run_network_sampler"; 		// fix a name for this logger
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

//' Simulate networks
//' @export
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

  std::string logname = "simulate_network"; 			// fix a name for this logger
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


// https://gallery.rcpp.org/articles/simulate-multivariate-normal/
// [[Rcpp::export]]
arma::vec proposal_mvrnorm(arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  arma::mat output = arma::repmat(mu, 1, 1).t() + Y * arma::chol(sigma);
  output = output.t();
  return output.col(0);
}

//' Run a Double Metropolis-Hastings sampler
//' @export
// [[Rcpp::export]]
arma::mat run_DMH_cpp(
    const arma::mat& adjacency_matrix,
    arma::mat sigma_proposal,
    arma::vec theta_init,
    arma::mat sigma_prior,
    arma::vec terms_included,
    int param_sample = 10000,
    int MCMC_network_interval = 100000,
    double p_one_node_swap = 0.001,
    double p_large_step = 0.001,
    double p_invert = 0.001,
    double lambda = 0.5,
    std::string level = "info"
) {
  // Necessary for R random number generator
  GetRNGstate();
  arma::mat adjmat = adjacency_matrix;

  std::string logname = "run_DMH"; 		// fix a name for this logger
  auto sp = spdlog::get(logname);       					// retrieve existing one
  if (sp == nullptr) sp = spdlog::r_sink_mt(logname);   	// or create new one if needed
  spdlog::set_default_logger(sp);                         // and set as default
  spdlog::stopwatch sw;       							// instantiate a stop watch
  spdlog::set_level(spdlog::level::from_str(level));

  // Instantiation
  NetworkSampler network_sampler(
      adjmat,
      theta_init[0],
      theta_init[1],
      theta_init[2],
      MCMC_network_interval,
      1, // MCMC_samplesize
      1, // MCMC_burnin
      p_one_node_swap,
      p_large_step,
      p_invert,
      lambda
  );

  arma::mat sampled_params = arma::zeros(param_sample, 3);
  arma::vec proposed_theta(3);
  arma::vec theta = theta_init % terms_included;
  arma::mat inv_sigma_prior = arma::inv(sigma_prior);

  arma::mat observed_stat = arma::zeros(1, 3);
  observed_stat(0, 0) = network_sampler.numOfEdges;
  observed_stat(0, 1) = network_sampler.numOfTwostars;
  observed_stat(0, 2) = network_sampler.numOfTriangles;

  arma::mat sampled_stat = arma::zeros(1, 3);
  double metropolis_ratio = 0;
  double acceptance_prob = 0;
  double random_number = 0;
  double n_param_accept = 0;

  for (int i = 0; i < param_sample; i++) {
    if (i % 100 == 0) {
      spdlog::info(
        "Parameter sampling interation: {0} / {1}",
        i, param_sample
      );
    }

    spdlog::trace(
      "theta_edges: {0}, theta_twostar: {1}, theta_triangle: {2}",
      theta[0], theta[1], theta[2]
    );

    // Propose a parameter vector
    proposed_theta = proposal_mvrnorm(theta, sigma_proposal) % terms_included;
    spdlog::trace(
      "p_edges: {0}, p_twostar: {1}, p_triangle: {2}",
      proposed_theta[0], proposed_theta[1], proposed_theta[2]
    );

    network_sampler.coefEdges = proposed_theta[0];
    network_sampler.coefTwostars = proposed_theta[1];
    network_sampler.coefTriangle = proposed_theta[2];

    network_sampler.run_simulation();
    sampled_stat = network_sampler.netstats;

    metropolis_ratio = dot(theta - proposed_theta, sampled_stat.row(0) - observed_stat.row(0));
    metropolis_ratio += - 0.5 * sum((proposed_theta - theta).t() * inv_sigma_prior * (proposed_theta - theta));

    acceptance_prob = metropolis_ratio >= 0 ? 0 : metropolis_ratio;
    random_number = log(unif_rand());
    spdlog::trace(
      "Metropolis ratio: {0}, Random number: {1}",
      acceptance_prob, random_number
    );

    // Accept or reject?
    if (random_number < acceptance_prob) {
      spdlog::trace("Accepted.");
      theta = proposed_theta;
      n_param_accept += 1;
    } else {
      spdlog::trace("Rejected.");
    }

    sampled_params.row(i) = theta.t();

    // Reset internal parameters
    network_sampler.reset_parameters(adjmat);
  }

  spdlog::debug(
    "Acceptance rate of parameter sampling: {0} / {1} = {2} %.",
    n_param_accept, param_sample, (n_param_accept / param_sample) * 100
  );

  return sampled_params;
}
