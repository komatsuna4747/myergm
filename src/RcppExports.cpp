// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// run_network_sampler
arma::mat run_network_sampler(const arma::mat& adjacency_matrix, double coefEdges, double coefTwostars, double coefTriangle, int MCMC_interval, int MCMC_samplesize, int MCMC_burnin, double p_one_node_swap, double p_large_step, double p_invert, double lambda, std::string level);
RcppExport SEXP _myergm_run_network_sampler(SEXP adjacency_matrixSEXP, SEXP coefEdgesSEXP, SEXP coefTwostarsSEXP, SEXP coefTriangleSEXP, SEXP MCMC_intervalSEXP, SEXP MCMC_samplesizeSEXP, SEXP MCMC_burninSEXP, SEXP p_one_node_swapSEXP, SEXP p_large_stepSEXP, SEXP p_invertSEXP, SEXP lambdaSEXP, SEXP levelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type adjacency_matrix(adjacency_matrixSEXP);
    Rcpp::traits::input_parameter< double >::type coefEdges(coefEdgesSEXP);
    Rcpp::traits::input_parameter< double >::type coefTwostars(coefTwostarsSEXP);
    Rcpp::traits::input_parameter< double >::type coefTriangle(coefTriangleSEXP);
    Rcpp::traits::input_parameter< int >::type MCMC_interval(MCMC_intervalSEXP);
    Rcpp::traits::input_parameter< int >::type MCMC_samplesize(MCMC_samplesizeSEXP);
    Rcpp::traits::input_parameter< int >::type MCMC_burnin(MCMC_burninSEXP);
    Rcpp::traits::input_parameter< double >::type p_one_node_swap(p_one_node_swapSEXP);
    Rcpp::traits::input_parameter< double >::type p_large_step(p_large_stepSEXP);
    Rcpp::traits::input_parameter< double >::type p_invert(p_invertSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< std::string >::type level(levelSEXP);
    rcpp_result_gen = Rcpp::wrap(run_network_sampler(adjacency_matrix, coefEdges, coefTwostars, coefTriangle, MCMC_interval, MCMC_samplesize, MCMC_burnin, p_one_node_swap, p_large_step, p_invert, lambda, level));
    return rcpp_result_gen;
END_RCPP
}
// simulate_network
Rcpp::List simulate_network(const arma::mat& adjacency_matrix, double coefEdges, double coefTwostars, double coefTriangle, int MCMC_interval, int MCMC_samplesize, int MCMC_burnin, double p_one_node_swap, double p_large_step, double p_invert, double lambda, std::string level);
RcppExport SEXP _myergm_simulate_network(SEXP adjacency_matrixSEXP, SEXP coefEdgesSEXP, SEXP coefTwostarsSEXP, SEXP coefTriangleSEXP, SEXP MCMC_intervalSEXP, SEXP MCMC_samplesizeSEXP, SEXP MCMC_burninSEXP, SEXP p_one_node_swapSEXP, SEXP p_large_stepSEXP, SEXP p_invertSEXP, SEXP lambdaSEXP, SEXP levelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type adjacency_matrix(adjacency_matrixSEXP);
    Rcpp::traits::input_parameter< double >::type coefEdges(coefEdgesSEXP);
    Rcpp::traits::input_parameter< double >::type coefTwostars(coefTwostarsSEXP);
    Rcpp::traits::input_parameter< double >::type coefTriangle(coefTriangleSEXP);
    Rcpp::traits::input_parameter< int >::type MCMC_interval(MCMC_intervalSEXP);
    Rcpp::traits::input_parameter< int >::type MCMC_samplesize(MCMC_samplesizeSEXP);
    Rcpp::traits::input_parameter< int >::type MCMC_burnin(MCMC_burninSEXP);
    Rcpp::traits::input_parameter< double >::type p_one_node_swap(p_one_node_swapSEXP);
    Rcpp::traits::input_parameter< double >::type p_large_step(p_large_stepSEXP);
    Rcpp::traits::input_parameter< double >::type p_invert(p_invertSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< std::string >::type level(levelSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_network(adjacency_matrix, coefEdges, coefTwostars, coefTriangle, MCMC_interval, MCMC_samplesize, MCMC_burnin, p_one_node_swap, p_large_step, p_invert, lambda, level));
    return rcpp_result_gen;
END_RCPP
}
// proposal_mvrnorm
arma::vec proposal_mvrnorm(arma::vec mu, arma::mat sigma);
RcppExport SEXP _myergm_proposal_mvrnorm(SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(proposal_mvrnorm(mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// run_DMH_cpp
arma::mat run_DMH_cpp(const arma::mat& adjacency_matrix, arma::mat sigma_proposal, arma::vec theta_init, arma::mat sigma_prior, arma::vec terms_included, int param_sample, int MCMC_network_interval, double p_one_node_swap, double p_large_step, double p_invert, double lambda, std::string level);
RcppExport SEXP _myergm_run_DMH_cpp(SEXP adjacency_matrixSEXP, SEXP sigma_proposalSEXP, SEXP theta_initSEXP, SEXP sigma_priorSEXP, SEXP terms_includedSEXP, SEXP param_sampleSEXP, SEXP MCMC_network_intervalSEXP, SEXP p_one_node_swapSEXP, SEXP p_large_stepSEXP, SEXP p_invertSEXP, SEXP lambdaSEXP, SEXP levelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type adjacency_matrix(adjacency_matrixSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_proposal(sigma_proposalSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta_init(theta_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma_prior(sigma_priorSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type terms_included(terms_includedSEXP);
    Rcpp::traits::input_parameter< int >::type param_sample(param_sampleSEXP);
    Rcpp::traits::input_parameter< int >::type MCMC_network_interval(MCMC_network_intervalSEXP);
    Rcpp::traits::input_parameter< double >::type p_one_node_swap(p_one_node_swapSEXP);
    Rcpp::traits::input_parameter< double >::type p_large_step(p_large_stepSEXP);
    Rcpp::traits::input_parameter< double >::type p_invert(p_invertSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< std::string >::type level(levelSEXP);
    rcpp_result_gen = Rcpp::wrap(run_DMH_cpp(adjacency_matrix, sigma_proposal, theta_init, sigma_prior, terms_included, param_sample, MCMC_network_interval, p_one_node_swap, p_large_step, p_invert, lambda, level));
    return rcpp_result_gen;
END_RCPP
}
// count_edges_cpp
double count_edges_cpp(const arma::sp_mat& adjmat);
RcppExport SEXP _myergm_count_edges_cpp(SEXP adjmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type adjmat(adjmatSEXP);
    rcpp_result_gen = Rcpp::wrap(count_edges_cpp(adjmat));
    return rcpp_result_gen;
END_RCPP
}
// count_triangle_cpp
double count_triangle_cpp(const arma::sp_mat& adjmat);
RcppExport SEXP _myergm_count_triangle_cpp(SEXP adjmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type adjmat(adjmatSEXP);
    rcpp_result_gen = Rcpp::wrap(count_triangle_cpp(adjmat));
    return rcpp_result_gen;
END_RCPP
}
// change_one_link
arma::sp_mat change_one_link(const arma::sp_mat& adjmat, double& numOfEdges, double& numOfTriangles, int i, int j, int verbose);
RcppExport SEXP _myergm_change_one_link(SEXP adjmatSEXP, SEXP numOfEdgesSEXP, SEXP numOfTrianglesSEXP, SEXP iSEXP, SEXP jSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type adjmat(adjmatSEXP);
    Rcpp::traits::input_parameter< double& >::type numOfEdges(numOfEdgesSEXP);
    Rcpp::traits::input_parameter< double& >::type numOfTriangles(numOfTrianglesSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(change_one_link(adjmat, numOfEdges, numOfTriangles, i, j, verbose));
    return rcpp_result_gen;
END_RCPP
}
// change_all_links_of_one_node
arma::sp_mat change_all_links_of_one_node(const arma::sp_mat& adjmat, int i);
RcppExport SEXP _myergm_change_all_links_of_one_node(SEXP adjmatSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type adjmat(adjmatSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(change_all_links_of_one_node(adjmat, i));
    return rcpp_result_gen;
END_RCPP
}
// change_multiple_links
arma::sp_mat change_multiple_links(const arma::sp_mat& adjmat, double lambda);
RcppExport SEXP _myergm_change_multiple_links(SEXP adjmatSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type adjmat(adjmatSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(change_multiple_links(adjmat, lambda));
    return rcpp_result_gen;
END_RCPP
}
// Metropolis_Hastings
void Metropolis_Hastings(arma::sp_mat& adjmat, double& n_accepted, double& n_one_node_swap, double& n_accepted_one_node_swap, double& n_large_step, double& n_accepted_large_step, double& n_invert, double& n_accepted_invert, double& numOfEdges, double& numOfTriangles, double numOfNodes, double coefEdges, double coefTriangle, double p_one_node_swap, double p_large_step, double p_invert, double lambda, int verbose);
RcppExport SEXP _myergm_Metropolis_Hastings(SEXP adjmatSEXP, SEXP n_acceptedSEXP, SEXP n_one_node_swapSEXP, SEXP n_accepted_one_node_swapSEXP, SEXP n_large_stepSEXP, SEXP n_accepted_large_stepSEXP, SEXP n_invertSEXP, SEXP n_accepted_invertSEXP, SEXP numOfEdgesSEXP, SEXP numOfTrianglesSEXP, SEXP numOfNodesSEXP, SEXP coefEdgesSEXP, SEXP coefTriangleSEXP, SEXP p_one_node_swapSEXP, SEXP p_large_stepSEXP, SEXP p_invertSEXP, SEXP lambdaSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type adjmat(adjmatSEXP);
    Rcpp::traits::input_parameter< double& >::type n_accepted(n_acceptedSEXP);
    Rcpp::traits::input_parameter< double& >::type n_one_node_swap(n_one_node_swapSEXP);
    Rcpp::traits::input_parameter< double& >::type n_accepted_one_node_swap(n_accepted_one_node_swapSEXP);
    Rcpp::traits::input_parameter< double& >::type n_large_step(n_large_stepSEXP);
    Rcpp::traits::input_parameter< double& >::type n_accepted_large_step(n_accepted_large_stepSEXP);
    Rcpp::traits::input_parameter< double& >::type n_invert(n_invertSEXP);
    Rcpp::traits::input_parameter< double& >::type n_accepted_invert(n_accepted_invertSEXP);
    Rcpp::traits::input_parameter< double& >::type numOfEdges(numOfEdgesSEXP);
    Rcpp::traits::input_parameter< double& >::type numOfTriangles(numOfTrianglesSEXP);
    Rcpp::traits::input_parameter< double >::type numOfNodes(numOfNodesSEXP);
    Rcpp::traits::input_parameter< double >::type coefEdges(coefEdgesSEXP);
    Rcpp::traits::input_parameter< double >::type coefTriangle(coefTriangleSEXP);
    Rcpp::traits::input_parameter< double >::type p_one_node_swap(p_one_node_swapSEXP);
    Rcpp::traits::input_parameter< double >::type p_large_step(p_large_stepSEXP);
    Rcpp::traits::input_parameter< double >::type p_invert(p_invertSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Metropolis_Hastings(adjmat, n_accepted, n_one_node_swap, n_accepted_one_node_swap, n_large_step, n_accepted_large_step, n_invert, n_accepted_invert, numOfEdges, numOfTriangles, numOfNodes, coefEdges, coefTriangle, p_one_node_swap, p_large_step, p_invert, lambda, verbose);
    return R_NilValue;
END_RCPP
}
// create_MCMC
arma::mat create_MCMC(const arma::sp_mat& adjmat, double coefEdges, double coefTriangle, int MCMC_interval, int MCMC_samplesize, int MCMC_burnin, double p_one_node_swap, double p_large_step, double p_invert, double lambda, bool full_sample, int verbose);
RcppExport SEXP _myergm_create_MCMC(SEXP adjmatSEXP, SEXP coefEdgesSEXP, SEXP coefTriangleSEXP, SEXP MCMC_intervalSEXP, SEXP MCMC_samplesizeSEXP, SEXP MCMC_burninSEXP, SEXP p_one_node_swapSEXP, SEXP p_large_stepSEXP, SEXP p_invertSEXP, SEXP lambdaSEXP, SEXP full_sampleSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type adjmat(adjmatSEXP);
    Rcpp::traits::input_parameter< double >::type coefEdges(coefEdgesSEXP);
    Rcpp::traits::input_parameter< double >::type coefTriangle(coefTriangleSEXP);
    Rcpp::traits::input_parameter< int >::type MCMC_interval(MCMC_intervalSEXP);
    Rcpp::traits::input_parameter< int >::type MCMC_samplesize(MCMC_samplesizeSEXP);
    Rcpp::traits::input_parameter< int >::type MCMC_burnin(MCMC_burninSEXP);
    Rcpp::traits::input_parameter< double >::type p_one_node_swap(p_one_node_swapSEXP);
    Rcpp::traits::input_parameter< double >::type p_large_step(p_large_stepSEXP);
    Rcpp::traits::input_parameter< double >::type p_invert(p_invertSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< bool >::type full_sample(full_sampleSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(create_MCMC(adjmat, coefEdges, coefTriangle, MCMC_interval, MCMC_samplesize, MCMC_burnin, p_one_node_swap, p_large_step, p_invert, lambda, full_sample, verbose));
    return rcpp_result_gen;
END_RCPP
}
// edit_spmat
arma::sp_mat edit_spmat(int N);
RcppExport SEXP _myergm_edit_spmat(SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(edit_spmat(N));
    return rcpp_result_gen;
END_RCPP
}
// edit_mat
arma::mat edit_mat(int N);
RcppExport SEXP _myergm_edit_mat(SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(edit_mat(N));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_myergm_run_network_sampler", (DL_FUNC) &_myergm_run_network_sampler, 12},
    {"_myergm_simulate_network", (DL_FUNC) &_myergm_simulate_network, 12},
    {"_myergm_proposal_mvrnorm", (DL_FUNC) &_myergm_proposal_mvrnorm, 2},
    {"_myergm_run_DMH_cpp", (DL_FUNC) &_myergm_run_DMH_cpp, 12},
    {"_myergm_count_edges_cpp", (DL_FUNC) &_myergm_count_edges_cpp, 1},
    {"_myergm_count_triangle_cpp", (DL_FUNC) &_myergm_count_triangle_cpp, 1},
    {"_myergm_change_one_link", (DL_FUNC) &_myergm_change_one_link, 6},
    {"_myergm_change_all_links_of_one_node", (DL_FUNC) &_myergm_change_all_links_of_one_node, 2},
    {"_myergm_change_multiple_links", (DL_FUNC) &_myergm_change_multiple_links, 2},
    {"_myergm_Metropolis_Hastings", (DL_FUNC) &_myergm_Metropolis_Hastings, 18},
    {"_myergm_create_MCMC", (DL_FUNC) &_myergm_create_MCMC, 12},
    {"_myergm_edit_spmat", (DL_FUNC) &_myergm_edit_spmat, 1},
    {"_myergm_edit_mat", (DL_FUNC) &_myergm_edit_mat, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_myergm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
