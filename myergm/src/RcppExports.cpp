// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

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
arma::sp_mat change_one_link(const arma::sp_mat& adjmat);
RcppExport SEXP _myergm_change_one_link(SEXP adjmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type adjmat(adjmatSEXP);
    rcpp_result_gen = Rcpp::wrap(change_one_link(adjmat));
    return rcpp_result_gen;
END_RCPP
}
// change_one_link_modified
arma::sp_mat change_one_link_modified(const arma::sp_mat& adjmat, double& numOfEdges, double& numOfTriangles, int i, int j, int verbose);
RcppExport SEXP _myergm_change_one_link_modified(SEXP adjmatSEXP, SEXP numOfEdgesSEXP, SEXP numOfTrianglesSEXP, SEXP iSEXP, SEXP jSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type adjmat(adjmatSEXP);
    Rcpp::traits::input_parameter< double& >::type numOfEdges(numOfEdgesSEXP);
    Rcpp::traits::input_parameter< double& >::type numOfTriangles(numOfTrianglesSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(change_one_link_modified(adjmat, numOfEdges, numOfTriangles, i, j, verbose));
    return rcpp_result_gen;
END_RCPP
}
// change_all_links_of_one_node
arma::sp_mat change_all_links_of_one_node(const arma::sp_mat& adjmat);
RcppExport SEXP _myergm_change_all_links_of_one_node(SEXP adjmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type adjmat(adjmatSEXP);
    rcpp_result_gen = Rcpp::wrap(change_all_links_of_one_node(adjmat));
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
void Metropolis_Hastings(arma::sp_mat& adjmat, double& n_accepted, double& numOfEdges, double& numOfTriangles, double numOfNodes, double coefEdges, double coefTriangle, double p_large_step, double lambda, int verbose);
RcppExport SEXP _myergm_Metropolis_Hastings(SEXP adjmatSEXP, SEXP n_acceptedSEXP, SEXP numOfEdgesSEXP, SEXP numOfTrianglesSEXP, SEXP numOfNodesSEXP, SEXP coefEdgesSEXP, SEXP coefTriangleSEXP, SEXP p_large_stepSEXP, SEXP lambdaSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat& >::type adjmat(adjmatSEXP);
    Rcpp::traits::input_parameter< double& >::type n_accepted(n_acceptedSEXP);
    Rcpp::traits::input_parameter< double& >::type numOfEdges(numOfEdgesSEXP);
    Rcpp::traits::input_parameter< double& >::type numOfTriangles(numOfTrianglesSEXP);
    Rcpp::traits::input_parameter< double >::type numOfNodes(numOfNodesSEXP);
    Rcpp::traits::input_parameter< double >::type coefEdges(coefEdgesSEXP);
    Rcpp::traits::input_parameter< double >::type coefTriangle(coefTriangleSEXP);
    Rcpp::traits::input_parameter< double >::type p_large_step(p_large_stepSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Metropolis_Hastings(adjmat, n_accepted, numOfEdges, numOfTriangles, numOfNodes, coefEdges, coefTriangle, p_large_step, lambda, verbose);
    return R_NilValue;
END_RCPP
}
// create_MCMC
arma::mat create_MCMC(const arma::sp_mat& adjmat, double coefEdges, double coefTriangle, int MCMC_interval, int MCMC_samplesize, int MCMC_burnin, double p_large_step, double lambda, bool full_sample, int verbose);
RcppExport SEXP _myergm_create_MCMC(SEXP adjmatSEXP, SEXP coefEdgesSEXP, SEXP coefTriangleSEXP, SEXP MCMC_intervalSEXP, SEXP MCMC_samplesizeSEXP, SEXP MCMC_burninSEXP, SEXP p_large_stepSEXP, SEXP lambdaSEXP, SEXP full_sampleSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type adjmat(adjmatSEXP);
    Rcpp::traits::input_parameter< double >::type coefEdges(coefEdgesSEXP);
    Rcpp::traits::input_parameter< double >::type coefTriangle(coefTriangleSEXP);
    Rcpp::traits::input_parameter< int >::type MCMC_interval(MCMC_intervalSEXP);
    Rcpp::traits::input_parameter< int >::type MCMC_samplesize(MCMC_samplesizeSEXP);
    Rcpp::traits::input_parameter< int >::type MCMC_burnin(MCMC_burninSEXP);
    Rcpp::traits::input_parameter< double >::type p_large_step(p_large_stepSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< bool >::type full_sample(full_sampleSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(create_MCMC(adjmat, coefEdges, coefTriangle, MCMC_interval, MCMC_samplesize, MCMC_burnin, p_large_step, lambda, full_sample, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_myergm_count_edges_cpp", (DL_FUNC) &_myergm_count_edges_cpp, 1},
    {"_myergm_count_triangle_cpp", (DL_FUNC) &_myergm_count_triangle_cpp, 1},
    {"_myergm_change_one_link", (DL_FUNC) &_myergm_change_one_link, 1},
    {"_myergm_change_one_link_modified", (DL_FUNC) &_myergm_change_one_link_modified, 6},
    {"_myergm_change_all_links_of_one_node", (DL_FUNC) &_myergm_change_all_links_of_one_node, 1},
    {"_myergm_change_multiple_links", (DL_FUNC) &_myergm_change_multiple_links, 2},
    {"_myergm_Metropolis_Hastings", (DL_FUNC) &_myergm_Metropolis_Hastings, 10},
    {"_myergm_create_MCMC", (DL_FUNC) &_myergm_create_MCMC, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_myergm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
