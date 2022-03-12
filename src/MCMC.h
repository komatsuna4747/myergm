#ifndef __MCMC__
#define __MCMC__
#include <RcppSpdlog>

double count_edges(const arma::mat& adjmat) {
  double edges = accu(adjmat);
  return edges/2;
}


double count_twostar(const arma::mat& adjmat) {
  arma::colvec colsum = arma::sum(adjmat, 1);
  return dot(colsum, colsum - 1) / 2;
}


double count_triangle(const arma::mat& adjmat) {
  double triangle = 0;
  int numOfNodes = adjmat.n_rows;
  for (int i = 0; i < numOfNodes-2; i++) {
    for (int j = i+1; j < numOfNodes-1; j++) {
      for (int k = j+1; k < numOfNodes; k++) {
        triangle += adjmat(i, j) * adjmat(j, k) * adjmat(k, i);
      }
    }
  }
  return triangle;
}


class NetworkSampler
{
public:
  arma::mat adjmat;
  double coefEdges;
  double coefTwostars;
  double coefTriangle;
  int MCMC_interval;
  int MCMC_samplesize;
  int MCMC_burnin;
  double p_one_node_swap ;
  double p_large_step;
  double p_invert;
  double lambda;
  // For internal use
  arma::mat netstats;
  double numOfEdges;
  double numOfEdges_next;
  double numOfTwostars;
  double numOfTwostars_next;
  double numOfTriangles;
  double numOfTriangles_next;
  double n_accepted = 0;
  double n_one_node_swap = 0;
  double n_accepted_one_node_swap = 0;
  double n_large_step = 0;
  double n_accepted_large_step = 0;
  double n_invert = 0;
  double n_accepted_invert = 0;
  int accept_flag = 0;

  NetworkSampler(
    arma::mat adjacency_matrix,
    double coef_edges,
    double coef_twostars,
    double coef_triangle,
    int mcmc_interval,
    int mcmc_samplesize,
    int mcmc_burnin,
    double prob_one_node_swap,
    double prob_large_step,
    double prob_invert,
    double lambda_for_large_step
  ){
    adjmat = adjacency_matrix;
    coefEdges = coef_edges;
    coefTwostars = coef_twostars;
    coefTriangle = coef_triangle;
    MCMC_interval = mcmc_interval;
    MCMC_samplesize = mcmc_samplesize;
    MCMC_burnin = mcmc_burnin;
    p_one_node_swap = prob_one_node_swap;
    p_large_step = prob_large_step;
    p_invert = prob_invert;
    lambda = lambda_for_large_step;
    numOfEdges = count_edges(adjmat);
    numOfTwostars = count_twostar(adjmat);
    numOfTriangles = count_triangle(adjmat);
    netstats = arma::zeros(MCMC_samplesize, 3);
  }

  void run_simulation() {
    // Total number of iterations
    int total_iter = MCMC_burnin + MCMC_interval * MCMC_samplesize;

    for (int i = 0; i < total_iter; i++) {
      run_update_stat();

      // Store the stats
      if (i >= MCMC_burnin - 1 && (i - MCMC_burnin + 1) % MCMC_interval == 0) {
        int index = (i - MCMC_burnin) / MCMC_interval;
        netstats(index, 0) = numOfEdges;
        netstats(index, 1) = numOfTwostars;
        netstats(index, 2) = numOfTriangles;
      }

      // Check the current state if necessary
      spdlog::trace("edges: {0}, twostars: {1}, triangle: {2}", numOfEdges, numOfTwostars, numOfTriangles);
    }
    display_acceptance_rate(total_iter);
  }

  void display_acceptance_rate(double total_iter) {
    double acceptance_rate = (n_accepted / total_iter) * 100;
    spdlog::debug(
      "Acceptance rate: {0} / {1} = {2} %.",
      n_accepted, total_iter, acceptance_rate
    );
    // Check large steps
    if (n_one_node_swap >= 1) {
      double one_node_swap_rate = (n_accepted_one_node_swap / n_one_node_swap) * 100;
      spdlog::debug(
        "Swapping all links of one node: {0} / {1} = {2} %.",
        n_accepted_one_node_swap, n_one_node_swap, one_node_swap_rate
      );
    }
    if (n_large_step >= 1) {
      double large_step_rate = (n_accepted_large_step / n_large_step) * 100;
      spdlog::debug(
        "Swapping multiple links in one step: {0} / {1} = {2} %.",
        n_accepted_large_step, n_large_step, large_step_rate
      );
    }
    if (n_invert >= 1) {
      double invert_rate = (n_accepted_invert / n_invert) * 100;
      spdlog::debug(
        "Inverting the adjacency matrix: {0} / {1} = {2} %.",
        n_accepted_invert, n_invert, invert_rate
      );
    }
  }

  void run_update_stat() {
    // Draw a number from (0, 1).
    double x = unif_rand();

    if (x < p_one_node_swap) {
      update_all_links_of_one_node();
      return;
    }
    if (x >= p_one_node_swap && x < p_one_node_swap + p_large_step) {
      change_multiple_links();
      return;
    }
    if (x >= p_one_node_swap + p_large_step && x < p_one_node_swap + p_large_step + p_invert) {
      invert_adjacency_matrix();
      return;
    }
    update_one_link();
  }

  void update_stats(double edges, double twostars, double triangle) {
    numOfEdges_next = edges;
    numOfTwostars_next = twostars;
    numOfTriangles_next = triangle;
  }

  void update_one_link() {
    // pick a random element
    int i = unif_rand() * adjmat.n_rows;
    // pick a random element from what's left (there is one fewer to choose from)...
    int j =  unif_rand() * (adjmat.n_rows- 1);
    // ...and adjust second choice to take into account the first choice
    if (j >= i)
    {
      ++j;
    }

    double edges_i_before = sum(adjmat.col(i));
    double edges_j_before = sum(adjmat.col(j));

    // Change the state of the selected dyad.
    adjmat(i, j) = 1 - adjmat(i, j);
    adjmat(j, i) = 1 - adjmat(j, i);

    // If added a link, add one to numOfEdges. Otherwise, subtract one.
    double edgediff = (adjmat(i, j) == 1) ? 1 : -1;

    // Changes in twostars
    double edges_i = sum(adjmat.col(i));
    double edges_j = sum(adjmat.col(j));
    double stardiff = edges_i * (edges_i - 1) / 2 - edges_i_before * (edges_i_before - 1) / 2 + edges_j * (edges_j - 1) / 2 - edges_j_before * (edges_j_before - 1) / 2;

    // Count the number of triangles that appear/disappear by adding/removing the link.
    double ij = dot(adjmat.col(i), adjmat.col(j));
    double trianglediff = (adjmat(i, j) == 1) ? ij : -ij;

    update_stats(numOfEdges + edgediff, numOfTwostars + stardiff, numOfTriangles + trianglediff);
    //update_stats(numOfEdges + edgediff, numOfTwostars + twostardiff, numOfTriangles + trianglediff);
    //update_stats(numOfEdges + edgediff, count_twostar(adjmat), count_triangle(adjmat));
    run_proposal_step();
    if (accept_flag == 1) {
      accept_flag = 0;
    } else {
      adjmat(i, j) = 1 - adjmat(i, j);
      adjmat(j, i) = 1 - adjmat(j, i);
    }
  }

  void update_all_links_of_one_node() {
    // pick a random element
    int i = unif_rand() * adjmat.n_rows;

    // Change all links of one node
    arma::vec G_i = adjmat.col(i);
    G_i = 1 - G_i;
    adjmat.col(i) = G_i;
    adjmat.row(i) = G_i.t();

    // Set diagonal entries to be zero
    adjmat.diag().zeros();

    update_stats(count_edges(adjmat), count_twostar(adjmat), count_triangle(adjmat));
    n_one_node_swap += 1;

    run_proposal_step();
    if (accept_flag == 1) {
      n_accepted_one_node_swap += 1;
      accept_flag = 0;
    } else {
      G_i = 1 - G_i;
      adjmat.col(i) = G_i;
      adjmat.row(i) = G_i.t();
      adjmat.diag().zeros();
    };
  }

  void change_multiple_links() {
    // Change [lambda * numOfNodes] links simultaneously.
    int numOfChanges = lambda * adjmat.n_rows;

    // Select which link to swap.
    // This part could be written in a more elegant and efficient way.
    arma::mat Dyads = arma::zeros(numOfChanges, 2);
    for (int l = 0; l < numOfChanges; l++) {
      // pick a random element
      int i = unif_rand() * adjmat.n_rows;
      // pick a random element from what's left (there is one fewer to choose from)...
      int j =  unif_rand() * (adjmat.n_rows- 1);
      // ...and adjust second choice to take into account the first choice
      if (j >= i)
      {
        ++j;
      }
      // Store i and j in Dyads.
      Dyads(l, 0) = i;
      Dyads(l, 1) = j;
    }

    // swap links of the selected dyads.
    for (int l = 0; l < numOfChanges; l++) {
      int i = Dyads(l, 0);
      int j = Dyads(l, 1);
      adjmat(i, j) = 1 - adjmat(i, j);
      adjmat(j, i) = 1 - adjmat(j, i);
    }
    update_stats(count_edges(adjmat), count_twostar(adjmat), count_triangle(adjmat));
    n_large_step += 1;

    run_proposal_step();
    if (accept_flag == 1) {
      n_accepted_large_step += 1;
      accept_flag = 0;
    } else {
      for (int l = 0; l < numOfChanges; l++) {
        int i = Dyads(l, 0);
        int j = Dyads(l, 1);
        adjmat(i, j) = 1 - adjmat(i, j);
        adjmat(j, i) = 1 - adjmat(j, i);
      }
    };
  }

  void invert_adjacency_matrix() {
    // Invert the adjacency matrix and set the diagonals to be zero
    adjmat = 1 - adjmat;
    adjmat.diag().zeros();
    update_stats(count_edges(adjmat), count_twostar(adjmat), count_triangle(adjmat));
    n_invert += 1;

    run_proposal_step();
    if (accept_flag == 1) {
      n_accepted_invert += 1;
      accept_flag = 0;
    } else {
      adjmat = 1 - adjmat;
      adjmat.diag().zeros();
    }
  }

  void run_proposal_step() {
    // Store the results.
    double changestat_edges = numOfEdges_next - numOfEdges;
    double changestat_twostars = numOfTwostars_next - numOfTwostars;
    double changestat_triangle = numOfTriangles_next - numOfTriangles;

    // Compute the logratio
    double logratio = changestat_edges * coefEdges + changestat_twostars * coefTwostars +changestat_triangle * coefTriangle;

    // Accept or reject?
    if (logratio >= 0.0 || log(unif_rand()) < logratio) {
      spdlog::trace("Accepted.");
      accept_flag = 1;
      n_accepted += 1;
      numOfEdges = numOfEdges_next;
      numOfTwostars = numOfTwostars_next;
      numOfTriangles = numOfTriangles_next;
    } else {
      spdlog::trace("Rejected.");
    }
  }

  void reset_parameters(arma::mat adjmat_original) {
    n_accepted = 0;
    n_accepted_invert = 0;
    n_accepted_large_step = 0;
    n_accepted_one_node_swap = 0;
    n_invert = 0;
    n_large_step = 0;
    n_one_node_swap = 0;
    adjmat = adjmat_original;
    numOfEdges = count_edges(adjmat_original);
    numOfTwostars = count_twostar(adjmat_original);
    numOfTriangles = count_triangle(adjmat_original);
    accept_flag = 0;
  }
};

#endif // __MCMC__
