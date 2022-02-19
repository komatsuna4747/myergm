#ifndef __MCMC__
#define __MCMC__

double count_edges(const arma::mat& adjmat) {
  double edges = accu(adjmat);
  return edges/2;
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
  double coefTriangle;
  int MCMC_interval;
  int MCMC_samplesize;
  int MCMC_burnin;
  double p_one_node_swap ;
  double p_large_step;
  double p_invert;
  double lambda;
  int verbose;
  // For internal use
  arma::mat netstats;
  arma::mat adjmat_next;
  double numOfEdges;
  double numOfEdges_next;
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
    double coef_triangle,
    int mcmc_interval,
    int mcmc_samplesize,
    int mcmc_burnin,
    double prob_one_node_swap,
    double prob_large_step,
    double prob_invert,
    double lambda_for_large_step,
    int log_verbose
  ){
    adjmat = adjacency_matrix;
    adjmat_next = adjacency_matrix;
    coefEdges = coef_edges;
    coefTriangle = coef_triangle;
    MCMC_interval = mcmc_interval;
    MCMC_samplesize = mcmc_samplesize;
    MCMC_burnin = mcmc_burnin;
    p_one_node_swap = prob_one_node_swap;
    p_large_step = prob_large_step;
    p_invert = prob_invert;
    lambda = lambda_for_large_step;
    verbose = log_verbose;
    numOfEdges = count_edges(adjmat);
    numOfTriangles = count_triangle(adjmat);
    netstats = arma::zeros(MCMC_samplesize, 2);
  }

  void run_simulation() {
    // Total number of iterations
    int total_iter = MCMC_burnin + MCMC_interval * MCMC_samplesize;

    for (int i = 0; i < total_iter; i++) {
      run_update_stat();

      // Store the stats
      if (i >= MCMC_burnin && (i - MCMC_burnin) % MCMC_interval == 0) {
        int index = (i - MCMC_burnin) / MCMC_interval;
        netstats(index, 0) = numOfEdges;
        netstats(index, 1) = numOfTriangles;
      }

      // Check the current state if necessary
      if (verbose >= 5) {
        Rcpp::Rcout << "edges:" << numOfEdges << ", triangle:" << numOfTriangles << "\n";
      }
    }

    if (verbose >= 1) {
      display_acceptance_rate(total_iter);
    }
  }

  void display_acceptance_rate(double total_iter) {
    double acceptance_rate = (n_accepted / total_iter) * 100;
    Rcpp::Rcout << "Acceptance rate: " << n_accepted << "/" << total_iter << " = " << acceptance_rate << " %." << "\n";
    // Check large steps
    if (n_one_node_swap >= 1) {
      double one_node_swap_rate = (n_accepted_one_node_swap / n_one_node_swap) * 100;
      Rcpp::Rcout << "Swapping all links of one node: " << n_accepted_one_node_swap << "/" << n_one_node_swap << " = " << one_node_swap_rate << " %." << "\n";
    }
    if (n_large_step >= 1) {
      double large_step_rate = (n_accepted_large_step / n_large_step) * 100;
      Rcpp::Rcout << "Swapping multiple links in one step: " << n_accepted_large_step << "/" << n_large_step << " = " << large_step_rate << " %." << "\n";
    }
    if (n_invert >= 1) {
      double invert_rate = (n_accepted_invert / n_invert) * 100;
      Rcpp::Rcout << "Inverting the adjacency matrix: " << n_accepted_invert << "/" << n_invert << " = " << invert_rate << " %." << "\n";
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

  void update_adjmat_and_stats(arma::mat G, double edges, double triangle) {
    adjmat_next = G;
    numOfEdges_next = edges;
    numOfTriangles_next = triangle;
  }

  void update_one_link() {
    arma::mat G = adjmat;
    double edges = numOfEdges;
    double triangle = numOfTriangles;

    // pick a random element
    int i = unif_rand() * G.n_rows;
    // pick a random element from what's left (there is one fewer to choose from)...
    int j =  unif_rand() * (G.n_rows- 1);
    // ...and adjust second choice to take into account the first choice
    if (j >= i)
    {
      ++j;
    }

    // Change the state of the selected dyad.
    G(i, j) = 1 - G(i, j);
    G(j, i) = 1 - G(j, i);

    // If added a link, add one to numOfEdges. Otherwise, subtract one.
    edges += G(i, j) == 1 ? 1 : -1;

    // Count the number of triangles that appear/disappear by adding/removing the link.
    arma::vec col_i = G.col(i);
    arma::vec col_j = G.col(j);
    double ij = dot(col_i, col_j);
    triangle += G(i, j) == 1 ? ij : -ij;

    update_adjmat_and_stats(G, edges, triangle);
    run_proposal_step();
    if (accept_flag == 1) {
      accept_flag = 0;
    };
  }

  void update_all_links_of_one_node() {
    arma::mat G = adjmat;

    // pick a random element
    int i = unif_rand() * G.n_rows;

    // Change all links of one node
    arma::vec G_i = adjmat.col(i);
    G_i = 1 - G_i;
    G.col(i) = G_i;
    G.row(i) = G_i.t();

    // Set diagonal entries to be zero
    G.diag().zeros();

    update_adjmat_and_stats(G, count_edges(G), count_triangle(G));
    n_one_node_swap += 1;

    run_proposal_step();
    if (accept_flag == 1) {
      n_accepted_one_node_swap += 1;
      accept_flag = 0;
    };
  }

  void change_multiple_links() {
    arma::mat G = adjmat;

    // Change [lambda * numOfNodes] links simultaneously.
    int numOfChanges = lambda * G.n_rows;

    // Select which link to swap.
    // This part could be written in a more elegant and efficient way.
    arma::mat Dyads = arma::zeros(numOfChanges, 2);
    for (int l = 0; l < numOfChanges; l++) {
      // pick a random element
      int i = unif_rand() * G.n_rows;
      // pick a random element from what's left (there is one fewer to choose from)...
      int j =  unif_rand() * (G.n_rows- 1);
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
      G(i, j) = 1 - G(i, j);
      G(j, i) = 1 - G(j, i);
    }
    update_adjmat_and_stats(G, count_edges(G), count_triangle(G));
    n_large_step += 1;

    run_proposal_step();
    if (accept_flag == 1) {
      n_accepted_large_step += 1;
      accept_flag = 0;
    };
  }

  void invert_adjacency_matrix() {
    // Invert the adjacency matrix and set the diagonals to be zero
    arma::mat G = 1 - adjmat;
    G.diag().zeros();
    update_adjmat_and_stats(G, count_edges(G), count_triangle(G));
    n_invert += 1;

    run_proposal_step();
    if (accept_flag == 1) {
      n_accepted_invert += 1;
      accept_flag = 0;
    };
  }

  void run_proposal_step() {
    // Store the results.
    double changestat_edges = numOfEdges_next - numOfEdges;
    double changestat_triangle = numOfTriangles_next - numOfTriangles;

    // Compute the logratio
    double logratio = changestat_edges * coefEdges + changestat_triangle * coefTriangle;

    // Accept or reject?
    if (logratio >= 0.0 || log(unif_rand()) < logratio) {
      if (verbose >= 5) {
        Rcpp::Rcout << "Accepted." << "\n";
      }
      accept_flag = 1;
      n_accepted += 1;
      adjmat = adjmat_next;
      numOfEdges = numOfEdges_next;
      numOfTriangles = numOfTriangles_next;
    } else {
      if (verbose >= 5) {
        Rcpp::Rcout << "Rejected." << "\n";
      }
    }
  }
};

#endif // __MCMC__
