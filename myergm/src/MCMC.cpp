#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Auxilirary functions
// [[Rcpp::export]]
double count_edges_cpp(const arma::sp_mat& adjmat) {
  double edges = 0;
  for (arma::sp_mat::const_iterator i = adjmat.begin(); i != adjmat.end(); ++i) {
    edges += *i;
  }
  return edges/2;
}


// [[Rcpp::export]]
double count_triangle_cpp(const arma::sp_mat& adjmat) {
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


// [[Rcpp::export]]
arma::sp_mat change_one_link(const arma::sp_mat& adjmat) {
  arma::sp_mat adjmat_next = adjmat;
  double numOfNodes = adjmat.n_rows;

  // Choose one dyad randomly and change the state.
  // pick a random element
  int i = unif_rand() * numOfNodes;
  // pick a random element from what's left (there is one fewer to choose from)...
  int j =  unif_rand() * (numOfNodes- 1);
  // ...and adjust second choice to take into account the first choice
  if (j >= i)
  {
    ++i;
  }
  adjmat_next(i, j) = 1 - adjmat_next(i, j);
  adjmat_next(j, i) = 1 - adjmat_next(j, i);

  return adjmat_next;
}


// [[Rcpp::export]]
arma::sp_mat change_all_links_of_one_node(const arma::sp_mat& adjmat) {
  arma::sp_mat adjmat_next = adjmat;
  double numOfNodes = adjmat.n_rows;

  // Choose one node randomly.
  // pick a random element
  int i = unif_rand() * numOfNodes;
  arma::sp_mat G_i = adjmat.col(i);
  G_i = 1 - G_i;
  adjmat_next.col(i) = G_i;
  adjmat_next.row(i) = G_i.t();

  // Set diagonal entries to be zero
  adjmat_next.diag().zeros();

  // Return the output
  return adjmat_next;
}


// [[Rcpp::export]]
arma::sp_mat change_multiple_links(const arma::sp_mat& adjmat, double lambda) {
  arma::sp_mat adjmat_next = adjmat;
  double numOfNodes = adjmat.n_rows;

  // Change [lambda * numOfNodes] links simultaneously.
  int numOfChanges = lambda * numOfNodes;

  // Select which link to flip.
  // This part could be written in a more elegant and efficient way.
  arma::mat Dyads = arma::zeros(numOfChanges, 2);
  for (int l = 0; l < numOfChanges; l++) {
    // pick a random element
    int i = unif_rand() * numOfNodes;
    // pick a random element from what's left (there is one fewer to choose from)...
    int j =  unif_rand() * (numOfNodes- 1);
    // ...and adjust second choice to take into account the first choice
    if (j >= i)
    {
      ++i;
    }
    // Store i and j in Dyads.
    Dyads(l, 0) = i;
    Dyads(l, 1) = j;
  }

  // Flip links of the selected dyads.
  for (int l = 0; l < numOfChanges; l++) {
    int i = Dyads(l, 0);
    int j = Dyads(l, 1);
    adjmat_next(i, j) = 1 - adjmat_next(i, j);
    adjmat_next(j, i) = 1 - adjmat_next(j, i);
  }

  // Return the output
  return adjmat_next;
}


// [[Rcpp::export]]
void Metropolis_Hastings(arma::sp_mat& adjmat,
                         double& n_accepted,
                         double coefEdges,
                         double coefTriangle,
                         double p_large_step,
                         double lambda,
                         int verbose) {
  // Draw a number from (0, 1).
  double x = unif_rand();

  // Initialization
  double numOfEdges = count_edges_cpp(adjmat);
  double numOfTriangles = count_triangle_cpp(adjmat);

  // If x < p_large_step, then select one node and flip its all links.
  if (x < p_large_step) {
    if (verbose >= 5) {
      Rcpp::Rcout << "Flipping all the links of one node..." << "\n";
      }
    // Update all the links of one node.
    arma::sp_mat adjmat_next = change_all_links_of_one_node(adjmat);

    // Store the results.
    double numOfEdges_next = count_edges_cpp(adjmat_next);
    double numOfTriangles_next = count_triangle_cpp(adjmat_next);
    double changestat_edges = numOfEdges_next - numOfEdges;
    double changestat_triangle = numOfTriangles_next - numOfTriangles;

    // Compute the logratio
    double logratio = changestat_edges * coefEdges + changestat_triangle * coefTriangle;

    // Accept or reject?
    if (logratio >= 0.0 || log(unif_rand()) < logratio) {
      if (verbose >= 5) {
        Rcpp::Rcout << "Accepted." << "\n";
      }
      n_accepted += 1;
      adjmat = adjmat_next;
    } else {
      if (verbose >= 5) {
        Rcpp::Rcout << "Rejected." << "\n";
      }
      adjmat_next = adjmat;
    }
  }

  // If p_large_step <= x < 2 * p_large_step, then flip [lambda * n] links simultaneously.
  else if (x >= p_large_step && x < 2 * p_large_step) {
    if (verbose >= 5) {
      Rcpp::Rcout << "Flipping lambda * N links simultaneously..." << "\n";
    }
    // Update [lambda * n] links.
    arma::sp_mat adjmat_next = change_multiple_links(adjmat, lambda);

    // Store the results.
    double numOfEdges_next = count_edges_cpp(adjmat_next);
    double numOfTriangles_next = count_triangle_cpp(adjmat_next);
    double changestat_edges = numOfEdges_next - numOfEdges;
    double changestat_triangle = numOfTriangles_next - numOfTriangles;

    // Compute the logratio
    double logratio = changestat_edges * coefEdges + changestat_triangle * coefTriangle;

    // Accept or reject?
    if (logratio >= 0.0 || log(unif_rand()) < logratio) {
      if (verbose >= 5) {
        Rcpp::Rcout << "Accepted." << "\n";
      }
      n_accepted += 1;
      adjmat = adjmat_next;
    } else {
      if (verbose >= 5) {
        Rcpp::Rcout << "Rejected." << "\n";
      }
    }
  }

  // Flip one link randomply.
  else {
    if (verbose >= 5) {
      Rcpp::Rcout << "Flipping one link..." << "\n";
    }
    // Store the results.
    arma::sp_mat adjmat_next = change_one_link(adjmat);

    // Store the results.
    double numOfEdges_next = count_edges_cpp(adjmat_next);
    double numOfTriangles_next = count_triangle_cpp(adjmat_next);
    double changestat_edges = numOfEdges_next - numOfEdges;
    double changestat_triangle = numOfTriangles_next - numOfTriangles;

    // Compute the logratio
    double logratio = changestat_edges * coefEdges + changestat_triangle * coefTriangle;

    // Accept or reject?
    double random_number = unif_rand();
    if (logratio >= 0.0 || log(random_number) < logratio) {
      if (verbose >= 5) {
        Rcpp::Rcout << "Accepted." << "\n";
      }
      n_accepted += 1;
      adjmat = adjmat_next;
    } else {
      if (verbose >= 5) {
        Rcpp::Rcout << "Rejected." << "\n";
      }
    }
  }
}


// Main wrapper function for MCMC.
// [[Rcpp::export]]
arma::mat create_MCMC(const arma::sp_mat& adjmat,
                      double coefEdges,
                      double coefTriangle,
                      int MCMC_interval = 1024,
                      int MCMC_samplesize = 1024,
                      double MCMC_burnin = 1024 * 16,
                      double p_large_step = 0.01,
                      double lambda = 0.5,
                      int verbose = 0) {
  // Total number of iterations
  int total_iter = MCMC_burnin + MCMC_interval * MCMC_samplesize;

  // Create a matrix in which network stats used for MC-MLE are stored.
  // arma::mat stats_for_MCMLE = arma::zeros(MCMC.samplesize, 2);

  // Create a matrix in which all network stats during MCMC are stored.
  arma::mat all_stats = arma::zeros(total_iter , 2);

  // Initialization
  arma::sp_mat G = adjmat;
  double n_accepted = 0;

  // Start MCMC
  for (int i = 0; i < total_iter; i++) {
    Metropolis_Hastings(G, n_accepted, coefEdges, coefTriangle, p_large_step, lambda, verbose);
    double edges = count_edges_cpp(G);
    double triangle = count_triangle_cpp(G);

    // Store the stats
    all_stats(i, 0) = edges;
    all_stats(i, 1) = triangle;
    if (verbose >= 2) {
      Rcpp::Rcout << "edges:" << edges << "\n";
      Rcpp::Rcout << "triangle:" << triangle << "\n";
    }
  }

  // Acceptance rate
  if (verbose >= 1) {
    double acceptance_rate = (n_accepted / total_iter) * 100;
    Rcpp::Rcout << "Acceptance rate: " << acceptance_rate << " %." << "\n";
  }

  // Return the output
  return all_stats;
}
