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
    ++j;
  }
  adjmat_next(i, j) = 1 - adjmat_next(i, j);
  adjmat_next(j, i) = 1 - adjmat_next(j, i);

  return adjmat_next;
}


// [[Rcpp::export]]
arma::sp_mat change_one_link_modified(const arma::sp_mat& adjmat, double& numOfEdges, double& numOfTriangles, int i, int j, int verbose) {
  arma::sp_mat G = adjmat;

  // Change the state of the selected dyad.
  G(i, j) = 1 - G(i, j);
  G(j, i) = 1 - G(j, i);

  // If added a link, add one to numOfEdges. Otherwise, subtract one.
  if (G(i, j) == 1) {
    numOfEdges += 1;
    if (verbose >= 5) {
      Rcpp::Rcout << "Added link: (" << i << ", " << j << ")"<<"\n";
    }
  } else {
    numOfEdges += -1;
    if (verbose >= 5) {
      Rcpp::Rcout << "Deleted link: (" << i << ", " << j << ")"<<"\n";
    }
  }

  // Count the number of triangles that appear/disappear by adding/removing the link.
  arma::sp_mat col_i = G.col(i);
  arma::sp_mat col_j = G.col(j);
  arma::sp_mat ij = col_i.t() * col_j;
  if (G(i, j) == 1) {
    numOfTriangles += ij(0, 0);
  } else {
    numOfTriangles += -ij(0, 0);
  }

  return G;
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
      ++j;
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
                         double& numOfEdges,
                         double& numOfTriangles,
                         double numOfNodes,
                         double coefEdges,
                         double coefTriangle,
                         double p_large_step,
                         double lambda,
                         int verbose) {
  // Draw a number from (0, 1).
  double x = unif_rand();

  // Initialization
  double numOfEdges_next = numOfEdges;
  double numOfTriangles_next = numOfTriangles;
  double changestat_edges = numOfEdges_next - numOfEdges;
  double changestat_triangle = numOfTriangles_next - numOfTriangles;
  double logratio = 0;

  // If x < p_large_step, then select one node and flip its all links.
  if (x < p_large_step) {
    if (verbose >= 5) {
      Rcpp::Rcout << "Flipping all the links of one node..." << "\n";
      }
    // Update all the links of one node.
    arma::sp_mat adjmat_next = change_all_links_of_one_node(adjmat);

    // Store the results.
    numOfEdges_next = count_edges_cpp(adjmat_next);
    numOfTriangles_next = count_triangle_cpp(adjmat_next);
    changestat_edges = numOfEdges_next - numOfEdges;
    changestat_triangle = numOfTriangles_next - numOfTriangles;

    // Compute the logratio
    logratio = changestat_edges * coefEdges + changestat_triangle * coefTriangle;

    // Accept or reject?
    if (logratio >= 0.0 || log(unif_rand()) < logratio) {
      if (verbose >= 5) {
        Rcpp::Rcout << "Accepted." << "\n";
      }
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

  // If p_large_step <= x < 2 * p_large_step, then flip [lambda * n] links simultaneously.
  else if (x >= p_large_step && x < 2 * p_large_step) {
    if (verbose >= 5) {
      Rcpp::Rcout << "Flipping lambda * N links simultaneously..." << "\n";
    }
    // Update [lambda * n] links.
    arma::sp_mat adjmat_next = change_multiple_links(adjmat, lambda);

    // Store the results.
    numOfEdges_next = count_edges_cpp(adjmat_next);
    numOfTriangles_next = count_triangle_cpp(adjmat_next);
    changestat_edges = numOfEdges_next - numOfEdges;
    changestat_triangle = numOfTriangles_next - numOfTriangles;

    // Compute the logratio
    logratio = changestat_edges * coefEdges + changestat_triangle * coefTriangle;

    // Accept or reject?
    if (logratio >= 0.0 || log(unif_rand()) < logratio) {
      if (verbose >= 5) {
        Rcpp::Rcout << "Accepted." << "\n";
      }
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

  // Flip one link randomply.
  else {
    if (verbose >= 5) {
      Rcpp::Rcout << "Flipping one link..." << "\n";
    }
    // Choose one dyad randomly and change the state.
    // pick a random element
    int i = unif_rand() * numOfNodes;
    // pick a random element from what's left (there is one fewer to choose from)...
    int j =  unif_rand() * (numOfNodes- 1);
    // ...and adjust second choice to take into account the first choice
    if (j >= i)
    {
      ++j;
    }
    // Store the results.
    arma::sp_mat adjmat_next = change_one_link_modified(adjmat, numOfEdges_next, numOfTriangles_next, i, j, verbose);

    // Store the results.
    changestat_edges = numOfEdges_next - numOfEdges;
    changestat_triangle = numOfTriangles_next - numOfTriangles;

    // Compute the logratio
    logratio = changestat_edges * coefEdges + changestat_triangle * coefTriangle;

        // Accept or reject?
    if (logratio >= 0.0 || log(unif_rand()) < logratio) {
      if (verbose >= 5) {
        Rcpp::Rcout << "Accepted." << "\n";
      }
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
}


// Main wrapper function for MCMC.
// [[Rcpp::export]]
arma::mat create_MCMC(const arma::sp_mat& adjmat,
                      double coefEdges,
                      double coefTriangle,
                      int MCMC_interval = 1024,
                      int MCMC_samplesize = 1024,
                      int MCMC_burnin = 1024 * 16,
                      double p_large_step = 0.01,
                      double lambda = 0.5,
                      bool full_sample = false,
                      int verbose = 0) {
  // Total number of iterations
  int total_iter = MCMC_burnin + MCMC_interval * MCMC_samplesize;

  int mat_size = MCMC_samplesize;
  if (full_sample == true) {
    mat_size = total_iter;
  }

  // Create a matrix in which network stats during MCMC are stored.
  arma::mat stats = arma::zeros(mat_size , 2);

  // Initialization
  arma::sp_mat G = adjmat;
  double n_accepted = 0;
  double numOfEdges = count_edges_cpp(adjmat);
  double numOfTriangles = count_triangle_cpp(adjmat);
  double numOfNodes = adjmat.n_rows;

  // Start MCMC
  if (full_sample == true) {
    for (int i = 0; i < total_iter; i++) {
      Metropolis_Hastings(G, n_accepted, numOfEdges, numOfTriangles, numOfNodes, coefEdges, coefTriangle, p_large_step, lambda, verbose);

      // Store the stats
      stats(i, 0) = numOfEdges;
      stats(i, 1) = numOfTriangles;

      // Check the current state if necessary
      if (verbose >= 2) {
        Rcpp::Rcout << "edges:" << numOfEdges << ", triangle:" << numOfTriangles << "\n";
      }
    }
  }
  else {
    for (int i = 0; i < total_iter; i++) {
      Metropolis_Hastings(G, n_accepted, numOfEdges, numOfTriangles, numOfNodes, coefEdges, coefTriangle, p_large_step, lambda, verbose);

      // Store the stats
      if (i >= MCMC_burnin && (i - MCMC_burnin) % 1024 == 0) {
        int index = (i - MCMC_burnin) / 1024;
        stats(index, 0) = numOfEdges;
        stats(index, 1) = numOfTriangles;
      }
      // Check the current state if necessary
      if (verbose >= 2) {
        Rcpp::Rcout << "edges:" << numOfEdges << ", triangle:" << numOfTriangles << "\n";
      }
    }
  }

  // Acceptance rate
  if (verbose >= 1) {
    double acceptance_rate = (n_accepted / total_iter) * 100;
    Rcpp::Rcout << "Acceptance rate: " << acceptance_rate << " %." << "\n";
  }

  // Return the output
  return stats;
}
