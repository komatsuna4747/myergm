#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Auxilirary functions
// [[Rcpp::export]]
double count_edges(const arma::mat& adjmat) {
  double edges = accu(adjmat);
  return edges/2;
}


// [[Rcpp::export]]
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


// [[Rcpp::export]]
arma::mat update_one_link(const arma::mat& adjmat, double& numOfEdges, double& numOfTriangles, int i, int j, int verbose) {
  arma::mat G = adjmat;

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
  arma::vec col_i = G.col(i);
  arma::vec col_j = G.col(j);
  double ij = dot(col_i, col_j);
  if (G(i, j) == 1) {
    numOfTriangles += ij;
  } else {
    numOfTriangles += -ij;
  }

  return G;
}


// [[Rcpp::export]]
arma::mat update_all_links_of_one_node(const arma::mat& adjmat, int i) {
  arma::mat adjmat_next = adjmat;

  // Change all links of one node
  arma::vec G_i = adjmat.col(i);
  G_i = 1 - G_i;
  adjmat_next.col(i) = G_i;
  adjmat_next.row(i) = G_i.t();

  // Set diagonal entries to be zero
  adjmat_next.diag().zeros();

  // Return the output
  return adjmat_next;
}
