FROM rocker/verse:3.6.3

RUN R -e "install.packages(c('Rcpp', 'RcppArmadillo', 'statnet', 'igraph', 'intergraph'), repos = c(CRAN = 'https://cloud.r-project.org'))"
