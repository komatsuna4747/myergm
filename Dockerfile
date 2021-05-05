FROM rocker/verse:3.6.3

# Installing necessary packages.
RUN R -e "install.packages(c('Rcpp', 'RcppArmadillo', 'statnet', 'igraph', 'intergraph'), repos = c(CRAN = 'https://cloud.r-project.org'))"

# Installing the package 'myergm'.
RUN R -e "devtools::install_github('komatsuna4747/myergm@main', subdir = 'myergm')"
