FROM rocker/verse:4.0.5

# Installing necessary packages.
RUN Rscript -e "install.packages(c('renv', 'devtools'), repos = 'https://cran.rstudio.com')"
WORKDIR ./
COPY ./myergm/renv.lock renv.lock
RUN Rscript -e "renv::restore(repos = 'https://cran.rstudio.com')"

# Installing the package 'myergm'.
COPY ./myergm myergm
RUN R CMD INSTALL --no-multiarch --with-keep.source myergm
