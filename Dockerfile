FROM rocker/verse:4.0.5

RUN echo "options(repos = c(CRAN = 'https://packagemanager.rstudio.com/all/__linux__/centos7/latest'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site

# Installing the package 'myergm'.
COPY ./myergm myergm
RUN R CMD INSTALL --no-multiarch --with-keep.source myergm
