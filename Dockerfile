FROM rocker/verse:3.6.3

# Installing necessary packages.
RUN echo "options(repos = c(CRAN = 'https://packagemanager.rstudio.com/all/__linux__/centos7/latest'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site

# Installing the package 'myergm'.
RUN R -e "devtools::install_github('komatsuna4747/myergm@main', subdir = 'myergm')"
