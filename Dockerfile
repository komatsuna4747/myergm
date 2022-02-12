FROM rocker/rstudio:4.0.5

# Installing necessary packages.
RUN echo "options(repos = c(CRAN = 'https://packagemanager.rstudio.com/all/__linux__/centos7/latest'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site
RUN Rscript -e "install.packages('renv')"

RUN mkdir -p /home/rstudio/.local/share/renv/cache \
    && chown -R rstudio:rstudio /home/rstudio
