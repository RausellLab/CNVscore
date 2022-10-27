FROM rocker/shiny-verse:4.2.1

#libv8-3.14-dev \


# install linux libraries
RUN apt-get update && apt-get install libcurl4-openssl-dev \
libbz2-dev \
zlib1g-dev \
libv8-dev \
libbz2-1.0 \
libbz2-ocaml \
libbz2-ocaml-dev \
libjpeg-dev \
libncurses5-dev \
libncursesw5-dev \
liblzma-dev  -y &&\
mkdir -p /var/lib/shiny-server/bookmarks/shiny

RUN apt-get install libnlopt-dev -y

RUN apt-get install cmake -y
# Install GitHub packages

RUN R -e "devtools::install_github('egenn/rtemis@c0148033c6c6557cd7b7efc3dc4be002891b282f')"

RUN R -e "devtools::install_github('rnabioco/valr@v0.6.3')"

RUN R -e "install.packages('tidyr', repos = 'https://cloud.r-project.org')"

RUN R -e "remotes::install_version('Cubist', '0.3.0', upgrade = 'always', repos = 'https://cloud.r-project.org')"

RUN R -e "remotes::install_cran('patchwork', repos = 'https://cloud.r-project.org')"


RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-12.tar.gz', repos=NULL, type='source')"

RUN R -e "devtools::install_version('xgboost', version = '1.4.1.1', repos = 'http://cran.us.r-project.org')"
RUN R -e "devtools::install_version('parsnip', version = '0.1.7', repos = 'http://cran.us.r-project.org')"

# RUN R -e "remotes::install_version('rstanarm', '2.21.1', repos = 'https://cran.wu.ac.at/')"


# Install CRAN packages
RUN R -e "install.packages(c('shiny', \
'DT', \
'yardstick', \
'glue', \
'patchwork', \
'shinycssloaders'), repos = 'https://cloud.r-project.org')"


# Copy the app to the image
COPY bancco_app /srv/shiny-server/

# test
RUN cat /srv/shiny-server/shiny-server.txt > /etc/shiny-server/shiny-server.conf

# Unzip local data
RUN echo 'allow_app_override;' >> /etc/shiny-server/shiny-server.conf

# Make all app files readable
RUN chmod -R +r /srv/shiny-server/

EXPOSE 3838

CMD ["/init"]

#CMD ["/usr/bin/shiny-server.sh"] 
