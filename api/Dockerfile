FROM rocker/shiny-verse:4.0.0

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

RUN R -e "devtools::install_github('rnabioco/valr')"

RUN R -e "install.packages('tidyr', repos = 'https://cloud.r-project.org')"

RUN R -e "remotes::install_version('Cubist', '0.3.0', upgrade = 'always', repos = 'https://cloud.r-project.org')"

RUN R -e "remotes::install_cran('patchwork', repos = 'https://cloud.r-project.org')"


RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-12.tar.gz', repos=NULL, type='source')"


RUN R -e "remotes::install_version('rstanarm', '2.21.1', repos = 'https://cran.wu.ac.at/')"

RUN R -e "devtools::install_cran('plumber')"

# Install CRAN packages
RUN R -e "install.packages(c('shiny', \
'DT', \
'yardstick', \
'parsnip', \
'glue', \
'patchwork', \
'shinycssloaders'), repos = 'https://cloud.r-project.org')"

RUN mkdir /cnvscore

COPY bancco_app/models/bayesian_clinvar_del_nohuman /cnvscore
COPY api/load_data.R /cnvscore
COPY api/main_api.R /cnvscore
COPY api/plumber.R /cnvscore

COPY bancco_app/local_data.RData.gz /cnvscore
COPY bancco_app/local_data_features1.RData.gz /cnvscore
COPY bancco_app/local_data_features2.RData /cnvscore
COPY bancco_app/local_data_features3.RData /cnvscore
COPY bancco_app/local_data_features4.RData.gz /cnvscore

COPY R/cnvscore_functions.R /cnvscore

EXPOSE 3838

ENTRYPOINT ["R", "-e", "getwd() ; list.files('/cnvscore') ; plumber::plumb('/cnvscore/main_api.R')"]
