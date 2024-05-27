#syntax=docker/dockerfile:1

FROM rocker/r-base:4.4.0

# github setup
ARG GITHUB_PAT

RUN apt-get update && apt-get install -y \
    binutils \
    sudo \
    pandoc \
    apt-utils \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    cmake \
    libharfbuzz-dev \
    libfribidi-dev

# last three are for scBio/CPM

# omnideconv and DeconvExplorer
#RUN install2.r --error textshaping

RUN R -e "install.packages('pak')"
RUN R -e "pak::pkg_install(c('car', 'matlib', 'remotes'))"
RUN R -e "remotes::install_github('omnideconv/omnideconv', dependencies=TRUE)"
RUN R -e "pak::pkg_install('omnideconv/DeconvExplorer', dependencies=TRUE)"

#RUN R -e "reticulate::install_miniconda()" # will not install automatically


RUN R -e "omnideconv::install_all_python()"

#RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))" > /usr/lib/R/etc/Rprofile.site

# in case set python environment for reticulate
# Sys.setenv...

# Cibersort Source Code
COPY CIBERSORT.R .
#RUN R -e "omnideconv::bseqsc_config(file = 'CIBERSORT.R')"

RUN groupadd  shiny \
&& useradd --gid shiny --shell /bin/bash --create-home shiny

# copy server configs
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY shiny-server.sh /usr/bin/shiny-server.sh
COPY /*.R /srv/shiny-server/
#RUN chmod -R +r /srv/shiny-server
RUN chown -R shiny:shiny /srv/shiny-server

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/', host = '0.0.0.0', port = 3838)"]
#CMD ["/usr/bin/shiny-server.sh"]
