#Base image
FROM rocker/shiny:4.4.0
LABEL authors="Lis Arend"

# system libraries of general use
## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    curl \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libglpk-dev \
    libssl-dev \
    libnetcdf-dev \
    cmake \
    libmpfr-dev

## update system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get clean
    
# Adjust compiler flags globally for R
RUN echo "CFLAGS += -Wno-format-security" > /usr/local/lib/R/etc/Makevars.site
RUN echo "CXXFLAGS += -Wno-format-security" >> /usr/local/lib/R/etc/Makevars.site

## copy shiny app to work dir
RUN mkdir /srv/PRONE_app
WORKDIR /srv/PRONE_app

# install R packages via renv
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

COPY renv.lock renv.lock
RUN R -e "renv::restore()"

RUN R -e "BiocManager::install('preprocessCore', configure.args = c(preprocessCore='--disable-threading'), force=TRUE, update=TRUE, type='source')"

# copy shiny app
COPY ./*.R ./
COPY ./data ./data
COPY ./ui ./ui
COPY ./server ./server
COPY ./www ./www
COPY ./www/PRONE_Workflow.png ./www/PRONE_Workflow.png


RUN ls -la /srv/PRONE_app/www

# expose port
EXPOSE 3838

# run app on container start
CMD ["R", "-e", "shiny::runApp(host = '0.0.0.0', port = 3838)"]


