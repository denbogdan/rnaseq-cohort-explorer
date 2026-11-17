FROM rocker/shiny:4.4.1
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev libcairo2-dev libglpk-dev libgsl-dev \
    libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /srv/app
COPY install_packages.R .
RUN Rscript install_packages.R
COPY R/ R/
COPY app/ app/
COPY config/ config/
RUN mkdir -p cache data && chown -R shiny:shiny /srv/app
USER shiny
EXPOSE 3838
CMD ["R", "-e", "rmarkdown::run('app/app.Rmd', shiny_args = list(port = 3838, host = '0.0.0.0'))"]
