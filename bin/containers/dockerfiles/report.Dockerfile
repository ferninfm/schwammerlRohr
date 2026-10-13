FROM rocker/r2u:24.04
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
RUN apt-get update && apt-get install -y --no-install-recommends \
    pandoc git ca-certificates libxml2-dev libcurl4-openssl-dev libssl-dev libgit2-dev \
    && rm -rf /var/lib/apt/lists/*
RUN R --quiet --vanilla -e 'if (requireNamespace("bspm", quietly=TRUE)) bspm::disable(); options(bspm.enable=FALSE, repos=c(CRAN="https://cloud.r-project.org")); install.packages("pak", repos="https://cloud.r-project.org")'
RUN R --quiet --vanilla -e 'pak::pak(c("ggplot2", "dplyr", "tidyr", "jsonlite", "knitr", "kableExtra", "purrr"), upgrade=FALSE)'
RUN R --quiet --vanilla -e 'pak::pak(c("taui11/mycoBinR"), upgrade=FALSE)'
