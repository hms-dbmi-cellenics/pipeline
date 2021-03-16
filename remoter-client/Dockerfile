FROM r-base:4.0.2 AS builder
WORKDIR /src

# install required debian packages for R installation
ENV DEBIAN_FRONTEND noninteractive
COPY requirements_debian.txt .
RUN apt-get update && apt-get -y upgrade
RUN cat requirements_debian.txt | xargs apt-get -y install --no-install-recommends
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# add GitHub PAT if required for GitHub installations.
ARG GITHUB_PAT
ENV GITHUB_PAT ${GITHUB_PAT}

# install pak for more streamlined installation of packages
RUN R -e 'utils::install.packages("pak", repos = "https://r-lib.github.io/p/pak/dev/")'

# install required packages
RUN R -e 'pak::pkg_install("remoter")'
RUN R -e 'pak::pkg_install("ids")'
RUN R -e 'pak::pkg_install("import")'
RUN R -e 'pak::pkg_install("RJSONIO")'

# ---------------------------------------------------
# PRODUCTION BUILD
# ---------------------------------------------------
FROM builder AS prod

# start app
COPY src ./

ENTRYPOINT ["Rscript", "init.r"]
