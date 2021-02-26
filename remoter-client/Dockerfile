# Create builder step
# pull official base image and use it as builder
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
ENV SERVER_HOSTNAME $SERVER_HOSTNAME:6969

# install pak for more streamlined installation of packages
RUN R -e 'utils::install.packages("pak", repos = "https://r-lib.github.io/p/pak/dev/")'

# install required packages
RUN R -e 'pak::pkg_install("remoter")'

# add app
COPY src ./

# ---------------------------------------------------
# PRODUCTION BUILD
# ---------------------------------------------------
FROM builder AS prod

# start app
ENTRYPOINT ["bash", "/var/lib/watchfile/entrypoint.sh"]
CMD ["R", "--no-save", "-s", "-e", "require('remoter');remoter::batch(addr='remoter-server', port=6969, file='rscript.r')" ]


# ---------------------------------------------------
# DEVELOPMENT BUILD
# ---------------------------------------------------
FROM builder AS dev

# install Radian for interactive R shell
# also install watchdog to automatically restart
# when source files change
RUN pip install -U jedi radian PyYAML watchdog[watchmedo]

# start app
CMD ["R", "--no-save", "-s", "-e", "require('remoter');remoter::batch(addr='remoter-server', port=6969, file='rscript.r')" ]
