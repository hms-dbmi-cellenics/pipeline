# Create builder step
# pull official base image and use it as builder
FROM r-base:4.0.5 AS builder
WORKDIR /src

# install required debian packages for R installation
ENV DEBIAN_FRONTEND noninteractive
COPY requirements_debian.txt .
RUN apt-get update && \
    apt-mark hold r-base-core r-base-dev && \
    apt-get -y upgrade
RUN cat requirements_debian.txt | xargs apt-get -y install --no-install-recommends

# default BLAS/LAPACK issue: Seurat::FindIntegrationAnchors finds ~30 vs ~13,000 anchors
RUN update-alternatives --set libblas.so.3-x86_64-linux-gnu /usr/lib/x86_64-linux-gnu/atlas/libblas.so.3

RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# add GitHub PAT if required for GitHub installations.
ARG GITHUB_PAT
ENV GITHUB_PAT $GITHUB_PAT
RUN R -e "if(Sys.getenv('GITHUB_PAT') == '') stop('need to export GITHUB_PAT')"

# install required R packages for the pipeline.
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    R -e "remotes::install_github('rstudio/renv@0.13.2')"

# an initial lockfile is used to avoid frequent re-installs
# use renv::snapshot(lockfile='renv.lock.init') if R dependency updates become slow to build
COPY renv.lock.init .
RUN R -e "renv::init(bare = TRUE); options(renv.consent = TRUE); renv::restore(lockfile='renv.lock.init')"

# use renv::snapshot() while R dependency updates are quick to build
COPY renv.lock .
RUN R -e "options(renv.consent = TRUE); renv::restore(lockfile='renv.lock', clean = TRUE)"


# ---------------------------------------------------
# PRODUCTION BUILD
# ---------------------------------------------------
FROM builder AS prod

# add jq and aws cli so we can get identities using AWS IRSA roles
#
# see at https://docs.aws.amazon.com/eks/latest/userguide/iam-roles-for-service-accounts.html
# for what IRSA roles are
#
# see at https://github.com/paws-r/paws/blob/main/docs/credentials.md for how PAWS
# handles this
#
RUN wget -O jq https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64 && chmod +x ./jq && cp jq /usr/bin
RUN pip install awscli
COPY aws_config /root/.aws/config

# start app
ADD R ./R
COPY init.R .
ENTRYPOINT ["Rscript", "init.R"]


# ---------------------------------------------------
# DEVELOPMENT BUILD
# ---------------------------------------------------
FROM builder AS dev

# install Radian for interactive R shell
# also install watchdog to automatically restart
# when source files change
RUN pip install -U jedi radian PyYAML watchdog[watchmedo]

# start app
ADD R ./R
COPY init.R .
ENTRYPOINT ["Rscript", "init.R"]
