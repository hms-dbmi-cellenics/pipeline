# Create builder step
# pull official base image and use it as builder
FROM rocker/r-ver:4.4.0 AS builder
WORKDIR /src/pipeline-runner

# install required debian packages to install R packages
COPY setup/install_debian_packages.sh .
COPY setup/sysdeps_build_debian.txt .
RUN cat sysdeps_build_debian.txt | xargs ./install_debian_packages.sh

# add GitHub PAT if required for GitHub installations.
ARG GITHUB_PAT
ENV GITHUB_PAT $GITHUB_PAT
RUN R -q -e "if(Sys.getenv('GITHUB_PAT') == '') stop('need to export GITHUB_PAT')"

# add renv library to .libPaths
# this is path with R package folders when renv deactivated or not used
ENV RENV_LIB=/src/lib
RUN echo ".libPaths(c('$RENV_LIB', .libPaths()))" >> $(R RHOME)/etc/Rprofile.site

# install renv to install required R packages
RUN R -q -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    R -q -e "remotes::install_github('rstudio/renv@v1.0.7')" && \
    R -q -e "renv::init(bare = TRUE, settings = list(use.cache = FALSE))"

# fast pre-restore with pkgbuild
# an initial lockfile is used to avoid frequent re-installs
# copy renv.lock.init to renv.lock if R dependency updates become slow to build
COPY ./renv.lock.init ./renv.lock
COPY setup/restore_fast.R .
RUN Rscript restore_fast.R

# restore anything pkgbuild missed
COPY setup/restore_renv.R .
RUN Rscript restore_renv.R

# deactivate so that .libPaths used (set above)
RUN R -q -e "renv::deactivate()"

# use renv::snapshot() while R dependency updates are quick to build
COPY ./renv.lock .
RUN Rscript restore_renv.R

# determine system run-time deps
COPY setup/get_sysdeps_run.R .
RUN Rscript get_sysdeps_run.R

# remove identified packages with potentially conflicting licenses
COPY setup/create_ghost_packages.R .
RUN Rscript create_ghost_packages.R

# check for any potentially problematic licenses
COPY setup/check_package_licenses.R .
RUN Rscript check_package_licenses.R

# ---------------------------------------------------
# COMMON MINIMAL BUILD
# ---------------------------------------------------
FROM rocker/r-ver:4.4.0 AS common
WORKDIR /src/pipeline-runner

# get source code and R packages
COPY --from=builder /src /src

# add renv library to .libPaths
ENV RENV_LIB=/src/lib
RUN echo ".libPaths(c('$RENV_LIB', .libPaths()))" >> $(R RHOME)/etc/Rprofile.site && \
    echo "python3-pip python3-venv python3-dev wget" >> sysdeps_run.txt && \
    cat sysdeps_run.txt | xargs ./install_debian_packages.sh && \
    rm -rf *

# install python packages in virtualenv
ENV WORKON_HOME=/src/.virtualenvs
RUN R -q -e "reticulate::virtualenv_create('r-reticulate')" && \
    R -q -e "reticulate::virtualenv_install('r-reticulate', c('geosketch==1.2', 'scanorama==1.7.3'), pip_options='--no-cache-dir')"

# ---------------------------------------------------
# PRODUCTION BUILD
# ---------------------------------------------------
FROM common AS prod

# add jq and aws cli so we can get identities using AWS IRSA roles
#
# see at https://docs.aws.amazon.com/eks/latest/userguide/iam-roles-for-service-accounts.html
# for what IRSA roles are
#
# see at https://github.com/paws-r/paws/blob/main/docs/credentials.md for how PAWS
# handles this
RUN wget -O jq https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64 && \
    chmod +x ./jq && \
    cp jq /usr/bin && \
    pip install PyYAML==3.10 awscli==1.25.0

COPY aws_config /root/.aws/config

# add R package files and runner
ADD R ./R
ADD tests ./tests
ADD inst ./inst
COPY DESCRIPTION NAMESPACE init.R ./

# start app
ENTRYPOINT ["Rscript", "init.R"]

# ---------------------------------------------------
# AWS BATCH BUILD
# ---------------------------------------------------
FROM common AS batch

# add jq and aws cli so we can get identities using AWS IRSA roles
#
# see at https://docs.aws.amazon.com/eks/latest/userguide/iam-roles-for-service-accounts.html
# for what IRSA roles are
#
# see at https://github.com/paws-r/paws/blob/main/docs/credentials.md for how PAWS
# handles this
COPY ./datadog-batch ./datadog-batch
RUN wget -O jq https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64 && \
    chmod +x ./jq && \
    cp jq /usr/bin && \
    pip install PyYAML==3.10 awscli==1.25.0 -r ./datadog-batch/requirements.txt

# add R package files and runner
ADD R ./R
ADD tests ./tests
ADD inst ./inst
COPY DESCRIPTION NAMESPACE init.R start-batch.sh ./


# ---------------------------------------------------
# DEVELOPMENT BUILD
# ---------------------------------------------------
FROM common AS dev

# install watchdog to automatically restart
# when source files change
RUN pip install -U jedi PyYAML watchdog[watchmedo] && \
    apt update && \
    apt -y install git && \
    rm -rf /var/cache/apt/archives /var/lib/apt/lists/*

# add R package files and runner
ADD R ./R
ADD tests ./tests
ADD inst ./inst
COPY DESCRIPTION NAMESPACE init.R ./

# start app
ENTRYPOINT Rscript init.R
