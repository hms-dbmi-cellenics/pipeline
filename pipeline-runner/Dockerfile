# Create builder step
# pull official base image and use it as builder
FROM rocker/r-ver:4.2.2 AS builder
WORKDIR /src/pipeline-runner

# install required debian packages to install R packages
COPY setup/install_debian_packages.sh .
COPY setup/sysdeps_build_debian.txt .
RUN cat sysdeps_build_debian.txt | xargs ./install_debian_packages.sh

# add GitHub PAT if required for GitHub installations.
ARG GITHUB_PAT
ENV GITHUB_PAT $GITHUB_PAT
RUN R -q -e "if(Sys.getenv('GITHUB_PAT') == '') stop('need to export GITHUB_PAT')"

# renv library
ENV RENV_LIB=/src/lib

# add renv library to .libPaths
RUN echo ".libPaths(c('$RENV_LIB', .libPaths()))" >> $(R RHOME)/etc/Rprofile.site

# install renv to install required R packages
RUN R -q -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" && \
    R -q -e "remotes::install_github('rstudio/renv@0.16.0')" && \
    R -q -e "renv::init(bare = TRUE, settings = list(use.cache = FALSE))"

# an initial lockfile is used to avoid frequent re-installs
# use renv::snapshot(lockfile='renv.lock.init') if R dependency updates become slow to build
# delete renv cache
# strip debug from shared libraries
# see http://dirk.eddelbuettel.com/blog/2017/08/20/#010_stripping_shared_libraries
COPY ./renv.lock.init .
RUN R -q -e "renv::restore(lockfile='renv.lock.init', library = '$RENV_LIB')" && \
    R -q -e 'root <- renv::paths$root(); unlink(root, recursive = TRUE)' && \
    strip --strip-debug $RENV_LIB/*/libs/*.so

RUN R -q -e "renv::deactivate()"

# install miniconda and geosketch
# clean conda
ENV RETICULATE_MINICONDA_PATH=/src/r-miniconda
RUN R -q -e "reticulate::install_miniconda()" && \
    R -q -e "reticulate::conda_install('r-reticulate', c('geosketch', 'scanorama'), pip = TRUE)" && \
    CONDA_PATH=$(R -q -s -e "cat(reticulate::conda_binary())") && \
    $CONDA_PATH clean --force-pkgs-dirs -y

# use renv::snapshot() while R dependency updates are quick to build
COPY ./renv.lock .
RUN R -q -e "renv::restore(lockfile='renv.lock', library = '$RENV_LIB', clean = TRUE)" && \
    R -q -e 'root <- renv::paths$root(); unlink(root, recursive = TRUE)' && \
    strip --strip-debug $RENV_LIB/*/libs/*.so

# determine system run-time deps
COPY setup/get_sysdeps_run.R .
RUN Rscript get_sysdeps_run.R

# ---------------------------------------------------
# COMMON MINIMAL BUILD
# ---------------------------------------------------
FROM rocker/r-ver:4.2.2 AS common
WORKDIR /src/pipeline-runner
ENV RETICULATE_MINICONDA_PATH=/src/r-miniconda

# get source code and R packages
COPY --from=builder /src /src

# add renv library to .libPaths
ENV RENV_LIB=/src/lib
RUN echo ".libPaths(c('$RENV_LIB', .libPaths()))" >> $(R RHOME)/etc/Rprofile.site && \
    echo "python3-pip" >> sysdeps_run.txt && \
    echo "wget" >> sysdeps_run.txt && \
    cat sysdeps_run.txt | xargs ./install_debian_packages.sh && \
    rm -rf *

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
