FROM condaforge/miniforge3:25.9.1-0

ARG conda_env=conda_rnaseq_tximport

LABEL authors="Gennadii Zakharov" \
  maintainer="Gennadii Zakharov <gz35@sanger.ac.uk>" \
  description="TxImport Docker image for WSTI-HGI RNA-Seq Nextflow pipeline"

# nuke cache dirs before installing pkgs; tip from Dirk E fixes broken img
RUN rm -f /var/lib/dpkg/available && rm -rf  /var/cache/apt/*
RUN apt-get update && \
  apt-get -y upgrade && \
  apt-get install -y --no-install-recommends \
  r-base \
  build-essential curl git python3-pip procps \
  && rm -rf /var/lib/apt/lists/*

# g++ gcc gfortran make autoconf automake libtool \
# zlib1g-dev liblzma-dev libbz2-dev lbzip2 libgsl-dev \

# libreadline-dev libxt-dev libpcre2-dev libcurl4-openssl-dev \

RUN which R && R --version
RUN Rscript -e ".libPaths()"
RUN which python3 && python3 --version

# update conda && install Conda env:
RUN mamba update -n base
ADD environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml

# Set installed Conda env as default:
ENV CONDA_DEFAULT_ENV=$conda_env
ENV PATH=/opt/conda/envs/$conda_env/bin:$PATH
RUN echo $PATH

# check bin  folder:
RUN ls -ltra /opt/conda/envs/$conda_env/bin

# test R libraries can  be loaded:
RUN Rscript -e "sessionInfo();.libPaths();library(AnnotationHub);library(ensembldb);library(tximport);library(magrittr);library(readr)"

## check software versions:
RUN which R >> /usr/conda_software_versions.txt && R --version >> /usr/conda_software_versions.txt
RUN Rscript -e ".libPaths()" >> /usr/conda_software_versions.txt
RUN which python >> /usr/conda_software_versions.txt && python --version >> /usr/conda_software_versions.txt
RUN cat /usr/conda_software_versions.txt

CMD ["/bin/sh"]
