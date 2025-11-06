FROM continuumio/miniconda3
# the folllowing ARG conda_env variable must match the conda env name defined in environment.yml:
ARG conda_env=conda_rna_seq

LABEL authors="Guillaume Noell" \
  maintainer="Guillaume Noell <gn5@sanger.ak>" \
  description="Docker image for WSTI-HGI RNA-Seq Nextflow pipeline"

# nuke cache dirs before installing pkgs; tip from Dirk E fixes broken img
RUN rm -f /var/lib/dpkg/available && rm -rf  /var/cache/apt/*
RUN apt-get update && \
  apt-get -y upgrade && \
  apt-get install -y --no-install-recommends \
  r-base \
  python2.7 \
  build-essential curl git python-pip procps \ 
  g++ gcc gfortran make autoconf automake libtool \
  zlib1g-dev liblzma-dev libbz2-dev lbzip2 libgsl-dev \
  libblas-dev libx11-dev \
  libreadline-dev libxt-dev libpcre2-dev libcurl4-openssl-dev \
  && rm -rf /var/lib/apt/lists/*

RUN which R && R --version
RUN Rscript -e ".libPaths()"
RUN which python3 && python3 --version
RUN which python2 && python2 --version

# update conda && install Conda env:
RUN conda update -n base -c defaults conda
ADD environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml

# Set installed Conda env as default:
ENV CONDA_DEFAULT_ENV $conda_env
ENV PATH /opt/conda/envs/$conda_env/bin:$PATH
RUN echo $PATH

## Add additional software using pip from Conda env:
## RUN /bin/bash -c "source activate $conda_env \
##    && pip install cellSNP \
##    && pip install vireoSNP \
##    && conda env list"

## Add software that is not available via conda:

## featureCounts binary:
RUN wget https://sourceforge.net/projects/subread/files/subread-2.0.2/subread-2.0.2-Linux-x86_64.tar.gz \
  && tar -zxvf subread-2.0.2-Linux-x86_64.tar.gz \
  && cp -r subread-2.0.2-Linux-x86_64/bin/* /opt/conda/envs/$conda_env/bin/
  

## QTLtools MBV:
RUN wget https://qtltools.github.io/qtltools/binaries/QTLtools_1.2_Ubuntu16.04_x86_64.tar.gz \ 
    && tar -zxvf QTLtools_1.2_Ubuntu16.04_x86_64.tar.gz \
    && chmod a+x QTLtools_1.2_Ubuntu16.04_x86_64/QTLtools_1.2_Ubuntu16.04_x86_64 \
    && mv QTLtools_1.2_Ubuntu16.04_x86_64/* /opt/conda/envs/$conda_env/bin/ \
    && ln -s /opt/conda/envs/$conda_env/bin/QTLtools_1.2_Ubuntu16.04_x86_64 /opt/conda/envs/$conda_env/bin/QTLtools

## leafcutter regtools:
## cf. http://davidaknowles.github.io/leafcutter/articles/Installation.html
RUN git clone https://github.com/davidaknowles/leafcutter /opt/conda/envs/$conda_env/bin/leafcutter
ENV PATH /opt/conda/envs/$conda_env/bin:/opt/conda/envs/$conda_env/bin/leafcutter/scripts:/opt/conda/envs/$conda_env/bin/leafcutter/clustering:$PATH
RUN echo $PATH

# check bin  folder:
RUN ls -ltra /opt/conda/envs/$conda_env/bin

# test R libraries can  be loaded:
RUN Rscript -e "sessionInfo();.libPaths();library(AnnotationHub);library(ensembldb);library(tximport);library(magrittr);library(readr)"

# test python libraries can be loaded:
RUN python -c 'import sys;print(sys.version_info)'

## check software versions:
RUN featureCounts -v >> /usr/conda_software_versions.txt 2>&1
RUN fastqc -version >> /usr/conda_software_versions.txt
RUN salmon --version >> /usr/conda_software_versions.txt
RUN echo STAR >> /usr/conda_software_versions.txt && STAR --version >> /usr/conda_software_versions.txt
RUN multiqc --version >> /usr/conda_software_versions.txt
# now picking R and python from conda env by default:
RUN which R >> /usr/conda_software_versions.txt && R --version >> /usr/conda_software_versions.txt
RUN Rscript -e ".libPaths()" >> /usr/conda_software_versions.txt
RUN which python >> /usr/conda_software_versions.txt && python --version >> /usr/conda_software_versions.txt
RUN cat /usr/conda_software_versions.txt

CMD /bin/sh
