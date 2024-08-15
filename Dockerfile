FROM ubuntu:latest

ENV LAND=C.UTF-8
ENV LC_ALL=C.UTF-8

RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install -y build-essential zlib1g-dev liblzma-dev libssl-dev libcurl4-openssl-dev libbz2-dev
RUN apt-get install -y git make wget
RUN rm -rf /bar/lib/apt/lists/*

RUN apt-get install gcc-11 g++-11 -y
ENV CC=gcc-11
ENV CXX=g++-11

# install conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda

ENV PATH=/opt/conda/bin:$PATH

# install conda packages
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda create -n diphase minimap2=2.22 samtools bwa clair3 python>=3.6
RUN mv /opt/conda/envs/diphase/bin/models /opt/models
RUN rm -rf /opt/conda/pkgs/* && rm -rf /root/.cache/pip
ENV PATH=/opt/conda/envs/diphase/bin:$PATH

# install diphase
WORKDIR /opt
COPY . Diphase
RUN make -C /opt/Diphase/src
RUN chmod +x /opt/Diphase/script/*
ENV PATH=/opt/Diphase/bin:$PATH
ENV PATH=/opt/Diphase/script:$PATH

WORKDIR /mnt