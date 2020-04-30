FROM ubuntu:18.04

SHELL ["/bin/bash", "--login", "-c"]

# apt packages
RUN apt update && apt install --yes \
    wget \
&& rm -rf /var/lib/apt/lists/*


ENV miniconda=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/opt/miniconda3/bin:$PATH"

RUN \
    wget --quiet --continue $miniconda && \
    bash Miniconda3-latest-Linux-x86_64.sh -bfp /opt/miniconda3 && \
    conda clean --all --yes && \
    rm Miniconda3-latest-Linux-x86_64.sh

RUN \
    conda config --add channels defaults && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda install --yes --channel bioconda snakemake-minimal=5.15 pandas && \
    conda clean --all --yes

RUN mkdir /.conda && chmod ugo+rwx /.conda
