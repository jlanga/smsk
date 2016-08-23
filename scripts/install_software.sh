#!/usr/bin/env bash

## 1. Pip packages
pip install \
    snakemake \
    pyyaml \
    docutils



## 2. Homebrew packages
brew update

brew install \
    graphviz

brew tap homebrew/science
brew install \
    homebrew/science/bwa \
    homebrew/science/samtools \
    homebrew/science/bcftools



## 3. Custom packages
## wget tar.gz; tar xvf tar.gz; make; make test; make install
