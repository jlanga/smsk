#!/usr/bin/env bash

set -euo pipefail

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda clean --all --yes
