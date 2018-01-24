#!/usr/bin/env bash

# Install miniconda

if [[ -d $HOME/miniconda3/bin ]]; then
    echo "miniconda3 already installed."
else
    echo "Installing miniconda3"
    mkdir -p $HOME/download
    url="https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    wget \
        --continue \
        --output-document $HOME/download/miniconda.sh \
        $url
    chmod +x $HOME/download/miniconda.sh
    $HOME/download/miniconda.sh \
        -u \
        -b \
        -p $HOME/miniconda3

    $HOME/miniconda3/bin/conda clean --all --yes
fi
