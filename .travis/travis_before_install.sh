if test -e $HOME/miniconda3/bin ; then
    echo "miniconda already installed."
else
    echo "Installing miniconda."
    rm -rf $HOME/miniconda3
    mkdir -p $HOME/download
    if [[ -d $HOME/download/miniconda.sh ]] ; then
        rm -rf $HOME/download/miniconda.sh
    fi
    url="https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    wget \
        --continue \
        --output-document $HOME/download/miniconda.sh \
        $url
    chmod +x $HOME/download/miniconda.sh
    $HOME/download/miniconda.sh -b
fi
