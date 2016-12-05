conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda env create --name smsk --file bin/install/requirements.txt
