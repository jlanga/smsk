# smsk: A Snakemake skeleton to jumpstart projects

## 1. Description

This is a small skeleton to create Snakemake workflows. [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) "is a workflow management system that aims to reduce the complexity of creating workflows by providing a fast and comfortable execution environment, together with a clean and modern specification language in python style."

The idea is to create a workflow with of snakefiles, resolve dependencies with pip and brew, and

## 2. First steps

1. Clone this repo:

    ```sh
    git clone https://github.com/jlanga/smsk.git my_project
    cd my_project
    ```

2. Create a Python(3) environment:
    ```sh
    virtualenv --python=python3 bin/py3
    ```


3. Execute `scripts/install_brew.sh` to download homebrew:

    ```sh
    bash scripts/install_brew.sh
    ```

4. Activate the environment:
    ```sh
    source bin/py3/bin/activate
    ```

5. Install software and packages via pip and homebrew :

    ```sh
    bash scripts/install_software.sh
    ```
6. Execute the test files:

    ```{sh}
    snakemake \
        --cores 8
    ```



## 3. File organization

The hierarchy of the folder is the one described in [A Quick Guide to Organizing Computational Biology Projects](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424):

```
smsk
├── .linuxbrew: brew files
├── bin: your binaries and virtualenv related files.
├── data: raw data, hopefully links to backuped data.
├── doc: reports and figures.
├── README.md
├── results: processed data.
├── scripts: python, R, etc scripts to porcess data.
└── src: additional source code, tarballs, etc.
```



## 4. Configuration file



## 5. Considerations when installing software

As a rule of thumb, download python packages with `pip`, use `brew` whenever possible, download binary tarballs into `src/`` and copy them to `bin/` or download the source tarball and compile it. Example:

   ```
   pip install \
       snakemake

   brew install \
       samtools

   wget \
       --continue \
       --output-document src/bin1.tar.gz \
       http://bin1.com/bin1.tar.gz
   tar xvf src/bin1.tar.gz
   cp src/bin1/bin1 bin/ # or link

   wget \
       --continue \
       --output-document src/bin2.tar.gz \
       http://bin2.com/bin2.tar.gz
   tar xvf src/bin2.tar.gz
   pushd src/bin2/
   make -j
   cp build/bin2 ../../bin/
   ```

## Bibliography

- [A Quick Guide to Organizing Computational Biology Projects](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)

- [Snakemake—a scalable bioinformatics workflow engine](http://bioinformatics.oxfordjournals.org/content/28/19/2520)
