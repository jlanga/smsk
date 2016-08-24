# smsk: A Snakemake skeleton to jumpstart projects

## 1. Description

This is a small skeleton to create Snakemake workflows. [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) "is a workflow management system that aims to reduce the complexity of creating workflows by providing a fast and comfortable execution environment, together with a clean and modern specification language in python style."

The idea is to create a workflow with of snakefiles, resolve dependencies with pip and brew, and

## 2. First steps

1. Installation

    ```sh
    git clone https://github.com/jlanga/smsk.git my_project # Clone
    virtualenv --python=python3 bin/py3                     # Create an environment
    bash scripts/install_brew.sh                            # Download linuxbrew
    ```

2. Activate the environment (`deactivate` to deactivate):
    ```sh
    source bin/py3/bin/activate
    ```

3. Install software and packages via pip and homebrew (edit whatever is necessary):

    ```sh
    bash scripts/install_software.sh
    ```
4. Execute the pipeline:

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



## 4. Writting workflows considerations

- The workflow should be written in the main `Snakefile` and all the subworkflows in `scripts/snakefiles`.

- Split into different snakefiles as much as possible. This way code supervision is more bearable and you can recycle them for other projects.

- Configuration for software, samples, etcetera, should be written in the `config.yaml` (instead of hardwritting then somewhere in a 1000 lines script).

- Install as many possible packages from `brew` and `pip` instead of using `apt`/`apt-get`: software is more recent this way, and you don't have to unzip tarballs. This way your workflow is more reproducible. The problem I see is that you cannot specify exact versions in `brew`.

- Use as much as possible `temp()` and `protected()` so you save space and also protect yourself from deleting everything.

- Pipe and compress as much as possible.

- End each subworkflow with a report for your own sanity.

- Use in command line applications long flags (`wget --continue $URL`): this way it is more readable. The computer does not care.



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
