# smsk: A Snakemake skeleton to jumpstart projects

## 1. Description

This is a small skeleton to create Snakemake workflows. [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) "is a workflow management system that aims to reduce the complexity of creating workflows by providing a fast and comfortable execution environment, together with a clean and modern specification language in python style."

The idea is to create a workflow with of snakefiles, resolve dependencies with pip and brew, and

## 2. First steps

1. Installation

    ```sh
    git clone https://github.com/jlanga/smsk.git my_project # Clone
    cd my_project
    virtualenv --python=python3 bin/py3                     # Create an environment
    git clone https://github.com/Linuxbrew/brew.git         # Download linuxbrew
    ```

2. Activate the environment (`deactivate` to deactivate):
    ```sh
    source bin/activate
    ```

3. Install software and packages via pip and homebrew (edit whatever is necessary):

    ```sh
    bash scripts/install_software.sh
    ```
4. Execute the pipeline:

    ```sh
    snakemake
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

- Start each rule name with the name of the subworkflow (`map`), and mark that it is executed over a item (`_sample`): `map_bowtie_sample`, `map_sort_sample`, `map_index_sample`.

- Use a snakefile to store all the folder names instead of typing them explicitelly (`scripts/snakefiles/folders.snakefile`), and using variables with the convention `SUBWORKFLOW_NAME`: `map_bwa_sample`, `map_sort_sample`, etc.

- End a workflow with a checkpoint rule: a rule that takes as input the result of the workflow (`map`). Use the subworkflow name as a folder name to store its results: `map` results go into `results/map/`.

- Log everything. Store it next to the results: `rule/rest_of_rule_name_sample.log`. Store also benchmarks in JSON format.

- End it also with a clean rule that deletes everything of the workflow (`clean_map`).

- Use the `scripts/snakefiles/raw.snakefile` to get/link your raw data, databases, etcetera. You should be careful when cleaning this folder.

- Configuration for software, samples, etcetera, should be written in the `config.yaml` (instead of hardwritting them somewhere in a 1000 line script).

- `shell.prefix("set -euo pipefail;")` in the first line of the Snakefile makes the entire workflow to stop in case of even a warning.

- If compressing, use `pigz`, `pbzip2` or `pxz` instead of `gzip`. Get them from `brew`.

- Install as many possible packages from `brew` and `pip` instead of using `apt`/`apt-get`: software is more recent this way, and you don't have to unzip tarballs. This way your workflow is more reproducible. The problem I see is that you cannot specify exact versions in `brew`.

- To install software from tarballs, download them into `src/` and copy them to `bin/` (and write the steps in `scripts/install_software.sh`):

    ```sh
    # Binaries are already compiled
    wget \
        --continue \
        --output-document src/bin1.tar.gz \
        http://bin1.com/bin1.tar.gz
    tar xvf src/bin1.tar.gz
    cp src/bin1/bin1 bin/ # or link
    
    # Tarball contains the source
    wget \
        --continue \
        --output-document src/bin2.tar.gz \
        http://bin2.com/bin2.tar.gz
    tar xvf src/bin2.tar.gz
    pushd src/bin2/
    make -j
    cp build/bin2 ../../bin/
    ```

- Use as much as possible `temp()` and `protected()` so you save space and also protect yourself from deleting everything.

- Pipe and compress as much as possible. Make use of the [process substitution](http://vincebuffalo.org/blog/2013/08/08/using-names-pipes-and-process-substitution-in-bioinformatics.html) feature in `bash`: `cmd <(gzip -dc fa.gz)` and `cmd >(gzip -9 > file.gz)`. The problem is that it is hard to estimate the CPU usage of each step of the workflow.

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
