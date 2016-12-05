rule raw_download_tarball:
    output:
        tarball = raw_dir + "snakemake-tutorial-data.tar.gz"
    threads:
        1
    params:
        url = config["urls"]["tarball"]
    log:
        raw_dir + "download_tarball.log"
    benchmark:
        raw_dir + "download_tarball.json"
    shell:
        "wget "
            "--continue "
            "--output-document {output.tarball} "
            "{params.url} "
        "2> {log} 1>&2"



rule raw_extract_genome:
    input:
        tarball = raw_dir + "snakemake-tutorial-data.tar.gz"
    output:
        touch(raw_dir + "genome.fa")
    threads:
        1
    log:
        raw_dir + "extract_genome.log"
    benchmark:
        raw_dir + "extract_genome.json"
    shell:
        "tar "
            "--extract "
            "--verbose "
            "--file {input.tarball} "
            "--directory {raw_dir} "
            "--strip 1 "
            "data/genome.fa "
        "2> {log} 1>&2"



rule raw_extract_samples:
    input:
        tarball =  raw_dir + "snakemake-tutorial-data.tar.gz"
    output:
        a = touch(raw_dir + "samples/A.fastq"),
        b = touch(raw_dir + "samples/B.fastq")
    threads:
        1
    params:
        a_path = "data/samples/A.fastq",
        b_path = "data/samples/B.fastq"
    log:
        raw_dir + "extract_genome.log"
    benchmark:
        raw_dir + "extract_genome.json"
    shell:
        "tar "
            "--extract "
            "--verbose "
            "--file {input.tarball} "
            "--directory {raw_dir} "
            "--strip 1 "
            "{params.a_path} {params.b_path} "
        "2> {log} 1>&2"
