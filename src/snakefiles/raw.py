rule raw_download_tarball:
    output:
        tarball = RAW + "snakemake-tutorial-data.tar.gz"
    threads:
        1
    params:
        url = config["urls"]["tarball"]
    log:
        RAW + "download_tarball.log"
    benchmark:
        RAW + "download_tarball.time"
    shell:
        "wget "
            "--continue "
            "--output-document {output.tarball} "
            "{params.url} "
        "2> {log} 1>&2"



rule raw_extract_genome:
    input:
        tarball = RAW + "snakemake-tutorial-data.tar.gz"
    output:
        touch(RAW + "genome.fa")
    log:
        RAW + "extract_genome.log"
    benchmark:
        RAW + "extract_genome.time"
    shell:
        "tar "
            "--extract "
            "--verbose "
            "--file {input.tarball} "
            "--directory {RAW} "
            "--strip 1 "
            "data/genome.fa "
        "2> {log} 1>&2"



rule raw_extract_samples:
    input:
        tarball =  RAW + "snakemake-tutorial-data.tar.gz"
    output:
        a = touch(RAW + "samples/A.fastq"),
        b = touch(RAW + "samples/B.fastq")
    params:
        a_path = "data/samples/A.fastq",
        b_path = "data/samples/B.fastq"
    log:
        RAW + "extract_genome.log"
    benchmark:
        RAW + "extract_genome.time"
    shell:
        "tar "
            "--extract "
            "--verbose "
            "--file {input.tarball} "
            "--directory {RAW} "
            "--strip 1 "
            "{params.a_path} {params.b_path} "
        "2> {log} 1>&2"
