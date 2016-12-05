rule map_bwa_index:
    input:
        raw_dir + "genome.fa"
    output:
        expand(
            raw_dir + "genome.fa.{extension}",
            extension = "amb ann bwt pac sa".split(" ")
        )
    threads:
        1
    log:
        map_dir + "bwa_index.log"
    benchmark:
        map_dir + "bwa_index.json"
    shell:
        "bwa index "
            "{input} "
        "2> {log}"



rule map_bwa_sample:
    input:
        genome = raw_dir + "genome.fa",
        sample = lambda wildcards: config["samples"][wildcards.sample],
        ref_files = expand(
            raw_dir + "genome.fa.{extension}",
            extension = "amb ann bwt pac sa".split(" ")
        )
    output:
        temp(
            map_dir + "{sample}.unsorted.bam"
        )
    params:
        rg="@RG\tID:{sample}\tSM:{sample}"
    log:
        map_dir + "bwa_{sample}.log"
    benchmark:
        map_dir + "bwa_{sample}.json"
    threads:
        8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input.genome} {input.sample} | "
        "samtools view -Sb - "
        "> {output}) "
        "2> {log}"



rule map_sort_sample:
    input:
        map_dir + "{sample}.unsorted.bam"
    output:
        protected("results/map/{sample}.sorted.bam")
    log:
        map_dir + "sort_{sample}.log"
    benchmark:
        map_dir + "sort_{sample}.json"
    shell:
        "samtools sort "
            "-T $(mktemp --dry-run) "
            "-O bam {input} "
        "> {output} "
        "2> {log}"



rule map_index_sample:
    input:
        map_dir + "{sample}.sorted.bam"
    output:
        protected(map_dir + "{sample}.sorted.bam.bai")
    log:
        map_dir + "index_{sample}.log"
    benchmark:
        map_dir + "index_{sample}.json"
    shell:
        "samtools index {input} > {log} 2>&1"


rule map:
    input:
        expand(
            map_dir + "{sample}.sorted.bam.bai",
            sample = config["samples"]
        )
