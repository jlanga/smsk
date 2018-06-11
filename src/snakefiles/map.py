rule map_bwa_index:
    input:
        RAW + "genome.fa"
    output:
        expand(
            RAW + "genome.fa.{extension}",
            extension="amb ann bwt pac sa".split(" ")
        )
    threads:
        1
    log:
        MAP + "bwa_index.log"
    benchmark:
        MAP + "bwa_index.time"
    conda:
        "map.yml"
    shell:
        "bwa index "
            "{input} "
        "2> {log}"



rule map_bwa_sample:
    input:
        genome = RAW + "genome.fa",
        sample = lambda wildcards: samples["samples"][wildcards.sample],
        ref_files = expand(
            RAW + "genome.fa.{extension}",
            extension="amb ann bwt pac sa".split(" ")
        )
    output:
        temp(MAP + "{sample}.unsorted.bam")
    params:
        rg="@RG\tID:{sample}\tSM:{sample}"
    log:
        MAP + "bwa_{sample}.log"
    benchmark:
        MAP + "bwa_{sample}.time"
    threads:
        MAX_THREADS
    conda:
        "map.yml"
    shell:
        "(bwa mem "
            "-R '{params.rg}' "
            "-t {threads} "
            "{input.genome} "
            "{input.sample} "
        "| samtools view "
            "-Sb "
            "- "
        "> {output}) "
        "2> {log}"



rule map_sort_sample:
    input:
        MAP + "{sample}.unsorted.bam"
    output:
        protected("results/map/{sample}.sorted.bam")
    log:
        MAP + "sort_{sample}.log"
    benchmark:
        MAP + "sort_{sample}.time"
    conda:
        "map.yml"
    shell:
        "samtools sort "
            "-T $(mktemp --dry-run) "
            "-O bam {input} "
        "> {output} "
        "2> {log}"



rule map_index_sample:
    input:
        MAP + "{sample}.sorted.bam"
    output:
        protected(MAP + "{sample}.sorted.bam.bai")
    log:
        MAP + "index_{sample}.log"
    benchmark:
        MAP + "index_{sample}.bmk"
    conda:
        "map.yml"
    shell:
        "samtools index {input} > {log} 2>&1"


rule map:
    input:
        expand(
            MAP + "{sample}.sorted.bam.bai",
            sample=samples["samples"]
        )
