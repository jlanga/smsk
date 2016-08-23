rule map_bwa_sample:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp("results/map/{sample}.unsorted.bam")
    params:
        rg="@RG\tID:{sample}\tSM:{sample}"
    log:
        "results/map/bwa_{sample}.log"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"



rule map_sort_sample:
    input:
        "results/map/{sample}.unsorted.bam"
    output:
        protected("results/map/{sample}.bam")
    shell:
        "samtools sort "
            "-T $(mktemp --dry-run) "
            "-O bam {input} > {output}"



rule map_index_sample:
    input:
        "results/map/{sample}.bam"
    output:
        "results/map/{sample}.bam.bai"
    shell:
        "samtools index {input}"


rule map:
    input:
        expand(
            "results/map/{sample}.bam.bai",
            sample = config["samples"]
        )
