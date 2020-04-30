rule map_bwa_index:
    """Index the reference genome with BWA."""
    input:
        RAW + "genome.fa"
    output:
        expand(
            RAW + "genome.fa.{extension}",
            extension="amb ann bwt pac sa".split(" ")
        )
    log: MAP + "bwa_index.log"
    benchmark: MAP + "bwa_index.bmk"
    conda: "map.yml"
    shell: "bwa index {input} 2> {log}"


def get_sample(wildcards):
    """Get the sample name inferred from the wildcards"""
    return samples["samples"][wildcards.sample]


rule map_bwa:
    """Map a sample to reference with BWA. Convert online from SAM to BAM."""
    input:
        genome = RAW + "genome.fa",
        sample = get_sample,
        ref_files = expand(
            RAW + "genome.fa.{extension}",
            extension="amb ann bwt pac sa".split(" ")
        )
    output: temp(MAP + "{sample}.unsorted.bam")
    log: MAP + "bwa_{sample}.log"
    benchmark: MAP + "bwa_{sample}.bmk"
    threads: MAX_THREADS
    conda: "map.yml"
    params:
        rg = "@RG\tID:{sample}\tSM:{sample}"
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input.genome} {input.sample} "
        "| samtools view -Sb - "
        "> {output}) "
        "2> {log}"


rule map_sort:
    """Sort BAM"""
    input: MAP + "{sample}.unsorted.bam"
    output: protected("results/map/{sample}.sorted.bam")
    log: MAP + "sort_{sample}.log"
    benchmark: MAP + "sort_{sample}.bmk"
    conda: "map.yml"
    shell:
        "samtools sort -T $(mktemp --dry-run) -O bam {input} > {output} 2> {log}"


rule map_index:
    """Index BAM to generate a .bai file"""
    input: MAP + "{sample}.sorted.bam"
    output: protected(MAP + "{sample}.sorted.bam.bai")
    log: MAP + "index_{sample}.log"
    benchmark: MAP + "index_{sample}.bmk"
    conda: "map.yml"
    shell: "samtools index {input} > {log} 2>&1"


rule map:
    """Map all samples to reference, sort and index."""
    input:
        expand(
            MAP + "{sample}.sorted.bam.bai",
            sample=samples["samples"]
        )
