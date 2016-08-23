rule call_bcftools:
    input:
        fa="data/genome.fa",
        bam=expand(
            "results/map/{sample}.bam",
            sample=config["samples"]
        ),
        bai=expand(
            "results/map/{sample}.bam.bai",
            sample=config["samples"]
        )
    output:
        "results/call/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"



rule call:
    input:
        "results/call/all.vcf"
