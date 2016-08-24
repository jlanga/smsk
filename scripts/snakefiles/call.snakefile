rule call_bcftools:
    input:
        fa="data/genome.fa",
        bam=expand(
            map_dir + "{sample}.bam",
            sample=config["samples"]
        ),
        bai=expand(
            map_dir + "{sample}.bam.bai",
            sample=config["samples"]
        )
    output:
        protected(call_dir + "all.vcf")
    log:
        call_dir + "bcftools.log"
    benchmark:
        call_dir + "bcftools.json"
    shell:
        "( samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output} ) 2> {log}"



rule call:
    input:
        call_dir + "all.vcf"



rule call_report:
    input:
        "results/call/all.vcf"
    output:
        protected(call_doc + "call.html")
    run:
        from snakemake.utils import report
        with open(input[0]) as vcf:
            n_calls = sum(1 for l in vcf if not l.startswith("#"))

        report("""
        An example variant calling workflow
        ===================================

        Reads were mapped to the Yeast
        reference genome and variants were called jointly with
        SAMtools/BCFtools.

        This resulted in {n_calls} variants (see Table T1_).
        """, output[0], T1=input[0])
