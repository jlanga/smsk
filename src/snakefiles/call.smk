rule call_bcftools:
    """
    Generate VCF with genotypes and call SNPs with bcftools
    """
    input:
        fa = RAW + "genome.fa",
        bam = expand(
            MAP + "{sample}.sorted.bam",
            sample=samples["samples"]
        ),
        bai = expand(
            MAP + "{sample}.sorted.bam.bai",
            sample=samples["samples"]
        )
    output: protected(CALL + "all.vcf")
    log: CALL + "bcftools.log"
    benchmark: CALL + "bcftools.bmk"
    conda: "call.yml"
    shell:
        "( samtools mpileup -g -f {input.fa} {input.bam} "
        "| bcftools call -mv - > {output} ) 2> {log}"


rule call:
    """
    Checkpoint rule. Generate all.vcf
    """
    input:
        CALL + "all.vcf"


rule call_report:
    """
    Make a simple report in MarkDown for all.vcf
    """
    input:
        CALL + "all.vcf"
    output:
        protected(REPORT_CALL + "call.html")
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
