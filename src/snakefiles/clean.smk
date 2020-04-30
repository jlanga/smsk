rule clean_map:
    """Just remove mappings"""
    shell:
        "rm -rf {MAP}"


rule clean_call:
    """Just remove the SNP calling"""
    shell:
        "rm -rf {CALL}"


rule clean_report:
    """Just remove the report"""
    shell:
        "rm -rf {REPORT_CALL}"


rule clean:
    """Delete everything"""
    shell:
        "rm -rf {MAP} {CALL} {REPORT_CALL}"
