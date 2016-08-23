rule clean_map:
    shell:
        "rm -rf results/map/"



rule clean_call:
    shell:
        "rm -rf results/call/"



rule clean_report:
    shell:
        "rm doc/report.html"



rule clean:
    shell:
        "rm -rf "
            "results/map/ "
            "results/call/ "
            "doc/report.html"
