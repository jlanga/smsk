rule clean_map:
    shell:
        "rm -rf " + MAP



rule clean_call:
    shell:
        "rm -rf " + CALL



rule clean_report:
    shell:
        "rm -rf " + REPORT_CALL



rule clean:
    shell:
        "rm -rf " +
            MAP + " " +
            CALL + " " +
            REPORT_CALL
