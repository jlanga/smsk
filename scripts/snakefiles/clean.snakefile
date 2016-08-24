rule clean_map:
    shell:
        "rm -rf " + map_dir



rule clean_call:
    shell:
        "rm -rf " + call_dir



rule clean_report:
    shell:
        "rm -rf " + call_doc



rule clean:
    shell:
        "rm -rf " +
            map_dir + " " +
            call_dir + " " +
            call_doc
