configfile: "config.yaml"

rule all:
    input:
        "doc/report.html"

include: "scripts/snakefiles/clean.snakefile"
include: "scripts/snakefiles/raw.snakefile"
include: "scripts/snakefiles/map.snakefile"
include: "scripts/snakefiles/call.snakefile"
include: "scripts/snakefiles/report.snakefile"
