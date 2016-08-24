shell.prefix("set -euo pipefail;")
configfile: "config.yaml"


snakefiles = "scripts/snakefiles/"

include: snakefiles + "folders.snakefile"
include: snakefiles + "clean.snakefile"
include: snakefiles + "raw.snakefile"
include: snakefiles + "map.snakefile"
include: snakefiles + "call.snakefile"

rule all:
    input:
        call_doc + "call.html"
