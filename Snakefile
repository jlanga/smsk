shell.prefix("set -euo pipefail;")
configfile: "config.yaml"


snakefiles = "bin/snakefiles/"

include: snakefiles + "folders"
include: snakefiles + "clean"
include: snakefiles + "raw"
include: snakefiles + "map"
include: snakefiles + "call"

rule all:
    input:
        call_doc + "call.html"
