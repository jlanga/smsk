shell.prefix("set -euo pipefail;")
configfile: "src/config.yaml"


snakefiles = "src/snakefiles/"

include: snakefiles + "folders"
include: snakefiles + "clean"
include: snakefiles + "raw"
include: snakefiles + "map"
include: snakefiles + "call"

rule all:
    input:
        call_doc + "call.html"
