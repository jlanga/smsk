shell.prefix("set -euo pipefail;")
configfile: "src/config.yaml"


snakefiles = "src/snakefiles/"

include: snakefiles + "folders.py"
include: snakefiles + "clean.py"
include: snakefiles + "raw.py"
include: snakefiles + "map.py"
include: snakefiles + "call.py"

rule all:
    input:
        call_doc + "call.html"
