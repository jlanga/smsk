shell.prefix("set -euo pipefail;")
configfile: "src/config.yaml"

MAX_THREADS = config["max_threads"]

snakefiles = "src/snakefiles/"

include: snakefiles + "folders.py"
include: snakefiles + "clean.py"
include: snakefiles + "raw.py"
include: snakefiles + "map.py"
include: snakefiles + "call.py"

rule all:
    input:
        REPORT_CALL + "call.html"
