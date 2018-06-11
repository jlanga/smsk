shell.prefix("set -euo pipefail;")

configfile: "features.yml"
features = config.copy()

configfile: "samples.yml"
samples = config.copy()

configfile: "params.yml"
params = config.copy()

del config

MAX_THREADS = params["max_threads"]

snakefiles = "src/snakefiles/"

include: snakefiles + "folders.py"
include: snakefiles + "clean.py"
include: snakefiles + "raw.py"
include: snakefiles + "map.py"
include: snakefiles + "call.py"

rule all:
    input:
        REPORT_CALL + "call.html"
