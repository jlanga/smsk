shell.prefix("set -euo pipefail;")


# Load configuration dicts
configfile: "features.yml"
features = config.copy()

configfile: "samples.yml"
samples = config.copy()

configfile: "params.yml"
params = config.copy()

del config


# Define cross-script variables
MAX_THREADS = params["max_threads"]

# Import ubworkflows
snakefiles = "src/snakefiles/"
include: snakefiles + "folders.smk"
include: snakefiles + "clean.smk"
include: snakefiles + "raw.smk"
include: snakefiles + "map.smk"
include: snakefiles + "call.smk"

rule all:
    """
    Execute the entire pipeline
    """
    input:
        REPORT_CALL + "call.html"
