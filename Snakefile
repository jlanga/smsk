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
include: snakefiles + "folders.py"
include: snakefiles + "clean.py"
include: snakefiles + "raw.py"
include: snakefiles + "map.py"
include: snakefiles + "call.py"

rule all:
    """
    Execute the entire pipeline
    """
    input:
        REPORT_CALL + "call.html"
