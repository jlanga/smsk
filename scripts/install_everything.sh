#!/usr/bin/env bash

virtualenv --python=python3 bin/py3

bash scripts/install_brew.sh

source bin/activate

bash scripts/install_software.sh
