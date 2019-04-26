#!/bin/bash

# Fetch the latest available version
wget https://storage.googleapis.com/ml4cvd/ml4cvd-master.zip
unzip ml4cvd-master.zip
cd ml4cvd-master/pyukbb

tools/run-after-git-clone
pip install -e .[dev]
