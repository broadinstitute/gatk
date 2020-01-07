#!/bin/bash

# Install the github repo version
git clone https://github.com/fastai/fastai
cd fastai

# Peg our version to a known-working SHA, since they make 
# post-1.0 breaking changes literally every day...
git reset --hard 14868ca69483afbaa8e28d4e281c148d1dad1c89

tools/run-after-git-clone
pip install -e .[dev]
