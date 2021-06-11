#!/bin/bash

# Other necessities
apt-get update
echo "ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true" | debconf-set-selections
apt-get install -y wget unzip curl python-pydot python-pydot-ng graphviz ttf-mscorefonts-installer git