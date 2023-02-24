# Docker image with a grab bag of utilities for Cromwell on Azure exploration spikes. Not currently optimized for size
# or anything else, this is currently just all the potentially useful things.
FROM ubuntu:20.04

# Azure CLI
# https://learn.microsoft.com/en-us/cli/azure/install-azure-cli-linux?pivots=apt#option-2-step-by-step-installation-instructions
RUN apt-get update
RUN apt-get install --assume-yes ca-certificates curl apt-transport-https lsb-release gnupg

RUN mkdir -p /etc/apt/keyrings
RUN curl -sLS https://packages.microsoft.com/keys/microsoft.asc | \
    gpg --dearmor | \
    tee /etc/apt/keyrings/microsoft.gpg > /dev/null
RUN chmod go+r /etc/apt/keyrings/microsoft.gpg

# ENV AZ_REPO=$(lsb_release -cs)
# Hardcode to focal/20.04 for consistency with the base image above and sqlcmd setup below
ENV AZ_REPO=focal
RUN echo "deb [arch=`dpkg --print-architecture` signed-by=/etc/apt/keyrings/microsoft.gpg] https://packages.microsoft.com/repos/azure-cli/ $AZ_REPO main" | \
    tee /etc/apt/sources.list.d/azure-cli.list

RUN apt-get update
RUN apt-get install --assume-yes azure-cli

# Install sqlcmd (Microsoft SQL client)
# https://learn.microsoft.com/en-us/sql/linux/sql-server-linux-setup-tools?view=sql-server-ver16&tabs=ubuntu-install%2Credhat-offline#install-tools-on-linux
# Also sneak in an installation of the driver for Microsoft databases via `msodbcsql18`.
RUN curl https://packages.microsoft.com/keys/microsoft.asc | \
    apt-key add -

RUN curl https://packages.microsoft.com/config/ubuntu/20.04/prod.list | \
    tee /etc/apt/sources.list.d/msprod.list

RUN apt-get update

# sneaky EULA "acceptance" https://stackoverflow.com/a/42383714
ENV ACCEPT_EULA=Y

# ODBC and Microsoft ODBC SQL driver
RUN apt-get install --assume-yes mssql-tools unixodbc-dev msodbcsql18
ENV PATH=$PATH:/opt/mssql-tools/bin

# Python
RUN apt-get install --assume-yes python3-pip
COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

# Temurin 11 JDK
# https://askubuntu.com/a/1386901
RUN apt-get install --assume-yes wget
RUN wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | apt-key add -
RUN echo "deb https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" | tee /etc/apt/sources.list.d/adoptium.list
RUN apt update && apt install --assume-yes temurin-11-jdk

# Coursier / Ammonite for scripting in the Java ecosystem
# https://get-coursier.io/docs/cli-installation#linux
#
# Use the statically linked version for now to get around a broken dynamically linked launcher
# https://github.com/coursier/coursier/issues/2624
# https://stackoverflow.com/a/75232986/21269164
RUN curl -fL "https://github.com/coursier/launchers/raw/master/cs-x86_64-pc-linux-static.gz" | gzip -d > /usr/local/bin/cs
RUN chmod +x /usr/local/bin/cs
RUN mkdir -p /coursier/bin
ENV PATH=/coursier/bin\:$PATH
RUN cs setup --install-dir /coursier/bin --yes
