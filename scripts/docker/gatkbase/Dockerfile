# Using OpenJDK 8
# This Dockerfile does not require any files that are in the GATK4 repo.
FROM ubuntu:22.04

#### Basic image utilities
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    python3 \
    wget \
    curl \
    bc \
    unzip \
    bzip2 \
    less \
    bedtools \
    samtools \
    openjdk-8-jdk \
    tabix \
    gpg-agent \
    build-essential \
    software-properties-common && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove

#### Specific for google cloud support
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends  \
    apt-transport-https && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove 

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    gnupg && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | \
    tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | \
    tee /usr/share/keyrings/cloud.google.gpg && \
    apt-get update -y && \
    apt-get install google-cloud-sdk -y && \
    apt-get -y autoremove && \
    apt-get -y clean
###########

# Set environment variables.
ENV HOME /root

# Define working directory.
WORKDIR /root

# Define default command.
CMD ["bash"]

ENV JAVA_LIBRARY_PATH /usr/lib/jni
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/

RUN java -version

# Install miniconda
ENV DOWNLOAD_DIR /downloads
ENV CONDA_URL=https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh
ENV CONDA_MD5=3190da6626f86eee8abf1b2fd7a5af492994eb2667357ee4243975cdbb175d7a
ENV CONDA_PATH /opt/miniconda
RUN mkdir $DOWNLOAD_DIR && \
    wget -nv -O $DOWNLOAD_DIR/miniconda.sh $CONDA_URL && \
    test "`md5sum $DOWNLOAD_DIR/miniconda.sh | awk -v FS='  ' '{print $1}'` = $CONDA_MD5" && \
    bash $DOWNLOAD_DIR/miniconda.sh -p $CONDA_PATH -b && \
    rm $DOWNLOAD_DIR/miniconda.sh

# Deleting unneeded caches
RUN rm -rf /var/lib/apt/lists/*
