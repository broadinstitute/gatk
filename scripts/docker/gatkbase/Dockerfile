# Using OpenJDK 11
# This Dockerfile does not require any files that are in the GATK4 repo.
FROM ubuntu:20.04

#### Basic image utilities
#Set noninteractive and timezone or tzdata install stops to ask where you live
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Boston
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    python \
    wget \
    curl \
    bc \
    unzip \
    bzip2 \
    less \
    bedtools \
    samtools \
    openjdk-11-jdk \
    tabix \
    gpg-agent \
    build-essential \
    software-properties-common && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove

#### install and configure gcloud
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get update -y && \
    apt-get install google-cloud-sdk -y && \ 
    apt-get -y clean && \
    apt-get -y autoclean && \
    apt-get -y autoremove 

RUN gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true 
      

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
ENV CONDA_URL https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh
ENV CONDA_MD5 = "81c773ff87af5cfac79ab862942ab6b3"
ENV CONDA_PATH /opt/miniconda
RUN mkdir $DOWNLOAD_DIR && \
    wget -nv -O $DOWNLOAD_DIR/miniconda.sh $CONDA_URL && \
    test "`md5sum $DOWNLOAD_DIR/miniconda.sh | awk -v FS='  ' '{print $1}'` = $CONDA_MD5" && \
    bash $DOWNLOAD_DIR/miniconda.sh -p $CONDA_PATH -b && \
    rm $DOWNLOAD_DIR/miniconda.sh

# Deleting unneeded caches
RUN rm -rf /var/lib/apt/lists/*
