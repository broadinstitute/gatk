# Using OpenJDK 8
# This Dockerfile does not require any files that are in the GATK4 repo.
FROM ubuntu:16.04

#### Basic image utilities
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y python && \
    apt-get install -y python3-pip && \
    apt-get install --no-install-recommends -y wget curl unzip gcc python-dev python-setuptools git less lynx ## hdfview && \
    apt-get install -y \
    aufs-tools \
    automake \
    bedtools \
    btrfs-tools \
    build-essential \
    dpkg-sig \
    iptables \
    samtools \
    software-properties-common && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove

#### Specific for google cloud support
RUN export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update -y && apt-get install google-cloud-sdk -y && \
    apt-get -y autoremove && \
    apt-get -y clean && \
    yes | gcloud components update && \
    yes | gcloud components update preview
###########

# Get OpenJDK8 on here.
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get clean

# Set environment variables.
ENV HOME /root

# Define working directory.
WORKDIR /root

# Define default command.
CMD ["bash"]

ENV JAVA_LIBRARY_PATH /usr/lib/jni
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/

RUN java -version

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 && \
    add-apt-repository "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" && \
    apt-get update && \
    apt-get install -y --force-yes \
        r-base-dev=3.2.5-1xenial \
        r-base-core=3.2.5-1xenial && \
    apt-get -y clean
    apt-get -y autoremove

COPY install_R_packages.R install_R_packages.R
RUN Rscript install_R_packages.R

RUN apt-get clean
RUN apt-get autoclean
RUN apt-get autoremove
