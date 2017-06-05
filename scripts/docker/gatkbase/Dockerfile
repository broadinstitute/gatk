# Using OpenJDK 8
# This Dockerfile does not require any files that are in the GATK4 repo.
FROM ubuntu:16.04

#### Basic image utilities
RUN apt-get update && \
     apt-get upgrade -y && \
     apt-get install -y python && \
     apt-get install -y python3-pip && \
     apt-get install -y wget curl unzip gcc python-dev python-setuptools git less lynx hdfview
##########

# Install some more useful tools
RUN apt-get install -y \
 aufs-tools \
 automake \
 bedtools \
 btrfs-tools \
 build-essential \
 dpkg-sig \
 iptables \
 samtools \
 software-properties-common

#### Specific for google cloud support
RUN wget https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.zip \
    && unzip google-cloud-sdk.zip \
    && rm google-cloud-sdk.zip

RUN google-cloud-sdk/install.sh --usage-reporting=true --path-update=true --bash-completion=true --rc-path=/.bashrc --disable-installation-options
VOLUME ["/root/.config"]
ENV PATH /google-cloud-sdk/bin:$PATH

RUN yes | gcloud components update
RUN yes | gcloud components update preview
###########

# Get OpenJDK8 on here.
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk

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
        r-base-core=3.2.5-1xenial

COPY install_R_packages.R install_R_packages.R
RUN Rscript install_R_packages.R
