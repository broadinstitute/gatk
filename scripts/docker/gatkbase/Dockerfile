# Using OpenJDK 8
# This Dockerfile does not require any files that are in the GATK4 repo.
FROM ubuntu:16.04

#### Basic image utilities
RUN apt-get update && \
     apt-get upgrade -y && \
     apt-get install -y python && \
     apt-get install -y wget curl unzip gcc python-dev python-setuptools git less \
     aufs-tools \
     openjdk-8-jdk \
     build-essential \
     dpkg-sig \
     software-properties-common && \
     apt-get -y autoremove && \
     apt-get -y clean 

#### Specific for google cloud support
RUN export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update -y && apt-get install google-cloud-sdk -y && \
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

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 && \
    add-apt-repository "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" && \
    apt-get update && \
    apt-get install -y --force-yes \
        r-base-dev=3.2.5-1xenial \
        r-base-core=3.2.5-1xenial && \
    apt-get -y autoremove && \
    apt-get -y clean

COPY install_R_packages.R install_R_packages.R
RUN Rscript install_R_packages.R
