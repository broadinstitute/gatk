# VERSION 0.2 -- very basic setup of an image.
FROM ubuntu:14.04

# Install some tools, including hdfview for HDF5 support and gradle for building the jar
RUN apt-get update && apt-get install -y \
 aufs-tools \
 automake \
 bedtools \
 btrfs-tools \
 build-essential \
 curl \
 dpkg-sig \
 git \
 iptables \
 hdfview \
 samtools \
 wget \
 curl \
 emacs \
 software-properties-common

# Set environment variables.
ENV HOME /root

# Define working directory.
WORKDIR /root

# Define default command.
CMD ["bash"]

# Installing Java 8.... dockerfile snippet from https://github.com/dockerfile/java/blob/master/oracle-java8/Dockerfile
RUN \
  echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | debconf-set-selections && \
  add-apt-repository -y ppa:webupd8team/java && \
  apt-get update && \
  apt-get install -y oracle-java8-installer && \
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /var/cache/oracle-jdk8-installer

# Define commonly used JAVA_HOME variable
ENV JAVA_HOME /usr/lib/jvm/java-8-oracle

# Added from GATK4 public (Getting R setup)
RUN mkdir -p /usr/local/lib/R/
RUN mkdir -p ~/site-library
RUN ln -sFv ~/site-library /usr/local/lib/R/site-library
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN add-apt-repository "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/"
RUN apt-get update
RUN apt-get install -y --force-yes r-base-dev=3.1.3-1trusty
RUN apt-get install -y --force-yes r-base-core=3.1.3-1trusty

# R setup complete...

# Installing GATK4 protected (from repo: gatk-protected)
#  This Dockerfile is getting the specified tag of gatk4-protected
#  This Dockerfile generates a jar file that will work without spark or in spark standalone.  Not on a spark cluster.
#  Unit tests are NOT being run

RUN git clone https://github.com/broadinstitute/gatk-protected.git

# Install custom R packages
RUN Rscript /root/gatk-protected/scripts/install_R_packages.R

# Build the shadowJar
ENV GATK_PROTECTED_VERSION 1.0.0.0-alpha1.3
WORKDIR /root/gatk-protected/
RUN git checkout tags/${GATK_PROTECTED_VERSION}
RUN ./gradlew clean shadowJar

WORKDIR /root

# Make sure we can see a help message
RUN ln -sFv /root/gatk-protected/build/libs/gatk-protected.jar
RUN java -jar gatk-protected.jar -h

ENV JAVA_LIBRARY_PATH /usr/lib/jni

WORKDIR /root/gatk-protected/
RUN echo This docker image is running gatk-protected `git describe --tags` > /root/GATK_PROTECTED_VERSION

WORKDIR /root
