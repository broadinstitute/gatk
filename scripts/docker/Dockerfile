# VERSION 0.1 -- very basic setup of an image.
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

RUN wget -N https://downloads.gradle.org/distributions/gradle-2.7-bin.zip
RUN unzip gradle-2.7-bin.zip
ENV PATH=/root/gradle-2.7/bin:$PATH
RUN echo "${PATH}"
RUN gradle --version

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

# Added from hellbender public (Getting R setup)
RUN mkdir -p /usr/local/lib/R/
RUN mkdir -p ~/site-library
RUN ln -sFv ~/site-library /usr/local/lib/R/site-library
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN add-apt-repository "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/"
RUN apt-get update
RUN apt-get install -y --force-yes r-base-dev=3.1.3-1trusty
RUN apt-get install -y --force-yes r-base-core=3.1.3-1trusty

# R setup complete...

# Installing GATK4 protected (hellbender-protected)
#  This Dockerfile is getting the latest from master branch of gatk4 protected (hellbender-protected)"
#  This Dockerfile generates a jar file that will work without spark or in spark standalone.  Not on a spark cluster."
#  Unit tests are NOT being run"

RUN git clone https://github.com/broadinstitute/hellbender-protected.git

# Install custom R packages
RUN Rscript /root/hellbender-protected/scripts/install_R_packages.R

# Build the shadowJar

WORKDIR /root/hellbender-protected/
RUN gradle clean shadowJar

WORKDIR /root

# Make sure we can see a help message
RUN ln -sFv /root/hellbender-protected/build/libs/hellbender-protected.jar
RUN java -jar hellbender-protected.jar -h

ENV JAVA_LIBRARY_PATH /usr/lib/jni
