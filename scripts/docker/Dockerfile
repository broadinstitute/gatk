FROM ubuntu:14.04

## Not supported by docker-compose

ARG GITHUB_TAG
ARG GITHUB_DIR=tags/

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
ENV JAVA_LIBRARY_PATH /usr/lib/jni

# Added from GATK4 public (Getting R setup)
RUN mkdir -p /usr/local/lib/R/ && \
    mkdir -p ~/site-library && \
    ln -sFv ~/site-library /usr/local/lib/R/site-library
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 && \
    add-apt-repository "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" && \
    apt-get update && \
    apt-get install -y --force-yes \
        r-base-dev=3.1.3-1trusty \
        r-base-core=3.1.3-1trusty

# R setup complete...

# Installing GATK4 protected (from repo: gatk-protected)
#  This Dockerfile is getting the specified tag of gatk4-protected
#  This Dockerfile generates a jar file that will work without spark or in spark standalone.  Not on a spark cluster.

RUN git clone https://github.com/broadinstitute/gatk-protected.git

# Install R dependencies and build the shadowJar
ENV GITHUB_TAG ${GITHUB_TAG}
ENV GITHUB_DIR ${GITHUB_DIR}
WORKDIR /root/gatk-protected/
RUN git checkout ${GITHUB_DIR}${GITHUB_TAG}
RUN Rscript /root/gatk-protected/scripts/install_R_packages.R
RUN ./gradlew clean compileTestJava shadowJar

WORKDIR /root

# Make sure we can see a help message
RUN ln -sFv /root/gatk-protected/build/libs/gatk-protected.jar
RUN java -jar gatk-protected.jar -h

# Install git lfs and get latest big files
WORKDIR /root/gatk-protected/
RUN bash scripts/install_git_lfs.sh

RUN echo This docker image is running gatk-protected `git describe --tags` > /root/GATK_PROTECTED_VERSION

# Create a simple unit test runner
ENV CI true
RUN echo "cd /root/gatk-protected/ && ./gradlew test" >/root/run_unit_tests.sh

WORKDIR /root