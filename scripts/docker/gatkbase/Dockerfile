# Using OpenJDK 8
# This Dockerfile does not require any files that are in the GATK4 repo.
FROM ubuntu:16.04

ENV DOWNLOAD_DIR /downloads
ENV CONDA_URL https://repo.continuum.io/miniconda/Miniconda3-4.3.30-Linux-x86_64.sh
ENV CONDA_MD5 = "0b80a152332a4ce5250f3c09589c7a81"
ENV CONDA_PATH /opt/miniconda

#### Install basic utilities / useful tools.
RUN apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y \
        python \
        python3-pip && \
    apt-get install --no-install-recommends -y \
        wget \
        curl \
        unzip \
        gcc \
        python-dev \
        python-setuptools \
        git \
        less \
        lynx \
        software-properties-common \
        openjdk-8-jdk && \
    echo "deb http://packages.cloud.google.com/apt cloud-sdk-xenial main" | \
        tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 && \
    add-apt-repository "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" && \
    apt-get update -y && \
    apt-get install -y \
        google-cloud-sdk \
        r-base-dev=3.2.5-1xenial \
        r-base-core=3.2.5-1xenial && \
    apt-get clean -y && \
    apt-get autoremove -y && \
    rm -r /var/lib/apt/lists/* && \
    mkdir $DOWNLOAD_DIR && \
    wget -nv -O $DOWNLOAD_DIR/miniconda.sh $CONDA_URL && \
    test "`md5sum $DOWNLOAD_DIR/miniconda.sh | awk -v FS='  ' '{print $1}'` = $CONDA_MD5" && \
    bash $DOWNLOAD_DIR/miniconda.sh -p $CONDA_PATH -b && \
    rm $DOWNLOAD_DIR/miniconda.sh

ENV PATH $CONDA_PATH/envs/gatk/bin:$CONDA_PATH/bin:$PATH

#### Specific for google cloud support
VOLUME ["/root/.config"]
###########

# Set environment variables.
ENV HOME /gatk
ENV JAVA_LIBRARY_PATH /usr/lib/jni
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/

# Define working directory.
WORKDIR /gatk

# Define default command.
CMD ["bash"]

RUN java -version

COPY install_R_packages.R install_R_packages.R
RUN Rscript install_R_packages.R
