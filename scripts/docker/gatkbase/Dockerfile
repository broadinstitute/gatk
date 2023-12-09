# Using OpenJDK 17
# This Dockerfile does not require any files that are in the GATK4 repo.
FROM ubuntu:20.04

# Avoid interactive prompts during apt installs/upgrades
ENV DEBIAN_FRONTEND noninteractive

#### Basic image utilities
RUN apt update && \
    apt full-upgrade -y && \
    apt install -y --no-install-recommends \
    python3 \
    python3-pip \
    wget \
    curl \
    bc \
    unzip \
    bzip2 \
    less \
    bedtools \
    samtools \
    bcftools \
    tabix \
    git \
    gpg-agent \
    build-essential \
    openjdk-17-jdk \
    vim \
    software-properties-common && \
    apt -y clean  && \
    apt -y autoclean  && \
    apt -y autoremove && \
    rm -rf /var/lib/apt/lists/*

RUN java -version

#### Specific for google cloud support
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" \
    | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
    | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt update -y && \
    apt install -y --no-install-recommends google-cloud-cli && \
    apt -y clean  && \
    apt -y autoclean  && \
    apt -y autoremove && \
    rm -rf /var/lib/apt/lists/* && \
#  Remove the anthos cli tool and related files since they are very large and we / anyone using the docker are unlikely to use them
#  Remove the bundled python because we have python installed separately
    rm -rf /usr/lib/google-cloud-sdk/bin/anthoscli /usr/lib/google-cloud-sdk/platform/anthoscli_licenses /usr/lib/google-cloud-sdk/platform/bundledpythonunix && \
    find / -wholename "*__pycache__/*.pyc" -exec rm {} + 

# Set environment variables.
ENV HOME /root

# Define working directory.
WORKDIR /root

# Define default command.
CMD ["bash"]

ENV JAVA_LIBRARY_PATH /usr/lib/jni

# Install miniconda
ENV DOWNLOAD_DIR /downloads
ENV CONDA_URL https://repo.anaconda.com/miniconda/Miniconda3-py310_23.1.0-1-Linux-x86_64.sh
ENV CONDA_SHA256 "32d73e1bc33fda089d7cd9ef4c1be542616bd8e437d1f77afeeaf7afdb019787"
ENV CONDA_PATH /opt/miniconda
ENV PATH $CONDA_PATH/bin:$PATH
RUN mkdir $DOWNLOAD_DIR && \
    wget -nv -O $DOWNLOAD_DIR/miniconda.sh $CONDA_URL && \
    test "$(sha256sum $DOWNLOAD_DIR/miniconda.sh | awk -v FS=' ' -v ORS='' '{print $1}')" = "$CONDA_SHA256" && \
    bash $DOWNLOAD_DIR/miniconda.sh -p $CONDA_PATH -b && \
    rm $DOWNLOAD_DIR/miniconda.sh && \
    conda clean -afy && \
    rm -rf /root/.cache/pip
