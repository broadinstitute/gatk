# Using OpenJDK 17
# This Dockerfile does not require any files that are in the GATK4 repo.
FROM ubuntu:18.04

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
    zlib \
    libbz2 \
    liblzma \
    libcurl \
    less \
    bedtools \
    gpg-agent \
    build-essential \
    openjdk-17-jdk \
    software-properties-common && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt/lists/*

RUN java -version

#### Specific for google cloud support
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" \
    | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
    | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get update -y && \
    apt-get install -y --no-install-recommends google-cloud-cli && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove && \
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
ENV CONDA_MD5 = "32d73e1bc33fda089d7cd9ef4c1be542616bd8e437d1f77afeeaf7afdb019787"
ENV CONDA_PATH /opt/miniconda
ENV PATH $CONDA_PATH/bin:$PATH
RUN mkdir $DOWNLOAD_DIR && \
    wget -nv -O $DOWNLOAD_DIR/miniconda.sh $CONDA_URL && \
    test "`md5sum $DOWNLOAD_DIR/miniconda.sh | awk -v FS='  ' '{print $1}'` = $CONDA_MD5" && \
    bash $DOWNLOAD_DIR/miniconda.sh -p $CONDA_PATH -b && \
    rm $DOWNLOAD_DIR/miniconda.sh && \
    conda clean -afy && \
    find /opt/miniconda/ -follow -type f -name '*.a' -delete && \
    find /opt/miniconda/ -follow -type f -name '*.pyc' -delete && \
    rm -rf /root/.cache/pip 

# Install samtools
ENV HTSLIB_VERSION "1.18"
ENV HTSLIB_VERSION_MD5 = "eef1378084dc3888f6dd59f522886d6c"
RUN wget -q https://github.com/samtools/samtools/releases/download/$HTSLIB_VERSION/samtools-$HTSLIB_VERSION.tar.bz2 && \
    test "`md5sum samtools-$HTSLIB_VERSION.tar.bz2 | awk -v FS='  ' '{print $1}'` = $HTSLIB_VERSION_MD5" && \
    tar xjf samtools-$HTSLIB_VERSION.tar.bz2 && \
    cd samtools-$HTSLIB_VERSION && \
    ./configure --quiet --enable-libcurl --prefix=$HOME/samtools && \
    make -s -j $(nproc --all) all all-htslib && \
    make -s install install-htslib && \
    cd .. && rm -r samtools-$HTSLIB_VERSION.tar.bz2 samtools-$HTSLIB_VERSION
RUN samtools --version
