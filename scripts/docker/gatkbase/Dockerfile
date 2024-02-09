# This Dockerfile does not require any files that are in the GATK4 repo.
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND noninteractive \
    HOME /root \
    JAVA_LIBRARY_PATH /usr/lib/jni \
    DOWNLOAD_DIR /downloads \
    CONDA_URL https://repo.anaconda.com/miniconda/Miniconda3-py310_23.10.0-1-Linux-x86_64.sh \
    CONDA_SHA256 "c7a34df472feb69805b64df6e8db58363c5ccab41cd3b40b07e3e6dfb924359a" \
    CONDA_PATH /opt/miniconda \
    PATH $CONDA_PATH/bin:$PATH

WORKDIR /root

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
    jq \
    git \
    gpg-agent \
    build-essential \
    openjdk-17-jdk \
    vim \
    software-properties-common && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" \
    | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg \
    | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt update -y && \
    apt install -y --no-install-recommends google-cloud-cli && \
    mkdir $DOWNLOAD_DIR && \
    wget -nv -O $DOWNLOAD_DIR/miniconda.sh $CONDA_URL && \
    test "$(sha256sum $DOWNLOAD_DIR/miniconda.sh | awk -v FS=' ' -v ORS='' '{print $1}')" = "$CONDA_SHA256" && \
    bash $DOWNLOAD_DIR/miniconda.sh -p $CONDA_PATH -b && \
    conda clean -afy && \
    conda config --set auto_update_conda false && \
    conda config --set solver libmamba && \
    apt -y clean && \
    apt -y autoclean && \
    apt -y autoremove && \
    rm -rf /var/lib/apt/lists/* /usr/lib/google-cloud-sdk/bin/anthoscli /usr/lib/google-cloud-sdk/platform/anthoscli_licenses /usr/lib/google-cloud-sdk/platform/bundledpythonunix /root/.cache/pip $DOWNLOAD_DIR && \
    find / -wholename "*__pycache__/*.pyc" -exec rm {} +

CMD ["bash"]
