# Using OpenJDK 8
# Location of the unzipped gatk bundle files
ARG BASEIMAGE
FROM $BASEIMAGE

# Location of the unzipped conda files
RUN mkdir /condaenv

ADD scripts/gatkcondaenv.yml.template /condaenv/gatkcondaenv.yml.template

RUN apt-get install -y python && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove

# Start GATK Python environment

ENV DOWNLOAD_DIR /downloads
ENV CONDA_URL https://repo.continuum.io/miniconda/Miniconda3-4.3.30-Linux-x86_64.sh
ENV CONDA_MD5 = "0b80a152332a4ce5250f3c09589c7a81"
ENV CONDA_PATH /opt/miniconda
RUN mkdir $DOWNLOAD_DIR && \
    wget -nv -O $DOWNLOAD_DIR/miniconda.sh $CONDA_URL && \
    test "`md5sum $DOWNLOAD_DIR/miniconda.sh | awk -v FS='  ' '{print $1}'` = $CONDA_MD5" && \
    bash $DOWNLOAD_DIR/miniconda.sh -p $CONDA_PATH -b && \
    rm $DOWNLOAD_DIR/miniconda.sh
WORKDIR /gatk
ENV PATH $CONDA_PATH/envs/gatk/bin:$CONDA_PATH/bin:$PATH
RUN conda env create -n gatk -f /condaenv/gatkcondaenv.yml && \
    conda clean -y -all && \
    rm -rf /root/.cache/pip

# End GATK Python environment
