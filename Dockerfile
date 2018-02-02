# Using OpenJDK 8
FROM broadinstitute/gatk:gatkbase-1.2.4-testA
ARG ZIPPATH

ADD $ZIPPATH /gatk

WORKDIR /gatk
RUN ln -s /gatk/$( find . -name "gatk*local.jar" ) gatk.jar
RUN ln -s /gatk/$( find . -name "gatk*spark.jar" ) gatk-spark.jar
# RUN /gatk/gradlew clean compileTestJava sparkJar localJar createPythonPackageArchive -Drelease=$DRELEASE

WORKDIR /root

# Make sure we can see a help message
RUN ln -sFv /gatk/gatk.jar
RUN mkdir /gatksrc
RUN mkdir .gradle

#Setup test data
WORKDIR /gatk

# Create a simple unit test runner
ENV CI true
RUN echo "source activate gatk" > /root/run_unit_tests.sh && \
    echo "export DOCKER_TEST=\"true\"" >> /root/run_unit_tests.sh && \
    echo "cd /gatk/ && /gatksrc/gradlew jacocoTestReport -a -p /gatksrc " >> /root/run_unit_tests.sh

WORKDIR /root
RUN cp -r /root/run_unit_tests.sh /gatk
RUN cp -r gatk.jar /gatk
RUN cp -r install_R_packages.R /gatk

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
ENV PATH $CONDA_PATH/envs/gatk/bin:$CONDA_PATH/bin:$PATH
RUN mv /gatk/gatkcondaenv.yml /gatk/scripts
WORKDIR /gatk
RUN conda-env create -n gatk -f /gatk/scripts/gatkcondaenv.yml && \
    echo "source activate gatk" >> /gatk/gatkenv.rc

CMD ["bash", "--init-file", "/gatk/gatkenv.rc"]

RUN conda clean -y -all
# End GATK Python environment

WORKDIR /gatk

ENV PATH /gatk:$PATH
