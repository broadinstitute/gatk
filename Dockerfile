# Using OpenJDK 8
FROM broadinstitute/gatk:gatkbase-1.2.3
ARG DRELEASE

ADD . /gatk

WORKDIR /gatk
RUN /gatk/gradlew clean compileTestJava sparkJar localJar condaEnvironmentDefinition -Drelease=$DRELEASE

WORKDIR /root

# Make sure we can see a help message
RUN ln -sFv /gatk/build/libs/gatk.jar
RUN java -jar gatk.jar -h

#Setup test data
WORKDIR /gatk
# Create link to where test data is expected
RUN ln -s /testdata src/test/resources

# Create a simple unit test runner
ENV CI true
RUN echo "source activate gatk" > /root/run_unit_tests.sh && \
    echo "cd /gatk/ && ./gradlew jacocoTestReport" >> /root/run_unit_tests.sh

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
WORKDIR /gatk/build
RUN conda env create -n gatk -f gatkcondaenv.yml && \
    echo "source activate gatk" >> /gatk/gatkenv.rc
WORKDIR /gatk

CMD ["bash", "--init-file", "/gatk/gatkenv.rc"]

# End GATK Python environment

WORKDIR /gatk

ENV PATH /gatk:$PATH
