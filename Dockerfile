# Using OpenJDK 8
FROM broadinstitute/gatk:gatkbase-2.0.1

# Location of the unzipped gatk bundle files
ARG ZIPPATH

ADD $ZIPPATH /gatk

WORKDIR /gatk

#Setup linked jars that may be needed for running gatk
RUN ln -s $( find /gatk -name "gatk*local.jar" ) gatk.jar
RUN ln -s $( find /gatk -name "gatk*local.jar" ) /root/gatk.jar
RUN ln -s $( find /gatk -name "gatk*spark.jar" ) gatk-spark.jar

WORKDIR /root

 # Make sure we can see a help message
RUN java -jar gatk.jar -h
RUN mkdir /gatkCloneMountPoint
RUN mkdir /jars
RUN mkdir .gradle

WORKDIR /gatk

# Create a simple unit test runner
ENV CI true
RUN echo "source activate gatk" > /root/run_unit_tests.sh && \
    echo "export TEST_JAR=\$( find /jars -name \"gatk*test.jar\" )" >> /root/run_unit_tests.sh && \
    echo "export TEST_DEPENDENCY_JAR=\$( find /jars -name \"gatk*testDependencies.jar\" )" >> /root/run_unit_tests.sh && \
    echo "export GATK_JAR=$( find /gatk -name "gatk*local.jar" )" >> /root/run_unit_tests.sh && \
    echo "cp -rp /gatkCloneMountPoint/src/main/java/* /gatk/srcdir" >> /root/run_unit_tests.sh && \
    echo "export SOURCE_DIR=/gatk/srcdir" >> /root/run_unit_tests.sh && \
    echo "export GRADLE_OPTS=\"-Xmx1024m -Dorg.gradle.daemon=false\"" /root/run_unit_tests.sh && \
    echo "export CP_DIR=/gatk/testClasses" /root/run_unit_tests.sh && \
    echo "ln -s /gatkCloneMountPoint/src/ /gatkCloneMountPoint/scripts/docker/src" >> /root/run_unit_tests.sh && \
    echo "ln -s /gatkCloneMountPoint/build/ /gatkCloneMountPoint/scripts/docker/build" >> /root/run_unit_tests.sh && \
    echo "cd /gatk/ && /gatkCloneMountPoint/gradlew -b /gatkCloneMountPoint/dockertest.gradle testOnPackagedReleaseJar jacocoTestReportOnPackagedReleaseJar -a -p /gatkCloneMountPoint" >> /root/run_unit_tests.sh

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
WORKDIR /gatk
ENV PATH $CONDA_PATH/envs/gatk/bin:$CONDA_PATH/bin:$PATH
RUN conda env create -n gatk -f /gatk/gatkcondaenv.yml && \
    echo "source activate gatk" >> /gatk/gatkenv.rc && \
    conda clean -y -all && \
    rm -rf /root/.cache/pip

CMD ["bash", "--init-file", "/gatk/gatkenv.rc"]

# End GATK Python environment

ENV PATH /gatk:$PATH
