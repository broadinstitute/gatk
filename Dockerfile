# Using OpenJDK 8
FROM broadinstitute/gatk:gatkbase-2.2.0

ARG UNAME=gatk
ARG UID=1000
ARG GID=1000
# Location of the unzipped gatk bundle files
ARG ZIPPATH

#RUN cp -r install_R_packages.R /gatk

ADD $ZIPPATH /gatk
RUN mkdir /gatk/gatkCloneMountPoint && \
    mkdir /gatk/jars && \
    mkdir /gatk/.gradle

RUN groupadd -g $GID -o $UNAME
RUN useradd -b / -u $UID -g $GID -o -s /bin/bash $UNAME
RUN chown -R $UID:$GID /gatk
USER $UNAME
ENV HOME /gatk

WORKDIR /gatk

#Setup linked jars that may be needed for running gatk
RUN ln -s $( find /gatk -name "gatk*local.jar" ) gatk.jar
RUN ln -s $( find /gatk -name "gatk*spark.jar" ) gatk-spark.jar

# Make sure we can see a help message
RUN java -jar gatk.jar -h

# Create a simple unit test runner
ENV CI true
RUN echo "source activate gatk" > /gatk/run_unit_tests.sh && \
    echo "export GATK_DOCKER_CONTAINER=true" >> /gatk/run_unit_tests.sh && \
    echo "export TEST_JAR=\$( find /gatk/jars -name \"gatk*test.jar\" )" >> /gatk/run_unit_tests.sh && \
    echo "export TEST_DEPENDENCY_JAR=\$( find /gatk/jars -name \"gatk*testDependencies.jar\" )" >> /gatk/run_unit_tests.sh && \
    echo "export GATK_JAR=$( find /gatk -name "gatk*local.jar" )" >> /gatk/run_unit_tests.sh && \
    echo "mkdir /gatk/srcdir" >> /gatk/run_unit_tests.sh && \
    echo "cp -rp /gatk/gatkCloneMountPoint/src/main/java/* /gatk/srcdir" >> /gatk/run_unit_tests.sh && \
    echo "export SOURCE_DIR=/gatk/srcdir" >> /gatk/run_unit_tests.sh && \
    echo "export GRADLE_OPTS=\"-Xmx1024m -Dorg.gradle.daemon=false\"" /gatk/run_unit_tests.sh && \
    echo "export CP_DIR=/gatk/testClasses" /gatk/run_unit_tests.sh && \
    echo "ln -s /gatk/gatkCloneMountPoint/src/ /gatk/gatkCloneMountPoint/scripts/docker/src" >> /gatk/run_unit_tests.sh && \
    echo "ln -s /gatk/gatkCloneMountPoint/build/ /gatk/gatkCloneMountPoint/scripts/docker/build" >> /gatk/run_unit_tests.sh && \
    echo "cd /gatk/ && /gatk/gatkCloneMountPoint/gradlew -b /gatk/gatkCloneMountPoint/dockertest.gradle testOnPackagedReleaseJar jacocoTestReportOnPackagedReleaseJar -a -p /gatk/gatkCloneMountPoint" >> /gatk/run_unit_tests.sh

ENV CLASSPATH /gatk/gatk.jar:$CLASSPATH

# Start GATK Python environment

ENV DOWNLOAD_DIR /gatk/downloads
ENV CONDA_URL https://repo.continuum.io/miniconda/Miniconda3-4.3.30-Linux-x86_64.sh
ENV CONDA_MD5 "0b80a152332a4ce5250f3c09589c7a81"
ENV CONDA_PATH /gatk/miniconda
RUN mkdir $DOWNLOAD_DIR && \
    wget -nv -O $DOWNLOAD_DIR/miniconda.sh $CONDA_URL && \
    test "`md5sum $DOWNLOAD_DIR/miniconda.sh | awk -v FS='  ' '{print $1}'` = $CONDA_MD5" && \
    bash $DOWNLOAD_DIR/miniconda.sh -p $CONDA_PATH -b && \
    rm -rf $DOWNLOAD_DIR
ENV PATH $CONDA_PATH/envs/gatk/bin:$CONDA_PATH/bin:$PATH
RUN conda env create -n gatk -f /gatk/gatkcondaenv.yml && \
    echo "source activate gatk" >> /gatk/gatkenv.rc && \
    echo "source /gatk/gatk-completion.sh" >> /gatk/gatkenv.rc && \
    conda clean -y -all && \
    rm -rf /gatk/.cache/pip

CMD ["bash", "--init-file", "/gatk/gatkenv.rc"]

# End GATK Python environment

ENV PATH /gatk:$PATH
