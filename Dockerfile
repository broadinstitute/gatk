# stage 1 for constructing the GATK zip
FROM broadinstitute/gatk:gatkbase-2.3.0 AS gradleBuild
LABEL stage=gatkIntermediateBuildImage
ARG RELEASE=false

RUN ls .
ADD . /gatk
WORKDIR /gatk

# Get an updated gcloud signing key, in case the one in the base image has expired
RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add -
RUN add-apt-repository universe && apt update
RUN apt-get --assume-yes install git-lfs
RUN git lfs install --force

RUN git lfs pull

RUN export GRADLE_OPTS="-Xmx4048m -Dorg.gradle.daemon=false" && /gatk/gradlew clean collectBundleIntoDir shadowTestClassJar shadowTestJar -Drelease=$RELEASE
RUN cp -r $( find /gatk/build -name "*bundle-files-collected" )/ /gatk/unzippedJar/
RUN unzip -o -j $( find /gatk/unzippedJar -name "gatkPython*.zip" ) -d /gatk/unzippedJar/scripts

# Using OpenJDK 8
FROM broadinstitute/gatk:gatkbase-2.3.0

WORKDIR /gatk

# Location of the unzipped gatk bundle files
COPY --from=gradleBuild /gatk/unzippedJar .

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
    echo "export GATK_DOCKER_CONTAINER=true" >> /root/run_unit_tests.sh && \
    echo "export TEST_JAR=\$( find /jars -name \"gatk*test.jar\" )" >> /root/run_unit_tests.sh && \
    echo "export TEST_DEPENDENCY_JAR=\$( find /jars -name \"gatk*testDependencies.jar\" )" >> /root/run_unit_tests.sh && \
    echo "export GATK_JAR=$( find /gatk -name "gatk*local.jar" )" >> /root/run_unit_tests.sh && \
    echo "export GATK_LAUNCH_SCRIPT=/gatk/gatk" >> /root/run_unit_tests.sh && \
    echo "mkdir /gatk/srcdir" >> /root/run_unit_tests.sh && \
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
ENV CLASSPATH /gatk/gatk.jar:$CLASSPATH

# Start GATK Python environment

WORKDIR /gatk
RUN chmod -R a+rw /gatk
ENV PATH $CONDA_PATH/envs/gatk/bin:$CONDA_PATH/bin:$PATH
RUN conda env create -n gatk -f /gatk/gatkcondaenv.yml && \
    echo "source activate gatk" >> /gatk/gatkenv.rc && \
    echo "source /gatk/gatk-completion.sh" >> /gatk/gatkenv.rc && \
    conda clean -afy && \
    find /opt/miniconda/ -follow -type f -name '*.a' -delete && \
    find /opt/miniconda/ -follow -type f -name '*.pyc' -delete && \
    rm -rf /root/.cache/pip

CMD ["bash", "--init-file", "/gatk/gatkenv.rc"]

# End GATK Python environment

ENV PATH /gatk:$PATH
