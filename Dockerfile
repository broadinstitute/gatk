ARG BASE_DOCKER=us.gcr.io/broad-dsde-methods/gatk-base-image-staging-area:gatkbase-3.3.0rc1

# stage 1 for constructing the GATK zip
FROM ${BASE_DOCKER} AS gradleBuild
LABEL stage=gatkIntermediateBuildImage
ARG RELEASE=false

RUN ls .
ADD . /gatk
WORKDIR /gatk

# Get an updated gcloud signing key, in case the one in the base image has expired
RUN rm /etc/apt/sources.list.d/google-cloud-sdk.list && \
    apt update &&\
    apt-key list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    add-apt-repository universe && apt update && \
    apt-get --assume-yes install git-lfs && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt/lists/*
RUN git lfs install --force

#Download only resources required for the build, not for testing
RUN git lfs pull --include src/main/resources/large

RUN export GRADLE_OPTS="-Xmx4048m -Dorg.gradle.daemon=false" && /gatk/gradlew clean collectBundleIntoDir shadowTestClassJar shadowTestJar -Drelease=$RELEASE
RUN cp -r $( find /gatk/build -name "*bundle-files-collected" )/ /gatk/unzippedJar/
RUN unzip -o -j $( find /gatk/unzippedJar -name "gatkPython*.zip" ) -d /gatk/unzippedJar/scripts
RUN chmod -R a+rw /gatk/unzippedJar

FROM ${BASE_DOCKER}

RUN rm /etc/apt/sources.list.d/google-cloud-sdk.list && \
    apt update && \
    apt-key list && \
    apt-get -y clean  && \
    apt-get -y autoclean  && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /gatk

RUN chmod -R a+rw /gatk
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
# "export GATK_DOCKER_CONTAINER=true" is used to allow tests to determine when the're running on the docker
# (some negative python tests use this to throw skip exceptions). See GATKBaseTest::isGATKDockerContainer.
RUN echo "source activate gatk" > /root/run_unit_tests.sh && \
    echo "export GATK_DOCKER_CONTAINER=true" >> /root/run_unit_tests.sh && \
    echo "export TEST_JAR=\$( find /jars -name \"gatk*test.jar\" )" >> /root/run_unit_tests.sh && \
    echo "export TEST_DEPENDENCY_JAR=\$( find /jars -name \"gatk*testDependencies.jar\" )" >> /root/run_unit_tests.sh && \
    echo "export GATK_JAR=$( find /gatk -name "gatk*local.jar" )" >> /root/run_unit_tests.sh && \
    echo "export GATK_LAUNCH_SCRIPT=/gatk/gatk" >> /root/run_unit_tests.sh && \
    echo "mkdir /gatk/srcdir" >> /root/run_unit_tests.sh && \
    echo "cp -rp /gatkCloneMountPoint/src/main/java/* /gatk/srcdir" >> /root/run_unit_tests.sh && \
    echo "export SOURCE_DIR=/gatk/srcdir" >> /root/run_unit_tests.sh && \
    echo "export GRADLE_OPTS=\"-Xmx1024m -Dorg.gradle.daemon=false --add-opens java.prefs/java.util.prefs=ALL-UNNAMED\"" >> /root/run_unit_tests.sh && \
    echo "export CP_DIR=/gatk/testClasses" >> /root/run_unit_tests.sh && \
    echo "ln -s /gatkCloneMountPoint/src/ /gatkCloneMountPoint/scripts/docker/src" >> /root/run_unit_tests.sh && \
    echo "ln -s /gatkCloneMountPoint/build/ /gatkCloneMountPoint/scripts/docker/build" >> /root/run_unit_tests.sh && \
    echo "cd /gatk/ && /gatkCloneMountPoint/gradlew -Dfile.encoding=UTF-8 -b /gatkCloneMountPoint/dockertest.gradle testOnPackagedReleaseJar jacocoTestReportOnPackagedReleaseJar -a -p /gatkCloneMountPoint" >> /root/run_unit_tests.sh

WORKDIR /root
RUN cp -r /root/run_unit_tests.sh /gatk
RUN cp -r gatk.jar /gatk
ENV CLASSPATH /gatk/gatk.jar:$CLASSPATH

# Start GATK Python environment

WORKDIR /gatk
ENV PATH $CONDA_PATH/envs/gatk/bin:$CONDA_PATH/bin:$PATH
RUN conda env create -vv -n gatk -f /gatk/gatkcondaenv.yml && \
    echo "source activate gatk" >> /gatk/gatkenv.rc && \
    echo "source /gatk/gatk-completion.sh" >> /gatk/gatkenv.rc && \
    conda clean -afy && \
    rm -rf /root/.cache/pip

CMD ["bash", "--init-file", "/gatk/gatkenv.rc"]

# End GATK Python environment

ENV PATH /gatk:$PATH
