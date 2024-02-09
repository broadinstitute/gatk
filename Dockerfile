ARG BASE_DOCKER=mattmcl.azurecr.io/gatk:latest
FROM ${BASE_DOCKER} AS gradleBuild
LABEL stage=gatkIntermediateBuildImage
ARG RELEASE=false

RUN ls . && \
    ADD . /gatk && \
    WORKDIR /gatk && \
    rm /etc/apt/sources.list.d/google-cloud-sdk.list && \
    apt update && \
    apt-key list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    add-apt-repository universe && apt update && \
    apt-get --assume-yes install git-lfs && \
    apt-get -y clean && \
    apt-get -y autoclean && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt/lists/* && \
    git lfs install --force && \
    git lfs pull --include src/main/resources/large && \
    export GRADLE_OPTS="-Xmx4048m -Dorg.gradle.daemon=false" && /gatk/gradlew clean collectBundleIntoDir shadowTestClassJar shadowTestJar -Drelease=$RELEASE && \
    cp -r $( find /gatk/build -name "*bundle-files-collected" )/ /gatk/unzippedJar/ && \
    unzip -o -j $( find /gatk/unzippedJar -name "gatkPython*.zip" ) -d /gatk/unzippedJar/scripts && \
    chmod -R a+rw /gatk/unzippedJar

FROM ${BASE_DOCKER}

RUN rm /etc/apt/sources.list.d/google-cloud-sdk.list && \
    apt update && \
    apt-key list && \
    apt-get -y clean && \
    apt-get -y autoclean && \
    apt-get -y autoremove && \
    rm -rf /var/lib/apt/lists/* && \
    WORKDIR /gatk && \
    chmod -R a+rw /gatk && \
    COPY --from=gradleBuild /gatk/unzippedJar . && \
    ln -s $( find /gatk -name "gatk*local.jar" ) gatk.jar && \
    ln -s $( find /gatk -name "gatk*local.jar" ) /root/gatk.jar && \
    ln -s $( find /gatk -name "gatk*spark.jar" ) gatk-spark.jar && \
    WORKDIR /root && \
    java -jar gatk.jar -h && \
    mkdir /gatkCloneMountPoint /jars .gradle && \
    WORKDIR /gatk && \
    echo "source activate gatk" > /root/run_unit_tests.sh && \
    echo "export GATK_DOCKER_CONTAINER=true" >> /root/run_unit_tests.sh && \
    echo "export TEST_JAR=\$( find /jars -name \"gatk*test.jar\" )" >> /root/run_unit_tests.sh && \
    echo "export TEST_DEPENDENCY_JAR=\$( find /jars -name \"gatk*testDependencies.jar\" )" >> /root/run_unit_tests.sh && \
    echo "export GATK_JAR=$( find /gatk -name \"gatk*local.jar\" )" >> /root/run_unit_tests.sh && \
    echo "export GATK_LAUNCH_SCRIPT=/gatk/gatk" >> /root/run_unit_tests.sh && \
    echo "mkdir /gatk/srcdir" >> /root/run_unit_tests.sh && \
    echo "cp -rp /gatkCloneMountPoint/src/main/java/* /gatk/srcdir" >> /root/run_unit_tests.sh && \
    echo "export SOURCE_DIR=/gatk/srcdir" >> /root/run_unit_tests.sh && \
    echo "export GRADLE_OPTS=\"-Xmx1024m -Dorg.gradle.daemon=false --add-opens java.prefs/java.util.prefs=ALL-UNNAMED\"" >> /root/run_unit_tests.sh && \
    echo "export CP_DIR=/gatk/testClasses" >> /root/run_unit_tests.sh && \
    echo "ln -s /gatkCloneMountPoint/src/ /gatkCloneMountPoint/scripts/docker/src" >> /root/run_unit_tests.sh && \
    echo "ln -s /gatkCloneMountPoint/build/ /gatkCloneMountPoint/scripts/docker/build" >> /root/run_unit_tests.sh && \
    echo "cd /gatk/ && /gatkCloneMountPoint/gradlew -Dfile.encoding=UTF-8 -b /gatkCloneMountPoint/dockertest.gradle testOnPackagedReleaseJar jacocoTestReportOnPackagedReleaseJar -a -p /gatkCloneMountPoint" >> /root/run_unit_tests.sh && \
    WORKDIR /root && \
    cp -r /root/run_unit_tests.sh /gatk && \
    cp -r gatk.jar /gatk && \
    ENV CLASSPATH /gatk/gatk.jar:$CLASSPATH

ENV PATH /gatk:$PATH
WORKDIR /gatk
RUN conda env create -n gatk -f /gatk/gatkcondaenv.yml && \
    echo "source activate gatk" >> /gatk/gatkenv.rc && \
    echo "source /gatk/gatk-completion.sh" >> /gatk/gatkenv.rc && \
    conda clean -afy && \
    rm -rf /root/.cache/pip

CMD ["bash", "--init-file", "/gatk/gatkenv.rc"]
