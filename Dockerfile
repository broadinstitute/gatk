# Using OpenJDK 8
FROM broadinstitute/gatk:gatkbase-1.2.1 AS builder
ARG DRELEASE

ADD . /gatk

WORKDIR /gatk
RUN /gatk/gradlew clean compileTestJava installAll localJar -Drelease=$DRELEASE

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
RUN echo "cd /gatk/ && ./gradlew jacocoTestReport" >/root/run_unit_tests.sh

WORKDIR /root

FROM broadinstitute/gatk:gatkbase-1.2.1
RUN mkdir /gatk
WORKDIR /gatk
COPY --from=builder /root/run_unit_tests.sh .
COPY --from=builder /gatk/gatk-launch .
COPY --from=builder /gatk/build/libs/gatk.jar .
COPY --from=builder /gatk/build/libs/gatk-spark.jar .
COPY --from=builder /root/install_R_packages.R .
ENV GATK_LOCAL_JAR /gatk/gatk.jar
ENV GATK_SPARK_JAR /gatk/gatk-spark.jar
WORKDIR /gatk
