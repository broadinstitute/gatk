# Using OpenJDK 8
FROM broadinstitute/gatk:gatkbase-1.1

ADD . /gatk

WORKDIR /gatk
RUN /gatk/gradlew clean compileTestJava installAll localJar

WORKDIR /root

# Make sure we can see a help message
RUN ln -sFv /gatk/build/libs/gatk.jar
RUN java -jar gatk.jar -h

#Setup test data
WORKDIR /gatk
# remove existing test data
RUN rm -rf src/test/resources
# Create link to where test data is expeced
RUN ln -s /testdata src/test/resources

# Create a simple unit test runner
ENV CI true
RUN echo "cd /gatk/ && ./gradlew jacocoTestReport" >/root/run_unit_tests.sh

WORKDIR /root
