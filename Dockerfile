# Using OpenJDK 8
FROM us.gcr.io/broad-gotc-prod/gatk:gatkbase-1.3.0

WORKDIR /gatk/scripts

ARG PYTHON_ZIP_PATH
ARG CONDA_ENV_PATH

ADD $PYTHON_ZIP_PATH /gatk/scripts/gatkPythonPackageArchive.zip
ADD $CONDA_ENV_PATH /gatk/scripts/gatkcondaenv.yml

RUN unzip -o -j /gatk/scripts/$PYTHON_ZIP_PATH -d /gatk/scripts && \
    conda-env create -n gatk -f /gatk/scripts/gatkcondaenv.yml && \
    echo "source activate gatk" >> /gatk/gatkenv.rc && \
    conda clean -y -all

WORKDIR /gatk

ARG GATK_JAR_PATH
ARG GATK_SPARK_JAR_PATH

ADD $GATK_JAR_PATH /gatk/gatk.jar
ADD $GATK_SPARK_JAR_PATH /gatk/gatk-spark.jar

ENV PATH /gatk:$PATH
ENTRYPOINT ["bash", "--init-file", "/gatk/gatkenv.rc"]
