FROM google/cloud-sdk:472.0.0-slim

# Use the stable URL to plink2 from the pgen-jni project.
ADD https://github.com/broadinstitute/pgen-jni/raw/de500ad710ce70a8171fd3d0806857a1f36d8a2e/plink2_assets/plink2_linux_x86_64_20240318.zip plink2.zip

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get -y --no-install-recommends install \
    unzip

RUN unzip plink2.zip -d /usr/local/bin &&\
    rm plink2.zip

ENV PATH /usr/local/bin/plink2:$PATH
