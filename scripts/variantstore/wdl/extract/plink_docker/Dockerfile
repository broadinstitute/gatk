FROM google/cloud-sdk:448.0.0

ADD https://s3.amazonaws.com/plink2-assets/alpha4/plink2_linux_avx2_20230922.zip plink2.zip

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get -y --no-install-recommends install \
    unzip

RUN unzip plink2.zip -d /usr/local/bin &&\
    rm plink2.zip

ENV PATH /usr/local/bin/plink2:$PATH