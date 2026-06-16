# This Dockerfile creates a "build-base" image with tools and libraries required to build the tools and libraries used
# in the Genomic Variant Store pipeline. The Alpine version of the Google Cloud SDK is used as the base image which is
# not only the most compact of the Google Cloud SDK Docker images, but is also the image currently used by Cromwell for
# (de)localization of files in Google Cloud Storage. Sharing the base image with Cromwell's GCS localization should
# result in reuse of a cached copy of this base (and by far largest) image layer when running GVS pipelines in Terra /
# Cromwell.
#
# Because this is an Alpine-based image it is more bare-bones than its Debian-based peers. Several tools used by GVS
# (htslib, bcftools, vcftools, bedtools) are not available as Alpine packages and must be compiled from source.
# Compiling these tools makes this a moderately expensive image to create. Since this image isn't expected to change
# often it's broken out into a separate "build-base" image that can effectively be globally cached and referenced from
# the main Dockerfile.
#
# Note: pyarrow was previously built from source here because Alpine's musl libc was incompatible with PyPI's
# manylinux wheels. Apache Arrow 24+ ships musllinux wheels on PyPI, so pyarrow is now installed via pip in the
# main Dockerfile's requirements.txt instead.
FROM gcr.io/google.com/cloudsdktool/cloud-sdk:524.0.0-alpine

RUN apk update && apk upgrade

# Add all required build tools. These are not added to the main stage as they are only needed at compile time.
RUN apk add autoconf bash g++ gcc make python3-dev git openssl-dev zlib-dev xz-dev bzip2-dev curl-dev
RUN python3 -m venv /localvenv && . /localvenv/bin/activate && python3 -m ensurepip --upgrade


ARG HTSLIB_VERSION=1.23.1
RUN mkdir /htslib /htslib-build && \
    cd /htslib-build && \
    curl -L -O https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    bunzip2 htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar xf htslib-${HTSLIB_VERSION}.tar && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure --enable-libcurl --enable-gcs --prefix=/htslib && \
    make && \
    make install && \
    cd / && \
    rm -rf /htslib-build


ARG BCFTOOLS_VERSION=1.23.1
RUN mkdir /bcftools /bcftools-build && \
    cd /bcftools-build && \
    curl -L -O https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    bunzip2 bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar xf bcftools-${BCFTOOLS_VERSION}.tar && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    ./configure --prefix /bcftools --with-htslib=/htslib && \
    make && \
    make install && \
    cd / && \
    rm -rf /bcftools-build

# Build vcf-tools following this example:
# https://github.com/overcookedfrog/vcftools/blob/master/Dockerfile
ARG VCFTOOLS_VERSION=0.1.17
RUN mkdir /vcftools /vcftools-build && \
    cd /vcftools-build && \
    curl -L -O https://github.com/vcftools/vcftools/releases/download/v$VCFTOOLS_VERSION/vcftools-$VCFTOOLS_VERSION.tar.gz && \
    tar zxf vcftools-$VCFTOOLS_VERSION.tar.gz && \
    cd vcftools-$VCFTOOLS_VERSION && \
    ./configure --prefix=/vcftools && \
    make && \
    make install && \
    cd / && \
    rm -rf /vcftools-build

ARG BEDTOOLS_VERSION=2.31.1
RUN mkdir -p /bedtools /bedtools-build/src && \
    cd /bedtools-build && \
    curl -L -O https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz && \
    tar -xzf bedtools-${BEDTOOLS_VERSION}.tar.gz -C src --strip-components=1 && \
    cd src && \
    make CXXFLAGS="-g -Wall -O2 -std=c++11 -include cstdint" && \
    make prefix=/bedtools install && \
    cd / && \
    rm -rf /bedtools-build

ENV PERL5LIB="/vcftools/share/perl5/site_perl/"
