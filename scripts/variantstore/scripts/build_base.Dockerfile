# This Dockerfile creates a "build-base" image with tools and libraries required to build the tools and libraries used
# in the Genomic Variant Store pipeline. The Alpine version of the Google Cloud SDK is used as the base image as it is
# the most compact of the Google Cloud SDK Docker images.
#
# Because this is an Alpine-based image it is more bare-bones than its Debian-based peers. Several tools used by GVS
# (htslib, bcftools, vcftools) are not available as Alpine packages and must be compiled from source.
# Compiling these tools makes this a moderately expensive image to create. Since this image isn't expected to change
# often it's broken out into a separate "build-base" image that can effectively be globally cached and referenced from
# the main Dockerfile.
#
# Note: pyarrow was previously built from source here because Alpine's musl libc was incompatible with PyPI's
# manylinux wheels. It has since been removed from the image entirely as no production code requires it.
#
# IMPORTANT: bumping the tag below only changes what this Dockerfile produces; it has no effect until the resulting
# image is rebuilt and pushed via `build_build_base_docker.sh` and the new tag is wired into Dockerfile's `build`
# stage FROM line. Until that happens, Dockerfile's `build` and `main` stages are on different Python versions, which
# silently breaks the venv copied from `build` into `main` (its `python3` symlink dangles, so `python3` falls back to
# the bare system interpreter with none of the pip-installed packages -- e.g. `ModuleNotFoundError: No module named
# 'google'`). Keep the version here and in Dockerfile's `build` stage FROM in lockstep.
FROM gcr.io/google.com/cloudsdktool/cloud-sdk:582.0.0-alpine

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


ENV PERL5LIB="/vcftools/share/perl5/site_perl/"
