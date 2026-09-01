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
# `python3` is listed explicitly (it's also a transitive dependency of python3-dev) because this image's cloud-sdk
# base bundles its own private Python interpreter for gcloud's internal use, at /usr/local/bin, ahead of Alpine's own
# apk-managed /usr/bin/python3 on PATH. That bundled interpreter's version moves independently with every cloud-sdk
# release (it jumped 3.12 -> 3.14 between the 524.0.0 and 582.0.0 tags) and is not meant to be used outside gcloud
# itself, so relying on a bare `python3` here would silently tie our venv to whatever Python gcloud happens to embed.
# Alpine's own python3 package tracks the (much more slowly moving) Alpine release instead, so the venv below is
# built against /usr/bin/python3 explicitly. Dockerfile's `main` stage must also `apk add python3` for the same
# reason: this venv is copied there via `COPY --from=build`, and needs the matching interpreter/shared libs present
# at the same path to not be a dangling symlink. See the notes in Dockerfile.
RUN apk add autoconf bash g++ gcc make python3 python3-dev git openssl-dev zlib-dev xz-dev bzip2-dev curl-dev
RUN /usr/bin/python3 -m venv /localvenv && . /localvenv/bin/activate && python3 -m ensurepip --upgrade


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
