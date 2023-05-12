# This Dockerfile creates a "build-base" image with tools and libraries required to build the tools and libraries used
# in the Genomic Variant Store pipeline. The Alpine version of the Google Cloud SDK is used as the base image which is
# not only the most compact of the Google Cloud SDK Docker images, but is also the image currently used by Cromwell for
# (de)localization of files in Google Cloud Storage. Sharing the base image with Cromwell's GCS localization should
# result in reuse of a cached copy of this base (and by far largest) image layer when running GVS pipelines in Terra /
# Cromwell.
#
# Because this is an Alpine-based image it is more bare-bones than its Debian-based peers. Key components missing here
# are the Apache Arrow library (a requirement for pyarrow which in turn is a requirement for the google-cloud-bigquery
# Python module) and bcftools. Compiling all these tools makes this a fairly expensive image to create (an hour or so
# under ideal circumstances, potentially much longer on low memory and/or non-x86 build hosts). Since this image isn't
# expected to change often it's broken out into a separate "build-base" image that can effectively be globally cached
# and referenced from the main Dockerfile.
FROM gcr.io/google.com/cloudsdktool/cloud-sdk:426.0.0-alpine

RUN apk update && apk upgrade
RUN python3 -m ensurepip --upgrade

# Add all required build tools. These will not be added to the main stage as they are only required to build PyArrow
# and bcftools but not to use them.
RUN apk add autoconf bash cmake g++ gcc make ninja python3-dev git openssl-dev zlib-dev xz-dev bzip2-dev curl-dev

# Unfortunately neither pyarrow nor google-cloud-bigquery will fetch or build Apache Arrow when `pip install`ed from
# this base image. Therefore we do the Apache Arrow build ourselves. In order to keep the final image size small this
# Dockerfile is set up to do a multi-stage build following the usual pattern of "build" stage / "main" stage.
# https://docs.docker.com/build/building/multi-stage/#use-multi-stage-builds
#
# The build stage installs the required development tools, downloads the Apache Arrow source bundle and builds all
# required components including Apache Arrow C++ libraries, pyarrow Python module, and all pyarrow dependencies
# including the numpy Python module. The main stage will then use the same base image and copy over the artifacts
# produced by the build stage without having to install development tools or clean up after a build.

ARG ARROW_VERSION=11.0.0
RUN cd / && \
    curl -O https://dlcdn.apache.org/arrow/arrow-$ARROW_VERSION/apache-arrow-$ARROW_VERSION.tar.gz && \
    tar xfz apache-arrow-$ARROW_VERSION.tar.gz

# Pyarrow build instructions from https://arrow.apache.org/docs/developers/python.html#python-development
# Modified slightly for the requirements of this installation:
# - Download a static source tarball rather than cloning the git repo.
# - Use `ninja` to build the C++ libraries as the `make` system doesn't seem to work as of Arrow 10.0.0.
# - Install PyArrow and its dependencies specifying the --user flag so all artifacts go to the /root/.local directory
#   which can easily be copied to the main stage below.
ARG ARROW_SRC_DIR=/apache-arrow-$ARROW_VERSION
RUN pip3 install --user -r $ARROW_SRC_DIR/python/requirements-build.txt

RUN mkdir /dist
RUN mkdir $ARROW_SRC_DIR/cpp/build && \
    cd $ARROW_SRC_DIR/cpp/build && \
    cmake .. --preset ninja-release-python && \
    cmake --build . && \
    cmake --install .

ARG PYARROW_WITH_PARQUET=1
ARG PYARROW_WITH_DATASET=1
ARG PYARROW_PARALLEL=4
RUN cd $ARROW_SRC_DIR/python && \
    python3 setup.py build_ext --inplace && \
    pip3 install wheel && \
    python3 setup.py build_ext --build-type=release \
              --bundle-arrow-cpp bdist_wheel && \
    pip3 install --user /apache-arrow-$ARROW_VERSION/python/dist/pyarrow-$ARROW_VERSION-*.whl

# Straightforward bcftools build following these instructions:
# https://github.com/samtools/bcftools/blob/develop/INSTALL
ARG BCFTOOLS_VERSION=1.17
RUN mkdir /bcftools bcftools-build && \
    cd bcftools-build && \
    git clone --recurse-submodules https://github.com/samtools/htslib.git && \
    git clone https://github.com/samtools/bcftools.git && \
    cd bcftools && \
    git checkout tags/$BCFTOOLS_VERSION -b $BCFTOOLS_VERSION && \
    autoheader && \
    autoconf && \
    ./configure --prefix /bcftools && \
    make && \
    make install
