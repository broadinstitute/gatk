# The image produced by this Dockerfile contains tools and libraries required to support the Genomic Variant Store
# pipeline. The Alpine version of the Google Cloud SDK is used as the base image which is not only the most compact of
# the Google Cloud SDK Docker images, but is also the image currently used by Cromwell for (de)localization of files in
# Google Cloud Storage. Sharing the base image with Cromwell's GCS localization should result in reuse of a cached copy
# of this base (and by far largest) image layer when running GVS pipelines in Terra / Cromwell.
#
# Because this is an Alpine-based image it is more bare-bones than its Debian-based peers. A key component missing here
# is the Apache Arrow library which is a requirement for pyarrow which is a requirement for the google-cloud-bigquery
# Python module.
#
# Unfortunately neither pyarrow nor google-cloud-bigquery will fetch or build Apache Arrow when `pip install`ed from
# this base image. Therefore we do the Apache Arrow build ourselves. In order to keep the final image size small this
# Dockerfile is set up to do a multi-stage build following the usual pattern of "build" stage / "main" stage.
# https://docs.docker.com/build/building/multi-stage/#use-multi-stage-builds
#
# The build stage installs the required development tools, downloads the Apache Arrow source bundle and builds all
# required components including Apache Arrow C++ libraries, pyarrow Python module, and all pyarrow dependencies
# including the numpy Python module. The main stage will then use the same base image and copy over the artifacts
# produced by the build stage without having to install development tools or clean up after a build.
FROM gcr.io/google.com/cloudsdktool/cloud-sdk:404.0.0-alpine as build

RUN apk update && apk upgrade
RUN python3 -m ensurepip --upgrade

# Add all required build tools. These will not be added to the main stage as they are only required to build PyArrow
# but not to use it.
RUN apk add autoconf bash cmake g++ gcc make ninja python3-dev git openssl-dev

# Build Apache Arrow version 8.0.0 as version 9.0.0 does not compile under Alpine:
# https://issues.apache.org/jira/browse/ARROW-17329
ARG ARROW_VERSION=8.0.0
RUN cd / && \
    curl -O https://dlcdn.apache.org/arrow/arrow-$ARROW_VERSION/apache-arrow-$ARROW_VERSION.tar.gz && \
    tar xfz apache-arrow-$ARROW_VERSION.tar.gz

# Pyarrow build instructions from https://arrow.apache.org/docs/developers/python.html#python-development
# Modified slightly for the requirements of this installation:
# - Download a static source tarball rather than cloning the git repo.
# - Use build type Release rather than Debug.
# - Do not build tests.
# - Install PyArrow and its dependencies specifying the --user flag so all artifacts go to the /root/.local directory
#   which can easily be copied to the main stage below.
RUN pip3 install --user -r /apache-arrow-$ARROW_VERSION/python/requirements-build.txt

RUN mkdir /dist
ARG ARROW_HOME=/dist
ARG LD_LIBRARY_PATH=/dist/lib:$LD_LIBRARY_PATH
RUN mkdir /apache-arrow-$ARROW_VERSION/cpp/build && \
    cd /apache-arrow-$ARROW_VERSION/cpp/build && \
    cmake -DCMAKE_INSTALL_PREFIX=$ARROW_HOME \
          -DCMAKE_INSTALL_LIBDIR=lib \
          -DCMAKE_BUILD_TYPE=Release \
          -DARROW_DATASET=ON \
          -DARROW_WITH_BZ2=ON \
          -DARROW_WITH_ZLIB=ON \
          -DARROW_WITH_ZSTD=ON \
          -DARROW_WITH_LZ4=ON \
          -DARROW_WITH_SNAPPY=ON \
          -DARROW_WITH_BROTLI=ON \
          -DARROW_PARQUET=ON \
          -DPARQUET_REQUIRE_ENCRYPTION=ON \
          -DARROW_PYTHON=ON \
          -DARROW_BUILD_TESTS=OFF \
          .. && \
    make -j4 && \
    make install

ARG PYARROW_WITH_PARQUET=1
ARG PYARROW_WITH_DATASET=1
ARG PYARROW_PARALLEL=4
RUN cd /apache-arrow-$ARROW_VERSION/python && \
    python3 setup.py build_ext --inplace && \
    pip3 install wheel && \
    python3 setup.py build_ext --build-type=release \
              --bundle-arrow-cpp bdist_wheel && \
    pip3 install --user /apache-arrow-$ARROW_VERSION/python/dist/pyarrow-$ARROW_VERSION-*.whl

# Install all of our variantstore Python requirements.
COPY requirements.txt requirements.txt
RUN pip3 install --user -r requirements.txt

# The main layer does not install development tools, instead copies artifacts from the build layer above.
FROM gcr.io/google.com/cloudsdktool/cloud-sdk:404.0.0-alpine as main

RUN apk update && apk upgrade
RUN python3 -m ensurepip --upgrade

# Add any generally useful Alpine packages here.
RUN apk add --no-cache jq

# The build stage generated Python artifacts to /root/.local via `pip install --user`, so grab all of those.
COPY --from=build /root/.local /root/.local

# Copy the application source code.
RUN mkdir /app
COPY *.py /app
COPY *.sql /app

# Copy the schema files needed for VAT creation.
COPY *.json /data/variant_annotation_table/schema/

WORKDIR /app
