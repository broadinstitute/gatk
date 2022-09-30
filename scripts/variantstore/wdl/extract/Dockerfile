# Build layer that installs all the development tools, downloads the Apache Arrow source bundle and builds the C++
# libraries, Python module, and all dependencies.
FROM gcr.io/google.com/cloudsdktool/cloud-sdk:404.0.0-alpine as build

RUN apk update && apk upgrade
RUN python3 -m ensurepip --upgrade

# Add all the build tools required for the build stage. These will not be added to the main stage as they are only
# required to build the PyArrow wheel.
RUN apk add autoconf bash cmake g++ gcc make ninja python3-dev git openssl-dev

# Build Apache Arrow version 8.0.0 as version 9.0.0 does not compile under Alpine:
# https://issues.apache.org/jira/browse/ARROW-17329
RUN cd / && \
    curl -O https://dlcdn.apache.org/arrow/arrow-8.0.0/apache-arrow-8.0.0.tar.gz && \
    tar xfz apache-arrow-8.0.0.tar.gz

# Pyarrow build instructions from https://arrow.apache.org/docs/developers/python.html#python-development
# Modified slightly for the requirements of this installation:
# - Download static source bundle rather than clone git repo.
# - Build Release rather than Debug.
# - Do not build tests.

RUN pip3 install -r /apache-arrow-8.0.0/python/requirements-build.txt

RUN mkdir /dist
ARG ARROW_HOME=/dist
ARG LD_LIBRARY_PATH=/dist/lib:$LD_LIBRARY_PATH
RUN mkdir /apache-arrow-8.0.0/cpp/build && \
    cd /apache-arrow-8.0.0/cpp/build && \
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
RUN cd /apache-arrow-8.0.0/python && \
    python3 setup.py build_ext --inplace && \
    pip3 install wheel && \
    python3 setup.py build_ext --build-type=release \
              --bundle-arrow-cpp bdist_wheel

# Main layer, does not install development tools, copies artifacts from the build layer instead.
FROM gcr.io/google.com/cloudsdktool/cloud-sdk:404.0.0-alpine as main

RUN apk update && apk upgrade
RUN python3 -m ensurepip --upgrade

RUN apk add openssl jq
# numpy is a dependency of PyArrow that was built as part of installing PyArrow's build requirements.
COPY --from=build /usr/lib/python3.9/site-packages/numpy /usr/lib/python3.9/site-packages/numpy
COPY --from=build /usr/lib/python3.9/site-packages/numpy-1.19.3.dist-info /usr/lib/python3.9/site-packages/numpy-1.19.3.dist-info
COPY --from=build /apache-arrow-8.0.0/python/dist/pyarrow-8.0.0-cp39-cp39-linux_x86_64.whl .
RUN pip3 install pyarrow-8.0.0-cp39-cp39-linux_x86_64.whl

COPY requirements.txt /app/requirements.txt
RUN pip3 install -r /app/requirements.txt

# Copy the application source code.
COPY *.py /app
COPY *.sql /app

# Copy the schema files needed for VAT creation.
COPY *.json /data/variant_annotation_table/schema/

WORKDIR /app
