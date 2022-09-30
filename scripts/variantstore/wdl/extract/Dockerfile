FROM gcr.io/google.com/cloudsdktool/cloud-sdk:404.0.0-alpine as build

RUN apk update && apk upgrade
RUN python3 -m ensurepip --upgrade

RUN apk add autoconf bash cmake g++ gcc make ninja python3-dev git openssl-dev

# Build Apache Arrow version 8.0.0 as version 9.0.0 does not compile under Alpine:
# https://issues.apache.org/jira/browse/ARROW-17329
RUN cd / && \
    git clone https://github.com/apache/arrow.git && \
    cd arrow && \
    git checkout tags/apache-arrow-8.0.0 -b arrow-8.0.0

# Pyarrow build instructions from https://arrow.apache.org/docs/developers/python.html#python-development
# Modified slightly for the purposes of building this image.

# Intentionally do not use venv, the changes here should affect the "global" installation in the container.
RUN pip3 install -r /arrow/python/requirements-build.txt

RUN mkdir /dist
ARG ARROW_HOME=/dist
ARG LD_LIBRARY_PATH=/dist/lib:$LD_LIBRARY_PATH
RUN mkdir /arrow/cpp/build && \
    cd /arrow/cpp/build && \
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
RUN cd /arrow/python && \
    python3 setup.py build_ext --inplace

RUN pip3 install wheel
RUN python3 setup.py build_ext --build-type=release \
              --bundle-arrow-cpp bdist_wheel


FROM gcr.io/google.com/cloudsdktool/cloud-sdk:404.0.0-alpine as main

RUN apk update && apk upgrade
RUN python3 -m ensurepip --upgrade

RUN apk add openssl jq
COPY --from=build /arrow/python/dist/pyarrow-8.0.0-cp39-cp39-linux_x86_64.whl .
RUN pip3 install pyarrow-8.0.0-cp39-cp39-linux_x86_64.whl

COPY requirements.txt /app/requirements.txt
RUN pip3 install -r /app/requirements.txt

# Copy the application source code.
COPY *.py /app
COPY *.sql /app

# Copy the schema files needed for VAT creation.
COPY *.json /data/variant_annotation_table/schema/

WORKDIR /app
