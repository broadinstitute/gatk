FROM gcr.io/google.com/cloudsdktool/cloud-sdk:402.0.0-alpine as build

# Edit /etc/apk/repositories to point to edge (required for py3-apache-arrow):
# https://wiki.alpinelinux.org/wiki/Repositories
RUN sed -i 's;v3.15;edge;' /etc/apk/repositories
RUN apk update
# Use prebuilt py3-apache-arrow package for convenience
# https://stackoverflow.com/a/73574662
RUN apk add --no-cache jq gcc musl-dev libffi-dev python3-dev cmake py3-numpy py3-apache-arrow

RUN curl -O https://bootstrap.pypa.io/get-pip.py
RUN python3 get-pip.py

# Copy the application's requirements.txt and run pip to install.
COPY requirements.txt /app/requirements.txt
RUN pip install --user -r /app/requirements.txt

# Build bcftools as there's no Alpine package.
RUN apk add --no-cache git make xz-dev zlib-dev bzip2-dev curl-dev
RUN git clone --recurse-submodules https://github.com/samtools/htslib.git
RUN git clone https://github.com/samtools/bcftools.git
RUN cd bcftools && make && make install

FROM gcr.io/google.com/cloudsdktool/cloud-sdk:402.0.0-alpine as main
WORKDIR /app

RUN sed -i 's;v3.15;edge;' /etc/apk/repositories
RUN apk update

# Copy Python modules from build stage to main stage.
COPY --from=build /root/.local /root/.local
# Upgrade to Python 3.10, required for compatibility with the build stage above.
RUN apk upgrade python3
RUN apk add --no-cache jq py3-numpy py3-apache-arrow xz zlib bzip2 curl

# Copy the application source code.
COPY *.py /app
COPY *.sql /app

# Copy the schema files needed for VAT creation.
COPY variant_annotation_table/schema/*.json /data/variant_annotation_table/schema/

ENV PATH=/root/.local/bin:$PATH
ENV PYTHONPATH=/root/.local/lib/python3.10:$PYTHONPATH

# Copy bcftools artifacts from build stage to main stage.
COPY --from=build /usr/local/bin /usr/local/bin
RUN mkdir -p /usr/local/libexec/bcftools
COPY --from=build /usr/local/libexec/bcftools /usr/local/libexec/bcftools
