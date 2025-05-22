# The image produced by this Dockerfile contains tools and libraries required to support the Genomic Variant Store
# pipeline. The Alpine version of the Google Cloud SDK is used as the base image which is not only the most compact of
# the Google Cloud SDK Docker images.
#
# This is a multi-stage build that uses a custom "build-base" image for the build stage. The build-base image is
# expensive to create and isn't expected to change often, while the steps in this Dockerfile are much less expensive and
# more likely to change. Using a build-base image essentially allows the expensive layers to be globally cached which
# should make building the final image much faster in most cases.
FROM us.gcr.io/broad-dsde-methods/variantstore:2025-05-21-alpine-build-base as build

# Install all of our variantstore Python requirements. We need that specific version of setuptools to avoid as it was
# the last to include the 'package_index' module required by our 'firecloud' (FISS) dependency.
COPY requirements.txt requirements.txt
RUN . /localvenv/bin/activate && \
    pip3 install setuptools==80.2.0 && \
    pip3 install -r requirements.txt

# The main layer does not install development tools, instead copies artifacts from the build layer above.
FROM gcr.io/google.com/cloudsdktool/cloud-sdk:522.0.0-alpine as main

RUN apk update && apk upgrade

# Add any generally useful Alpine packages here.
RUN apk add --no-cache perl jq

# The build stage generated Python artifacts to /root/.local via `pip install --user`, so grab all of those.
COPY --from=build /localvenv /localvenv
# Copy bcftools artifacts
COPY --from=build /bcftools /bcftools
# Copy vcftools artifacts
COPY --from=build /vcftools /vcftools

# Copy the application source code.
RUN mkdir /app
COPY *.py /app/
COPY *.sql /app/

# Copy the schema files needed for VAT creation.
RUN mkdir -p /data/variant_annotation_table/schema
COPY ./variant_annotation_table/schema/*.json /data/variant_annotation_table/schema/

ENV PERL5LIB=/vcftools/share/perl5/site_perl/:$PERL5LIB
ENV PATH=/bcftools/bin\:/vcftools/bin\:$PATH

WORKDIR /app
