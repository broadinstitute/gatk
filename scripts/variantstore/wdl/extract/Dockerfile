# The image produced by this Dockerfile contains tools and libraries required to support the Genomic Variant Store
# pipeline. The Alpine version of the Google Cloud SDK is used as the base image which is not only the most compact of
# the Google Cloud SDK Docker images, but is also the image currently used by Cromwell for (de)localization of files in
# Google Cloud Storage. Sharing the base image with Cromwell's GCS localization should result in reuse of a cached copy
# of this base (and by far largest) image layer when running GVS pipelines in Terra / Cromwell.
#
# This is a multi-stage build that uses a custom "build-base" image for the build stage. The build-base image is
# expensive to create and isn't expected to change often, while the steps in this Dockerfile are much less expensive and
# more likely to change. Using a build-base image essentially allows the expensive layers to be globally cached which
# should make building the final image much faster in most cases.
FROM us.gcr.io/broad-dsde-methods/variantstore:2023-04-13-alpine-build-base as build

# Install all of our variantstore Python requirements.
COPY requirements.txt requirements.txt
RUN pip3 install --user -r requirements.txt

# The main layer does not install development tools, instead copies artifacts from the build layer above.
FROM gcr.io/google.com/cloudsdktool/cloud-sdk:426.0.0-alpine as main

RUN apk update && apk upgrade
RUN python3 -m ensurepip --upgrade

# Add any generally useful Alpine packages here.
RUN apk add --no-cache jq

# The build stage generated Python artifacts to /root/.local via `pip install --user`, so grab all of those.
COPY --from=build /root/.local /root/.local
# Copy bcftools artifacts
COPY --from=build /bcftools /bcftools

# Copy the application source code.
RUN mkdir /app
COPY *.py /app/
COPY *.sql /app/

# Copy the schema files needed for VAT creation.
RUN mkdir -p /data/variant_annotation_table/schema
COPY ./variant_annotation_table/schema/*.json /data/variant_annotation_table/schema/

ENV PATH=/bcftools/bin\:$PATH
WORKDIR /app
