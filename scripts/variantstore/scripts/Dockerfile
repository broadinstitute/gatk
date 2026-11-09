# The image produced by this Dockerfile contains tools and libraries required to support the Genomic Variant Store
# pipeline. The Alpine version of the Google Cloud SDK is used as the base image which is not only the most compact of
# the Google Cloud SDK Docker images.
#
# This is a multi-stage build that uses a custom "build-base" image for the build stage. The build-base image is
# expensive to create and isn't expected to change often, while the steps in this Dockerfile are much less expensive and
# more likely to change. Using a build-base image essentially allows the expensive layers to be globally cached which
# should make building the final image much faster in most cases.
FROM us.gcr.io/broad-dsde-methods/variantstore:2026-09-02-alpine-build-base AS build

# Install all of our variantstore Python requirements.
# requirements-prod.txt pulls in requirements.txt via `-r`, so both must be present.
COPY requirements.txt requirements.txt
COPY requirements-prod.txt requirements-prod.txt
RUN . /localvenv/bin/activate && \
    pip3 install -r requirements-prod.txt

# Patch `terra-notebook-utils` for a urllib3 2.x incompatibility.
# A PR has been opened in `terra-notebook-utils` to fix this upstream, but until that lands we patch this here.
# https://github.com/DataBiosphere/terra-notebook-utils/pull/426
COPY terra_notebook_utils_urllib3.patch /terra_notebook_utils_urllib3.patch
RUN . /localvenv/bin/activate && \
    cd "$(python3 -c 'import terra_notebook_utils as m, pathlib; print(pathlib.Path(m.__file__).parent.parent)')" && \
    git apply /terra_notebook_utils_urllib3.patch

# The main layer does not install development tools, instead copies artifacts from the build layer above.
# IMPORTANT: must stay on the same Python version as the `build` stage's FROM above -- see the note there.
FROM gcr.io/google.com/cloudsdktool/cloud-sdk:565.0.0-alpine AS main

RUN apk update && apk upgrade

# Add any generally useful Alpine packages here.
RUN apk add --no-cache perl jq

# The build stage generated Python artifacts to /root/.local via `pip install --user`, so grab all of those.
COPY --from=build /localvenv /localvenv
# Copy htslib artifacts
COPY --from=build /htslib /htslib
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

ENV PERL5LIB=/vcftools/share/perl5/site_perl/
ENV PATH=/bcftools/bin:/vcftools/bin:$PATH
ENV LD_LIBRARY_PATH=/htslib/lib

# Have the bash shell of the command script automatically activate the Python virtual environment.
COPY activate_localvenv.sh /root/.bashrc
ENV BASH_ENV=/root/.bashrc

WORKDIR /app
