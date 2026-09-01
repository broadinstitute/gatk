# The image produced by this Dockerfile contains tools and libraries required to support the Genomic Variant Store
# pipeline. The Alpine version of the Google Cloud SDK is used as the base image which is not only the most compact of
# the Google Cloud SDK Docker images.
#
# This is a multi-stage build that uses a custom "build-base" image for the build stage. The build-base image is
# expensive to create and isn't expected to change often, while the steps in this Dockerfile are much less expensive and
# more likely to change. Using a build-base image essentially allows the expensive layers to be globally cached which
# should make building the final image much faster in most cases.
#
# IMPORTANT: this tag is built and pushed separately by build_build_base_docker.sh from build_base.Dockerfile, and
# must be on the same Google Cloud SDK / Python version as the `main` stage's FROM below. If they diverge, the venv
# built here and copied into `main` via `COPY --from=build /localvenv /localvenv` gets a dangling `python3` symlink,
# so `python3` in `main` silently falls back to the bare system interpreter with none of the pip-installed packages.
# Bumping either FROM requires rebuilding build_base.Dockerfile, pushing a new tag, and updating both lines together.
FROM us.gcr.io/broad-dsde-methods/variantstore:2026-09-01-alpine-build-base AS build

# Install all of our variantstore Python requirements.
# We constrain the version of setuptools to sidestep an issue in `terra-notebook-util`'s FISS dependency.
# See https://github.com/broadinstitute/fiss/issues/192.
# requirements-prod.txt pulls in requirements.txt via `-r`, so both must be present.
COPY requirements.txt requirements.txt
COPY requirements-prod.txt requirements-prod.txt
RUN . /localvenv/bin/activate && \
    pip3 install 'setuptools<80.3.0' && \
    pip3 install -r requirements-prod.txt

# Install FISS from source and patch for a `ConfigParser` issue.
# See https://github.com/broadinstitute/fiss/issues/194.
ARG FISS_VERSION=0.16.38
COPY fiss_configparser.patch /fiss_configparser.patch
RUN mkdir /fiss-build && \
    cd /fiss-build && \
    curl -L -O https://github.com/broadinstitute/fiss/archive/refs/tags/v$FISS_VERSION.zip && \
    unzip v$FISS_VERSION.zip && \
    cd fiss-${FISS_VERSION} && \
    git apply /fiss_configparser.patch && \
    . /localvenv/bin/activate && \
    python setup.py install

# The main layer does not install development tools, instead copies artifacts from the build layer above.
# IMPORTANT: must stay on the same Python version as the `build` stage's FROM above -- see the note there.
FROM gcr.io/google.com/cloudsdktool/cloud-sdk:582.0.0-alpine AS main

RUN apk update && apk upgrade

# Add any generally useful Alpine packages here.
# `python3` (Alpine's own package, not this base image's bundled gcloud-internal interpreter -- see the note on the
# `build` stage FROM above) must match what /localvenv was built against in the `build` stage, since the venv copied
# in below is not self-contained and needs a matching /usr/bin/python3 and its shared libs present at the same path.
RUN apk add --no-cache perl jq python3

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
COPY /variant_annotation_table/left_alignment_fixups/*.py /app/
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
