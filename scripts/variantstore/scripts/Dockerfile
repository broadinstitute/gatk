# The image produced by this Dockerfile contains tools and libraries required to support the Genomic Variant Store
# pipeline. The Alpine version of the Google Cloud SDK is used as the base image which is not only the most compact of
# the Google Cloud SDK Docker images.
#
# This is a multi-stage build that uses a custom "build-base" image for the build stage. The build-base image is
# expensive to create and isn't expected to change often, while the steps in this Dockerfile are much less expensive and
# more likely to change. Using a build-base image essentially allows the expensive layers to be globally cached which
# should make building the final image much faster in most cases.
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

# Patch `terra-notebook-utils` for a urllib3 2.x incompatibility.
# Its `table.py` still passes urllib3 1.x's `method_whitelist` keyword to `Retry`; urllib3 2.0 removed that
# keyword in favor of `allowed_methods`. Upstream migrated its own `http.py` but missed `table.py`, so every
# Terra entity-table call raises
# `TypeError: Retry.__init__() got an unexpected keyword argument 'method_whitelist'`.
# This was latent until the Cloud SDK base image moved to Python 3.14: on 3.12 pip happened to resolve
# urllib3 1.26.x (where the keyword still exists but is deprecated), so nothing failed. That was luck, not a
# constraint anyone wrote down, and it stopped holding on 3.14, where pip resolves urllib3 2.x.
# Patching is preferred over pinning `urllib3 < 2`, which would hold the whole image back to an end-of-life
# urllib3. `git apply` fails the build if the upstream line ever changes, which is the signal to re-check
# whether this workaround can be dropped.
COPY terra_notebook_utils_urllib3.patch /terra_notebook_utils_urllib3.patch
RUN . /localvenv/bin/activate && \
    cd "$(python3 -c 'import terra_notebook_utils as m, pathlib; print(pathlib.Path(m.__file__).parent.parent)')" && \
    git apply /terra_notebook_utils_urllib3.patch && \
    python3 -c "from urllib3.util.retry import Retry; Retry(total=10, status_forcelist=[429], allowed_methods=['GET'])"

# The main layer does not install development tools, instead copies artifacts from the build layer above.
# IMPORTANT: must stay on the same Python version as the `build` stage's FROM above -- see the note there.
# 565.0.0 is a deliberate ceiling -- see the note above this Dockerfile's `build` stage FROM, in
# build_base.Dockerfile, and in GvsUtils.wdl's `cloud_sdk_docker_decl`. 568.0.0 and later bundle Python 3.14,
# which hail cannot be installed on.
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
