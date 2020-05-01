# The suggested base images are:
#   - ufoym/deepo:all-py36-jupyter for GPU-enabled machines
#   - ufoym/deepo:all-py36-jupyter-cpu for CPU-only (non-GPU-enabled) machines
# BASE_IMAGE can be specified at build time by adding the following argument:
#   --build_arg BASE_IMAGE="some_other_image"

ARG BASE_IMAGE
FROM ${BASE_IMAGE}

LABEL maintainer="James Pirruccello <jamesp@broadinstitute.org>"

# Setup time zone (or else docker build hangs)
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

COPY ./config/* /app/
WORKDIR /app

# Note that some layers are kept separate to encourage layer re-use and to try
# to minimize full recompilation where possible.

# Basic setup
RUN ./ubuntu.sh

# FastAI. See the Developer Install under https://github.com/fastai/fastai/ to
# understand this odd sequence of installing then uninstalling fastai before 
# installing it from github. (Basically, to get its deps.)
# RUN pip3 install -r fastai-requirements.txt
# RUN pip3 uninstall -y fastai
# RUN ./fastai.sh

RUN apt-get install python3-tk libgl1-mesa-glx libxt-dev -y

# Requirements for the tensorflow project
RUN pip3 install --upgrade pip
RUN pip3 install -r pre_requirements.txt
RUN pip3 install -r tensorflow-requirements.txt
