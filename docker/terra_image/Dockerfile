FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-r:0.0.7
# https://github.com/DataBiosphere/terra-docker/blob/master/terra-jupyter-r/CHANGELOG.md

USER root

# Fix annoying pip warnings in the source image.
RUN chmod -R a+rwx $HOME/.cache

# Add to the PATH so that the tensorflow upgrade does not emit warnings.
ENV PATH $PATH:$HOME/.local/bin

# Temporarily make pip install to the system directory.
ENV PIP_USER=false

# Add minimal ml4cvd dependencies.
RUN apt-get install libxt-dev -y \
  # TODO(deflaux) update and use https://github.com/broadinstitute/ml/blob/master/docker/vm_boot_images/config/tensorflow-requirements.txt
  # instead, after the Tensorflow 2.0 pull request is merged.
  && pip3 install \
  apache-beam[gcp]==2.12.0 \
  google-cloud-storage==1.25.0 \
  h5py==2.9.0 \
  imageio==2.6.1 \
  keras==2.2.5 \
  nibabel==2.5.0 \
  seaborn==0.9.0 \
  vtk==8.1.2

# Add visualization libraries not installed by default on Terra and upgrade
# to newer versions for other libraries that are installed by default.
RUN pip3 install --upgrade \

  # Install signal processing libraries.
  biosppy \

  # Upgrade to latest version of this so that the tensorflow upgrade does not
  # emit warnings.
  setuptools \

  # Upgrade to latest version of this so that %%bigquery does not print the
  # first few rows of the downloaded dataframe of query results.
  # Pin version due to https://github.com/googleapis/google-cloud-python/issues/9965
  google-cloud-bigquery[pandas]==1.22.0 \

  # Upgrade to latest version of tensorflow.
  tensorflow \

  # Upgrade to latest version of statsmodels to resolve an issue with plotnine.
  statsmodels \

  # Install data visualization libraries.
  facets-overview \
  # Pin version due to Terra's more strict 'Content Security Policy'.
  altair==3.3.0 \
  plotnine \
  pydicom \
  vega \
  vega_datasets \

  # Configure notebook extensions.
  && jupyter nbextension install --py vega \
  && jupyter nbextension enable --py vega

# TODO(deflaux): this also needs libgl1-mesa-glx and/or python3-tk to fix error
# "ModuleNotFoundError: No module named 'vtkOpenGLKitPython'"
# but neither of those have installation candidate for this base image.

# Terra notebook Content Security Policy prohibits importing this HTML from a
# remote location, so we download a local copy instead. The instructions in
# https://github.com/PAIR-code/facets#setup do not apply to Terra due to the way
# URLs are constructed. In the python library code that refers to this HTML
# file, we create a symlink to it so that the notebook-relative URL
# ./facets-jupyter.html will load.
RUN cd $JUPYTER_HOME/custom \
  && wget https://raw.githubusercontent.com/PAIR-code/facets/1.0.0/facets-dist/facets-jupyter.html

# Make pip install to a user directory, instead of a system directory which
# requires root. This is useful so `pip install` commands can be run in the
# context of a notebook.
ENV PIP_USER=true

# Install custom private python package from Puneet Batra's team. Enable easy
# local edits by placing it in the home directory with write access for
# jupyter-user.
USER $USER
RUN mkdir -p $HOME/ml4cvd_pkg
COPY --chown=jupyter-user:users ml4cvd $HOME/ml4cvd_pkg/ml4cvd
ENV PYTHONPATH $PYTHONPATH:$HOME/ml4cvd_pkg

