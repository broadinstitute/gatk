FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-r:0.0.7
# https://github.com/DataBiosphere/terra-docker/blob/master/terra-jupyter-r/CHANGELOG.md

USER root

# Fix annoying pip warnings in the source image.
RUN chmod -R a+rwx $HOME/.cache

# Add to the PATH so that the tensorflow upgrade does not emit warnings.
ENV PATH $PATH:$HOME/.local/bin

# Temporarily make pip install to the system directory.
ENV PIP_USER=false

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

USER $USER

# Install custom private python package from Puneet Batra's team. Due to the
# preceeding USER command, this is done as the notebook user so that the file
# permissions allow for editing of the code.
RUN mkdir -p $HOME/ml4cvd/ecg_demo
COPY ecg_demo $HOME/ml4cvd/ecg_demo
ENV PYTHONPATH $PYTHONPATH:$HOME/ml4cvd

