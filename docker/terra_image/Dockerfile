FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-gatk:1.0.0
# https://github.com/DataBiosphere/terra-docker/blob/master/terra-jupyter-gatk/CHANGELOG.md

USER root
# Address errors due to enum34 package https://stackoverflow.com/a/45716067/4138705
RUN pip3 uninstall -y enum34

# Install custom private python package from Puneet Batra's team. Enable easy
# local edits by placing it in the home directory with write access for
# jupyter-user.
USER $USER
RUN mkdir -p $HOME/ml4cvd_pkg
COPY --chown=jupyter-user:users ml4cvd $HOME/ml4cvd_pkg/ml4cvd
COPY --chown=jupyter-user:users config $HOME/ml4cvd_pkg/config
ENV PYTHONPATH $PYTHONPATH:$HOME/ml4cvd_pkg

RUN pip3 install --user -r $HOME/ml4cvd_pkg/config/tensorflow-requirements.txt \
  # Upgrade to a newer version so that %%bigquery does not print the
  # first few rows of the downloaded dataframe of query results.
  # Pin version due to https://github.com/googleapis/google-cloud-python/issues/9965
  && pip3 install --upgrade --user google-cloud-bigquery[pandas]==1.22.0 \
  # Configure notebook extensions.
  && jupyter nbextension install --user --py vega \
  && jupyter nbextension enable --user --py vega
