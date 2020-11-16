FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-gatk:1.0.9
# https://github.com/DataBiosphere/terra-docker/blob/master/terra-jupyter-gatk/CHANGELOG.md

USER root
# TODO(deflaux) remove this after https://github.com/DataBiosphere/terra-docker/issues/175
RUN pip3 uninstall -y enum34

USER $USER
RUN pip3 install git+https://github.com/broadinstitute/ml4h.git \
  # Configure notebook extensions.
  && jupyter nbextension install --user --py vega \
  && jupyter nbextension enable --user --py vega \
  && jupyter nbextension install --user --py ipycanvas \
  && jupyter nbextension enable --user --py ipycanvas
