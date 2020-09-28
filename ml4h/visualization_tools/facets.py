"""Methods for Facets notebook integration."""
import base64
import os

import pandas as pd
from facets_overview.generic_feature_statistics_generator import GenericFeatureStatisticsGenerator

FACETS_DEPENDENCIES = {
    'facets_html': 'https://raw.githubusercontent.com/PAIR-code/facets/1.0.0/facets-dist/facets-jupyter.html',
    'webcomponents_js': 'https://cdnjs.cloudflare.com/ajax/libs/webcomponentsjs/1.3.3/webcomponents-lite.js',
}

if 'GOOGLE_PROJECT' in os.environ:  # This is Terra.
  # Terra notebook Content Security Policy prohibits pulling these files from
  # a remote location, so this code depends on the fact we can refer to it
  # from a location relative to the notebook.
  for dep, url in FACETS_DEPENDENCIES.items():
    if not os.path.exists(os.path.basename(url)):
      # If this is our Terra Docker image, the file is available locally.
      docker_image_path = os.path.join(os.environ['JUPYTER_HOME'], 'custom', os.path.basename(url))
      if os.path.exists(docker_image_path):
        os.system('ln -s ' + docker_image_path)
      else:
        os.system('wget --no-clobber ' + url)
    # Update dictionary to replace absolute url with relative url.
    FACETS_DEPENDENCIES[dep] = os.path.basename(url)


class FacetsOverview():
  """Methods for Facets Overview notebook integration."""

  def __init__(self, data: pd.DataFrame):
    # This takes the dataframe and computes all the inputs to the Facets
    # Overview plots such as:
    # - numeric variables: histogram bins, mean, min, median, max, etc..
    # - categorical variables: num unique, counts per category for bar chart,
    #     top category, etc.
    gfsg = GenericFeatureStatisticsGenerator()
    self._proto = gfsg.ProtoFromDataFrames(
        [{'name': 'data', 'table': data}],
    )

  def _repr_html_(self) -> str:
    """Html representation of Facets Overview for use in a Jupyter notebook."""
    protostr = base64.b64encode(self._proto.SerializeToString()).decode('utf-8')
    html_template = '''
        <script src="{webcomponents_js}"></script>
        <link rel="import" href="{facets_html}">
        <facets-overview id="overview_elem"></facets-overview>
        <script>
          document.querySelector("#overview_elem").protoInput = "{protostr}";
        </script>'''
    html = html_template.format(
        facets_html=FACETS_DEPENDENCIES['facets_html'],
        webcomponents_js=FACETS_DEPENDENCIES['webcomponents_js'],
        protostr=protostr,
    )
    return html


class FacetsDive():
  """Methods for Facets Dive notebook integration."""

  def __init__(self, data: pd.DataFrame, height: int = 1000):
    self._data = data
    self.height = height

  def _repr_html_(self) -> str:
    """Html representation of Facets Dive for use in a Jupyter notebook."""
    html_template = """
        <script src="{webcomponents_js}"></script>
        <link rel="import" href="{facets_html}">
        <facets-dive id="dive_elem" height="{height}"></facets-dive>
        <script>
          document.querySelector("#dive_elem").data = {data};
        </script>"""
    html = html_template.format(
        facets_html=FACETS_DEPENDENCIES['facets_html'],
        webcomponents_js=FACETS_DEPENDENCIES['webcomponents_js'],
        data=self._data.to_json(orient='records'),
        height=self.height,
    )
    return html
