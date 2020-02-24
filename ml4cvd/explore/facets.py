"""Methods for Facets notebook integration."""
import base64
import os

from facets_overview.generic_feature_statistics_generator import GenericFeatureStatisticsGenerator

FACETS_ABSOLUTE_URL = 'https://raw.githubusercontent.com/PAIR-code/facets/1.0.0/facets-dist/facets-jupyter.html'
FACETS_RELATIVE_URL = 'facets-jupyter.html'

FACETS_HTML = FACETS_ABSOLUTE_URL
if 'GOOGLE_PROJECT' in os.environ:  # This is Terra.
  # Terra notebook Content Security Policy prohibits importing this HTML from
  # a remote location, so this code depends on the fact we can refer to it
  # from a location relative to the notebook.
  if not os.path.exists(FACETS_RELATIVE_URL):
    # If this is our Terra Docker image, the file is available locally.
    docker_image_path = os.path.join(os.environ['JUPYTER_HOME'],
                                     'custom',
                                     FACETS_RELATIVE_URL)
    if os.path.exists(docker_image_path):
      os.system('ln -s ' + docker_image_path)
    else:
      os.system('wget --no-clobber ' + FACETS_ABSOLUTE_URL)
  FACETS_HTML = FACETS_RELATIVE_URL


class FacetsOverview(object):
  """Methods for Facets Overview notebook integration."""

  def __init__(self, data):
    # This takes the dataframe and computes all the inputs to the Facets
    # Overview plots such as:
    # - numeric variables: histogram bins, mean, min, median, max, etc..
    # - categorical variables: num unique, counts per category for bar chart,
    #     top category, etc.
    gfsg = GenericFeatureStatisticsGenerator()
    self._proto = gfsg.ProtoFromDataFrames(
        [{'name': 'data', 'table': data}])

  def _repr_html_(self):
    """Html representation of Facets Overview for use in a Jupyter notebook."""
    protostr = base64.b64encode(self._proto.SerializeToString()).decode('utf-8')
    html_template = '''
        <link rel="import" href="{facets_html}">
        <facets-overview id="overview_elem"></facets-overview>
        <script>
          document.querySelector("#overview_elem").protoInput = "{protostr}";
        </script>'''
    html = html_template.format(facets_html=FACETS_HTML, protostr=protostr)
    return html


class FacetsDive(object):
  """Methods for Facets Dive notebook integration."""

  def __init__(self, data, height=1000):
    self._data = data
    self.height = height

  def _repr_html_(self):
    """Html representation of Facets Dive for use in a Jupyter notebook."""
    html_template = """
        <link rel="import" href="{facets_html}">
        <facets-dive id="dive_elem" height="{height}"></facets-dive>
        <script>
          document.querySelector("#dive_elem").data = {data};
        </script>"""
    html = html_template.format(facets_html=FACETS_HTML,
                                data=self._data.to_json(orient='records'),
                                height=self.height)
    return html
