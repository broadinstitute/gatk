"""Methods for batch annotations of images stored as 3D tensors, such as MRIs, from within notebooks."""

import json
import os
import socket
import tempfile
from typing import Any, Dict, List, Tuple, Type, Union

from IPython.display import display
import numpy as np
import pandas as pd
import h5py
from ipyannotations import PolygonAnnotator
from ipyannotations.images.annotator import Annotator
import ipywidgets as widgets
from ml4h.visualization_tools.hd5_mri_plots import MRI_TMAPS
from ml4h.visualization_tools.annotation_storage import AnnotationStorage
from ml4h.visualization_tools.annotation_storage import TransientAnnotationStorage
from PIL import Image
import tensorflow as tf


class BatchImageAnnotator():
  """Annotate batches of images with shapes drawn over regions of interest."""

  SUBMIT_BUTTON_DESCRIPTION = 'Submit annotations, goto next sample'
  EXPECTED_COLUMN_NAMES = ['sample_id', 'tmap_name', 'instance_number', 'folder']
  DEFAULT_ANNOTATION_CLASSNAME = 'region_of_interest'
  CSS = '''
      <style>
      table, th, td {
        border: 2px solid black;
        border-collapse: collapse;
      }
      th, td {
        padding-left: 5px;
      }
      </style>
      '''

  def __init__(
      self, samples: pd.DataFrame, annotation_categories: List[str] = None,
      zoom: float = 1.5, annotation_storage: AnnotationStorage = TransientAnnotationStorage(),
      annotator: Type[Annotator] = PolygonAnnotator,
  ):
    """Initializes an instance of BatchImageAnnotator.

    Args:
      samples: A dataframe of samples to annotate. Columns must include those
               in BatchImageAnnotator.EXPECTED_COLUMN_NAMES.
      annotation_categories: A list of one or more strings to serve as tags for the annotations.
      zoom: Desired zoom level for the image. Defaults to 1.5.
      annotation_storage: An instance of AnnotationStorage. This faciltates the use of a user-provided
          strategy for the storage and processing of annotations. Defaults to  TransientAnnotationStorage.
      annotator: An instance of Annotator from package ipyannotations. Concrete classes include
        ipyannotations.BoxAnnotator, ipyannotations.PointAnnotator, and ipyannotations.PolygonAnnotator.
        Defaults to  ipyannotations.PolygonAnnotator.

    Raises:
      ValueError: The provided dataframe does not contain the expected columns.
    """
    if not set(self.EXPECTED_COLUMN_NAMES).issubset(samples.columns):
      raise ValueError(f'samples Dataframe must contain columns {self.EXPECTED_COLUMN_NAMES}')
    self.samples = samples
    self.current_sample = 0
    # TODO(deflaux) remove this after https://github.com/janfreyberg/ipyannotations/issues/11
    self.zoom = zoom
    self.annotation_storage = annotation_storage
    if annotation_categories is None:
      annotation_categories = [self.DEFAULT_ANNOTATION_CLASSNAME]

    self.annotation_widget = annotator(
        options=annotation_categories,
        canvas_size=(900, 280 * self.zoom),
    )
    self.annotation_widget.on_submit(self._store_annotations)
    self.annotation_widget.submit_button.description = self.SUBMIT_BUTTON_DESCRIPTION
    self.annotation_widget.submit_button.layout = widgets.Layout(width='300px')

    # Restructure the use instructions from the pydoc into a form that displays well as HTML.
    self.use_instructions = (
        '<ul><li>' +
        annotator.__doc__[0:annotator.__doc__.find('Parameters')].strip().replace('\n\n', '</li><li>') +
        '</li></ul>'
    )

    self.title_widget = widgets.HTML('')
    self.results_widget = widgets.HTML('')

  def _store_annotations(self, data: Dict[Any, Any]) -> None:
    """Transfer widget state to the annotation storage and advance to the next sample."""
    if self.current_sample >= self.samples.shape[0]:
      self.results_widget.value = '<h1>Annotation batch complete!</h1>Thank you for making the model better.'
      return

    # Convert canvas coordinates to tensor coordinates.
    image_canvas_position = self.annotation_widget.canvas.image_extent
    x_offset, y_offset, _, _ = image_canvas_position
    annotations = []
    for item in data:
      annotation: Dict[str, Union[List[Tuple[int, int]], Tuple[int, int], Tuple[int, int, int, int], str]] = {}
      annotations.append(annotation)
      for key in item.keys():
        if key == 'points':  # Polygons from PolygonAnnotator
          annotation[key] = [(
              int((p[0] - x_offset) / self.zoom),
              int((p[1] - y_offset) / self.zoom),
          ) for p in item[key]]
        elif key == 'coordinates':  # Points from PointAnnotator
          annotation[key] = (
              int((item[key][0] - x_offset) / self.zoom),
              int((item[key][1] - y_offset) / self.zoom),
          )
        elif key == 'xyxy':  # Rectangles from BoxAnnotator
          annotation[key] = (
              int((item[key][0] - x_offset) / self.zoom),
              int((item[key][1] - y_offset) / self.zoom),
              int((item[key][2] - x_offset) / self.zoom),
              int((item[key][3] - y_offset) / self.zoom),
          )
        else:
          # Pass all other values through unchanged.
          annotation[key] = item[key]

    # Store the annotation using the provided annotation storage strategy.
    self.annotation_storage.submit_annotation(
        sample_id=self.samples.loc[self.current_sample, 'sample_id'],
        annotator=os.getenv('OWNER_EMAIL') if os.getenv('OWNER_EMAIL') else socket.gethostname(),
        key=self.samples.loc[self.current_sample, 'tmap_name'],
        value_numeric=self.samples.loc[self.current_sample, 'instance_number'],
        value_string=self.samples.loc[self.current_sample, 'folder'],
        comment=json.dumps(annotations),
    )

    # Display this annotation at the bottom of the widget.
    results = f'''
        <hr>
        <h2>Prior sample's submitted annotations</h2>
        The <b>{self.SUBMIT_BUTTON_DESCRIPTION}</b> button is both printing out the annotations below and storing the annotations
        via strategy {self.annotation_storage.__class__.__name__}.<br>
        Details: <i>{self.annotation_storage.describe()}</i>
        <h3>sample info</h3>
        {self._format_info_for_current_sample()}
        <h3>canvas coordinates</h3>
        image extent {image_canvas_position}
        {[f'<pre>{json.dumps(x)}</pre>' for x in data]}
        <h3>source tensor coordinates</h3>
        {[f'<pre>{json.dumps(x)}</pre>' for x in annotations]}
        <hr>
      '''
    self.results_widget.value = results

    # Advance to the next sample.
    self.current_sample += 1
    self._annotate_image_for_current_sample()

  def _format_info_for_current_sample(self) -> str:
    """Convert information about the current sample to an HTML table for display within the widget."""
    headings = ' '.join([f'<th>{c}</th>' for c in self.EXPECTED_COLUMN_NAMES] + ['<th>TMAP shape</th>'])
    values = ' '.join([f'<td>{self.samples.loc[self.current_sample, c]}</td>' for c in self.EXPECTED_COLUMN_NAMES]
                      + [f'<td>{MRI_TMAPS[self.samples.loc[self.current_sample, "tmap_name"]].shape}</td>'])
    return f'''
        <table style="width:100%">
        <tr>{headings}</tr>
        <tr>{values}</tr>
        </table>
        '''

  def _annotate_image_for_current_sample(self) -> None:
    """Retrieve the data for the current sample and display its image in the annotation widget.

    If all samples have been processed, display the completion message.
    """
    if self.current_sample >= self.samples.shape[0]:
      self.annotation_widget.canvas.clear()
      self.title_widget.value = '<h1>Annotation batch complete!</h1>Thank you for making the model better.'
      return

    sample_id = self.samples.loc[self.current_sample, 'sample_id']
    tmap_name = self.samples.loc[self.current_sample, 'tmap_name']
    instance_number = self.samples.loc[self.current_sample, 'instance_number']
    folder = self.samples.loc[self.current_sample, 'folder']

    with tempfile.TemporaryDirectory() as tmpdirname:
      sample_hd5 = str(sample_id) + '.hd5'
      local_path = os.path.join(tmpdirname, sample_hd5)
      try:
        tf.io.gfile.copy(src=os.path.join(folder, sample_hd5), dst=local_path)
        hd5 = h5py.File(local_path, mode='r')
      except (tf.errors.NotFoundError, tf.errors.PermissionDeniedError) as e:
        self.annotation_widget.canvas.clear()
        self.title_widget.value = f'''
            <div class="alert alert-block alert-danger">
            <H2>Warning: MRI HD5 file not available for sample {sample_id} in folder {folder}</h2>
            Use the <kbd>folder</kbd> parameter to read HD5s from a different local directory or Cloud Storage bucket.
            <hr><p><pre>{e.message}</pre></p>
            </div>
            '''
        return

    tensor = MRI_TMAPS[tmap_name].tensor_from_file(MRI_TMAPS[tmap_name], hd5)
    tensor_instance = tensor[:, :, instance_number]
    if self.zoom > 1.0:
      # TODO(deflaux) remove this after https://github.com/janfreyberg/ipyannotations/issues/11
      img = Image.fromarray(tensor_instance)
      zoomed_img = img.resize([int(self.zoom * s) for s in img.size], Image.LANCZOS)
      tensor_instance = np.asarray(zoomed_img)

    self.annotation_widget.display(tensor_instance)
    self.title_widget.value = f'''
        {self.CSS}
        <div class='alert alert-block alert-info'>
        <h1>Batch annotation of {self.samples.shape[0]} samples</h1>
        {self.use_instructions}
        </div>
        <h2>Current sample</h2>
        {self._format_info_for_current_sample()}
        <hr>
    '''

  def annotate_images(self) -> None:
    """Begin the batch annotation task by displaying the annotation widget populated with the first sample.

    The submit button is used to proceed to the next sample until all samples have been processed.
    """
    self._annotate_image_for_current_sample()
    display(widgets.VBox([self.title_widget, self.annotation_widget, self.results_widget]))

  def view_recent_submissions(self, count: int = 10) -> pd.DataFrame:
    """View a dataframe of up to [count] most recent submissions.

    Args:
      count: The number of the most recent submissions to return.

    Returns:
      A dataframe of the most recent annotations.
    """
    return self.annotation_storage.view_recent_submissions(count=count)
