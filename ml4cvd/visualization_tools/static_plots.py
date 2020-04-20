"""Methods for integration of static plots within notebooks."""
import os
import tempfile

from IPython.display import HTML
from IPython.display import SVG
from ml4cvd.runtime_data_defines import get_resting_ecg_svg_folder
import numpy as np
import tensorflow as tf


def display_resting_ecg(sample_id, folder=None):
  """Retrieve and display the SVG of the resting ECG.

  Args:
    sample_id: The id of the ECG SVG to retrieve.
    folder: The local or Cloud Storage path under which the files reside.

  Returns:
    An IPython SVG object or a notebook-friendly error.
  """
  if folder is None:
    folder = get_resting_ecg_svg_folder(sample_id)

  with tempfile.TemporaryDirectory() as tmpdirname:
    sample_svg = str(sample_id) + '.svg'
    local_path = os.path.join(tmpdirname, sample_svg)
    try:
      tf.io.gfile.copy(src=os.path.join(folder, sample_svg), dst=local_path)
    except (tf.errors.NotFoundError, tf.errors.PermissionDeniedError) as e:
      return HTML(f'''
      <div class="alert alert-block alert-danger">
      <b>Warning:</b> Resting ECG image not available for sample {sample_id}:
      <hr><p><pre>{e.message}</pre></p>
      </div>''')

    return SVG(filename=local_path)


def major_breaks_x_resting_ecg(limits):
  """Method to compute breaks for plotnine plots of ECG resting data.

  Args:
    limits: The approximate limits.

  Returns:
    The desired limits.
  """
  step = 0.2
  if limits[0] <= 0:
    min_break = 0.0
    max_break = 2.5
  elif limits[0] <= 2.5:
    min_break = 2.5
    max_break = 5.0
  elif limits[0] <= 5.0:
    min_break = 5.0
    max_break = 7.5
  else:
    min_break = 7.5
    max_break = 10.0
  return np.arange(min_break, max_break + step, step)
