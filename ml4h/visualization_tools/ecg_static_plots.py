"""Methods for integration of static plots within notebooks."""
import os
import tempfile

from IPython.display import HTML
from IPython.display import SVG
from ml4h.plots import plot_ecg_rest
from ml4h.runtime_data_defines import get_resting_ecg_hd5_folder
from ml4h.runtime_data_defines import get_resting_ecg_svg_folder
import numpy as np
import tensorflow as tf


def display_resting_ecg(sample_id, folder=None):
  """Retrieve (or render) and display the SVG of the resting ECG.

  Args:
    sample_id: The id of the ECG SVG to retrieve.
    folder: The local or Cloud Storage path under which the files reside.

  Returns:
    An IPython SVG object or a notebook-friendly error.
  """
  if folder is None:
    svg_folder = get_resting_ecg_svg_folder(sample_id)
    hd5_folder = get_resting_ecg_hd5_folder(sample_id)
  else:
    svg_folder = folder
    hd5_folder = folder

  with tempfile.TemporaryDirectory() as tmpdirname:
    # First, see if we already have one rendered.
    sample_svg = str(sample_id) + '.svg'
    local_path = os.path.join(tmpdirname, sample_svg)
    try:
      tf.io.gfile.copy(src=os.path.join(svg_folder, sample_svg), dst=local_path)
      return SVG(filename=local_path)
    except (tf.errors.NotFoundError, tf.errors.PermissionDeniedError) as e:
      pass
    # If not, dynamically render a SVG
    sample_hd5 = str(sample_id) + '.hd5'
    local_path = os.path.join(tmpdirname, sample_hd5)
    try:
      tf.io.gfile.copy(src=os.path.join(hd5_folder, sample_hd5), dst=local_path)
    except (tf.errors.NotFoundError, tf.errors.PermissionDeniedError) as e:
      return HTML(f'''
      <div class="alert alert-block alert-danger">
      <b>Warning:</b> Resting ECG not available for sample {sample_id} in {svg_folder} or {hd5_folder}:
      <hr><p><pre>{e.message}</pre></p>
      Use the <kbd>folder</kbd> parameter to read from a different local directory or Cloud Storage bucket.
      </div>''')

    try:
      # We don't need the resulting SVG, so send it to a temporary directory.
      with tempfile.TemporaryDirectory() as tmpdirname:
        plot_ecg_rest(tensor_paths = [local_path], rows=[0], out_folder=tmpdirname, is_blind=False)
    except Exception as e:
      return HTML(f'''
        <div class="alert alert-block alert-danger">
        <b>Warning:</b> Unable to render static plot of resting ECG for sample {sample_id} from {hd5_folder}:
        <hr><p><pre>{e}</pre></p>
        </div>''')


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
