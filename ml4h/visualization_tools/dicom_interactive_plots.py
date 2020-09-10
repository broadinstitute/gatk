"""Methods for integration of interactive dicom plots within notebooks.

TODO:
* Continue to *pragmatically* improve this to make the visualization controls
  and the information shown more similar to that in full-featured DICOM viewers.
"""

import collections
import os
import tempfile
import zipfile

from IPython.display import display
from IPython.display import HTML
import ipywidgets as widgets
import matplotlib.pyplot as plt
from ml4h.runtime_data_defines import get_mri_folders
import numpy as np
import pydicom
import tensorflow as tf

MIN_IMAGE_WIDTH = 8
DEFAULT_IMAGE_WIDTH = 12
MAX_IMAGE_WIDTH = 24

MIN_COLOR_RANGE = 0
MAX_COLOR_RANGE = 6000


def choose_mri(sample_id, folder=None):
  """Render widget to choose the MRI to plot.

  Args:
    sample_id: The id of the sample to retrieve.
    folder: The local or Cloud Storage folder under which the files reside.

  Returns:
    ipywidget or HTML upon error.
  """
  if folder is None:
    folders = get_mri_folders(sample_id)
  else:
    folders = [folder]

  sample_mris = []
  sample_mri_glob = str(sample_id) + '_*.zip'
  try:
    for folder in folders:
      sample_mris.extend(tf.io.gfile.glob(pattern=os.path.join(folder, sample_mri_glob)))
  except (tf.errors.NotFoundError, tf.errors.PermissionDeniedError) as e:
    return HTML(f'''
    <div class="alert alert-block alert-danger">
    <b>Warning:</b> MRI not available for sample {sample_id} in {folders}:
    <hr><p><pre>{e.message}</pre></p>
    Use the <kbd>folder</kbd> parameter to read DICOMs from a different local directory or Cloud Storage bucket.
    </div>''')

  if not sample_mris:
    return HTML(f'''
    <div class="alert alert-block alert-danger">
    <b>Warning:</b> MRI DICOMs not available for sample {sample_id} in {folders}.<br>
    Use the <kbd>folder</kbd> parameter to read DICOMs from a different local directory or Cloud Storage bucket.
    </div>''')

  mri_chooser = widgets.Dropdown(
      options=sample_mris,
      value=sample_mris[0],
      description=f'Choose an MRI to visualize for sample {sample_id}:',
      style={'description_width': 'initial'},
      layout=widgets.Layout(width='800px'),
  )
  file_controls_ui = widgets.VBox(
      [widgets.HTML('<h3>File controls</h3>'), mri_chooser],
      layout=widgets.Layout(width='auto', border='solid 1px grey'),
  )
  file_controls_output = widgets.interactive_output(choose_mri_series, {'sample_mri': mri_chooser})
  display(file_controls_ui, file_controls_output)


def choose_mri_series(sample_mri):
  """Render widgets and interactive plots for MRIs.

  Args:
    sample_mri: The local or Cloud Storage path to the MRI file.

  Returns:
    ipywidget or HTML upon error.
  """
  with tempfile.TemporaryDirectory() as tmpdirname:
    local_path = os.path.join(tmpdirname, os.path.basename(sample_mri))
    try:
      tf.io.gfile.copy(src=sample_mri, dst=local_path)
      with zipfile.ZipFile(local_path, 'r') as zip_ref:
        zip_ref.extractall(tmpdirname)
    except (tf.errors.NotFoundError, tf.errors.PermissionDeniedError) as e:
      return HTML(f'''
      <div class="alert alert-block alert-danger">
      <b>Warning:</b> Cardiac MRI not available for sample {os.path.basename(sample_mri)}:
      <hr><p><pre>{e.message}</pre></p>
      </div>''')

    unordered_dicoms = collections.defaultdict(dict)
    for dcm_file in os.listdir(tmpdirname):
      if not dcm_file.endswith('.dcm'):
        continue
      dcm = pydicom.read_file(os.path.join(tmpdirname, dcm_file))
      key1 = (dcm.SeriesDescription.lower(), int(dcm.SeriesNumber))
      key2 = int(dcm.InstanceNumber) - 1
      if key2 in unordered_dicoms[key1]:
        # Notice invalid input, but don't throw an error.
        print(f'WARNING: Duplicate instances: {dcm.SeriesDescription} {dcm.SeriesNumber} {dcm.InstanceNumber}.')
      unordered_dicoms[key1][key2] = dcm

  if not unordered_dicoms:
    print(f'\n\nNo series available in MRI for sample {os.path.basename(sample_mri)}\n\nTry a different MRI.')
    return None

  # Convert from dict of dicts to dict of ordered lists.
  dicoms = {}
  for series in unordered_dicoms.keys():
    dicoms[series] = [None] * (max(unordered_dicoms[series]) + 1)
    for idx, val in unordered_dicoms[series].items():
      dicoms[series][idx] = val

  default_series_value = sorted(list(dicoms.keys()))[0]
  # Display the middle instance by default.
  default_instance_value, max_instance_value = compute_instance_range(dicoms, default_series_value)
  default_vmin_value, default_vmax_value = compute_color_range(dicoms, default_series_value)

  series_name_chooser = widgets.Dropdown(
      options=[(str(k), k) for k in sorted(dicoms.keys())],
      value=default_series_value,
      description='Choose the MRI series to visualize:',
      style={'description_width': 'initial'},
      layout=widgets.Layout(width='800px'),
  )
  # Slide through dicom image instances using a slide bar.
  instance_chooser = widgets.IntSlider(
      continuous_update=True,
      value=default_instance_value,
      min=1,
      max=max_instance_value,
      description='Image instance to display '
      + '(click on slider, then use left/right arrows):',
      style={'description_width': 'initial'},
      layout=series_name_chooser.layout,
  )
  vmin_chooser = widgets.IntSlider(
      continuous_update=True,
      value=default_vmin_value,
      min=MIN_COLOR_RANGE,
      max=MAX_COLOR_RANGE,
      description='Color range minimum:',
      style={'description_width': 'initial'},
      layout=widgets.Layout(width='300px'),
  )
  vmax_chooser = widgets.IntSlider(
      continuous_update=True,
      value=default_vmax_value,
      min=MIN_COLOR_RANGE,
      max=MAX_COLOR_RANGE,
      description='Color range maximum:',
      style={'description_width': 'initial'},
      layout=vmin_chooser.layout,
  )
  transpose_chooser = widgets.Checkbox(
      description='Whether to transpose the image.',
      style={'description_width': 'initial'},
      layout=vmin_chooser.layout,
  )
  fig_width_chooser = widgets.IntSlider(
      continuous_update=False,
      value=DEFAULT_IMAGE_WIDTH,
      min=MIN_IMAGE_WIDTH,
      max=MAX_IMAGE_WIDTH,
      description='Width of figure (height will be computed using input data):',
      style={'description_width': 'initial'},
      layout=vmin_chooser.layout,
  )

  viz_controls_ui = widgets.VBox(
      [
          widgets.HTML('<h3>Visualization controls</h3>'),
          series_name_chooser, instance_chooser,
          widgets.HBox([vmin_chooser, vmax_chooser]),
          widgets.HBox([transpose_chooser, fig_width_chooser]),
      ],
      layout=widgets.Layout(width='auto', border='solid 1px grey'),
  )
  viz_controls_output = widgets.interactive_output(
      dicom_animation,
      {
          'dicoms': widgets.fixed(dicoms),
          'series_name': series_name_chooser,
          'instance': instance_chooser,
          'vmin': vmin_chooser,
          'vmax': vmax_chooser,
          'transpose': transpose_chooser,
          'fig_width': fig_width_chooser,
          'title_prefix': widgets.fixed(os.path.basename(sample_mri)),
      },
  )

  def on_value_change(change):
    """Inner function to capture state being observed."""
    vmin_chooser.value, vmax_chooser.value = compute_color_range(dicoms, change['new'])
    instance_chooser.value, instance_chooser.max = compute_instance_range(dicoms, change['new'])

  # When the series changes, update the widgets to the proper ranges
  # for the series.
  series_name_chooser.observe(on_value_change, names='value')
  display(viz_controls_ui, viz_controls_output)


def compute_color_range(dicoms, series_name):
  """Compute the mean values for the color ranges of instances in the series."""
  vmin = np.mean([np.min(d.pixel_array) for d in dicoms[series_name]])
  vmax = np.mean([np.max(d.pixel_array) for d in dicoms[series_name]])
  return(vmin, vmax)


def compute_instance_range(dicoms, series_name):
  """Compute middle and max instances."""
  middle_instance = int(len(dicoms[series_name]) / 2)
  max_instance = len(dicoms[series_name])
  return(middle_instance, max_instance)


def dicom_animation(
    dicoms, series_name, instance, vmin, vmax, transpose,
    fig_width, title_prefix='',
):
  """Render one frame of a dicom animation.

  Args:
    dicoms: the dictionary DICOM series and instances lists
    series_name: the name of the series to be displayed
    instance: the particular instance to display
    vmin: minimum value for the color range
    vmax: maximum value for the color range
    transpose: whether or not to transpose the image
    fig_width: the desired width of the figure, note that height computed as
      the proportion of the width based on the data to be plotted
    title_prefix: text to display as the initial portion of the plot title
  """
  if len(dicoms[series_name]) < instance:
    dcm = dicoms[series_name][-1]
    print(f'Instance {str(instance)} not available for {series_name}, using final instance instead.')
  else:
    dcm = dicoms[series_name][instance - 1]
    if instance != dcm.InstanceNumber:
      # Notice invalid input, but don't throw an error.
      print(f'WARNING: Instance parameter {str(instance)} and dicom instance number {str(dcm.InstanceNumber)} do not match.')

  if transpose:
    height = dcm.pixel_array.T.shape[0]
    width = dcm.pixel_array.T.shape[1]
  else:
    height = dcm.pixel_array.shape[0]
    width = dcm.pixel_array.shape[1]

  fig_height = int(np.ceil(fig_width * (height/width)))

  _, ax = plt.subplots(figsize=(fig_width, fig_height), facecolor='beige')
  ax.imshow(dcm.pixel_array.T if transpose else dcm.pixel_array, cmap='gray', vmin=vmin, vmax=vmax)
  ax.set_title(
      title_prefix
      + ', Series: ' + dcm.SeriesDescription
      + ', Series Number: ' + str(dcm.SeriesNumber)
      + ', Instance: ' + str(dcm.InstanceNumber)
      + '\nColor range: ' + str(vmin) + '-' + str(vmax)
      + ', Transpose: ' + str(transpose)
      + ', Figure size:' + str(fig_width) + 'x' + str(fig_height),
      fontsize=fig_width,
  )
  ax.set_yticklabels([])
  ax.set_xticklabels([])
