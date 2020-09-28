"""Methods for integration of DICOM plots within notebooks."""

import collections
import os
import tempfile
from typing import Dict, List, Optional, Tuple, Union
import zipfile

from IPython.display import display
from IPython.display import HTML
import numpy as np
import ipywidgets as widgets
import matplotlib.pyplot as plt
from ml4h.runtime_data_defines import get_cardiac_mri_folder
import pydicom
from scipy.ndimage.morphology import binary_closing
from scipy.ndimage.morphology import binary_erosion
import tensorflow as tf

# Constants for use with 'CINE_segmented_SAX_InlineVF'
# TODO(deflaux) move these constants into ml4h/defines.py and then import the
# ml4h package.
MRI_FRAMES = 50
MRI_MIN_RADIUS = 2
MRI_MAX_MYOCARDIUM = 20
MRI_BIG_RADIUS_FACTOR = 0.9
MRI_SMALL_RADIUS_FACTOR = 0.19
MRI_SEGMENTED_CHANNEL_MAP = {'background': 0, 'ventricle': 1, 'myocardium': 2}


def _is_mitral_valve_segmentation(d: pydicom.FileDataset) -> bool:
  """Determine whether a DICOM has mitral valve segmentation.

  This is used for visualization of CINE_segmented_SAX_InlineVF.

  Args:
    d: the DICOM file

  Returns:
    Whether or not the DICOM has mitral valve segmentation
  """
  return d.SliceThickness == 6


def _get_overlay_from_dicom(d: pydicom.FileDataset) -> Tuple[int, int, int]:
  """Get an overlay from a DICOM file.

  Morphological operators are used to transform the pixel outline of the
  myocardium to the labeled pixel masks for myocardium and left ventricle. This
  is used for visualization of CINE_segmented_SAX_InlineVF.

  Args:
    d: the DICOM file

  Returns:
    Raw overlay array with myocardium outline, anatomical mask (a pixel
    mask with 0 for background 1 for myocardium and 2 for ventricle), and
    ventrical pixels.
  """
  i_overlay = 0
  dicom_tag = 0x6000 + 2 * i_overlay
  overlay_raw = d[dicom_tag, 0x3000].value
  rows = d[dicom_tag, 0x0010].value  # rows = 512
  cols = d[dicom_tag, 0x0011].value  # cols = 512
  overlay_frames = d[dicom_tag, 0x0015].value
  bits_allocated = d[dicom_tag, 0x0100].value

  np_dtype = np.dtype('uint8')
  length_of_pixel_array = len(overlay_raw)
  expected_length = rows * cols
  if bits_allocated == 1:
    expected_bit_length = expected_length
    bit = 0
    overlay = np.ndarray(shape=(length_of_pixel_array * 8), dtype=np_dtype)
    for byte in overlay_raw:
      for bit in range(bit, bit + 8):
        overlay[bit] = byte & 0b1
        byte >>= 1
      bit += 1
    overlay = overlay[:expected_bit_length]
  if overlay_frames != 1:
    raise ValueError(f'DICOM has {overlay_frames} overlay frames, but only one expected.')
  overlay = overlay.reshape(rows, cols)
  idx = np.where(overlay == 1)
  min_pos = (np.min(idx[0]), np.min(idx[1]))
  max_pos = (np.max(idx[0]), np.max(idx[1]))
  short_side = min((max_pos[0] - min_pos[0]), (max_pos[1] - min_pos[1]))
  small_radius = max(MRI_MIN_RADIUS, short_side * MRI_SMALL_RADIUS_FACTOR)
  big_radius = max(MRI_MIN_RADIUS+1, short_side * MRI_BIG_RADIUS_FACTOR)
  small_structure = _unit_disk(small_radius)
  m1 = binary_closing(overlay, small_structure).astype(np.int)
  big_structure = _unit_disk(big_radius)
  m2 = binary_closing(overlay, big_structure).astype(np.int)
  anatomical_mask = m1 + m2
  ventricle_pixels = np.count_nonzero(anatomical_mask == MRI_SEGMENTED_CHANNEL_MAP['ventricle'])
  myocardium_pixels = np.count_nonzero(anatomical_mask == MRI_SEGMENTED_CHANNEL_MAP['myocardium'])
  if ventricle_pixels == 0 and myocardium_pixels > MRI_MAX_MYOCARDIUM:
    erode_structure = _unit_disk(small_radius*1.5)
    anatomical_mask = anatomical_mask - binary_erosion(m1, erode_structure).astype(np.int)
    ventricle_pixels = np.count_nonzero(anatomical_mask == MRI_SEGMENTED_CHANNEL_MAP['ventricle'])
  return overlay, anatomical_mask, ventricle_pixels


def _unit_disk(r: int) -> np.ndarray:
  """Get the unit disk for a radius.

  This is used for visualization of CINE_segmented_SAX_InlineVF.

  Args:
    r: the radius

  Returns:
    The unit disk.
  """
  y, x = np.ogrid[-r: r + 1, -r: r + 1]
  return (x ** 2 + y ** 2 <= r ** 2).astype(np.int)


def plot_cardiac_long_axis(
    b_series: List[pydicom.FileDataset], sides: int = 7, fig_width: int = 18, title_prefix: str = '',
) -> None:
  """Visualize CINE_segmented_SAX_InlineVF series.

  Args:
    b_series: the DICOM
    sides: the number of sides to display
    fig_width: the desired width of the figure, note that height computed as
      the proportion of the width based on the data to be plotted
    title_prefix: text to display as the initial portion of the plot title
  """
  height = b_series[0].pixel_array.shape[0]
  width = b_series[0].pixel_array.shape[1]
  fig_height = int(np.ceil(fig_width * (height/width)))
  fig, axes = plt.subplots(sides, sides, figsize=(fig_width, fig_height), facecolor='beige')
  for dcm in b_series:
    idx = (dcm.InstanceNumber-1) % MRI_FRAMES
    if idx >= sides*sides:
      continue
    if _is_mitral_valve_segmentation(dcm):
      axes[idx%sides, idx//sides].imshow(dcm.pixel_array)
    else:
      try:
        _, anatomical_mask, _ = _get_overlay_from_dicom(dcm)
        axes[idx%sides, idx//sides].imshow(
            np.ma.masked_where(anatomical_mask == 2, dcm.pixel_array),
            cmap='gray',
            vmin=np.min(dcm.pixel_array),
            vmax=np.max(dcm.pixel_array),
        )
      except KeyError:
        axes[idx%sides, idx//sides].imshow(
            dcm.pixel_array,
            cmap='gray',
            vmin=np.min(dcm.pixel_array),
            vmax=np.max(dcm.pixel_array),
        )
    axes[idx%sides, idx//sides].set_yticklabels([])
    axes[idx%sides, idx//sides].set_xticklabels([])

  fig.suptitle(
      title_prefix + ', Number of sides: ' + str(sides)
      + ', Figure size:' + str((fig_width, fig_height)),
      fontsize=fig_width,
  )
  fig.subplots_adjust(
      top=0.96,    # the top of the subplots of the figure
      wspace=0.1,  # the amount of width reserved for space between subplots,
                   # expressed as a fraction of the average axis width
      hspace=0.1,  # the amount of height reserved for space between subplots,
                   # expressed as a fraction of the average axis height
  )


def plot_cardiac_short_axis(
    series: List[pydicom.FileDataset], transpose: bool = False, fig_width: int = 18,
    title_prefix: str = '',
) -> None:
  """Visualize CINE_segmented_LAX series.

  Args:
    series: the DICOM
    transpose: whether or not to transpose the image
    fig_width: the desired width of the figure, note that height computed as
      the proportion of the width based on the data to be plotted
    title_prefix: text to display as the initial portion of the plot title
  """
  cols = 5
  rows = 10
  if transpose:
    height = series[0].pixel_array.T.shape[0]
    width = series[0].pixel_array.T.shape[1]
  else:
    height = series[0].pixel_array.shape[0]
    width = series[0].pixel_array.shape[1]

  fig_height = int(np.ceil(fig_width * ((rows * height)/(cols * width))))
  fig, axes = plt.subplots(rows, cols, figsize=(fig_width, fig_height), facecolor='beige')
  for dcm in series:
    col = (dcm.InstanceNumber-1)%cols
    row = (dcm.InstanceNumber-1)//cols
    if transpose:
      axes[row, col].imshow(
          dcm.pixel_array.T,
          cmap='gray',
          vmin=np.min(dcm.pixel_array),
          vmax=np.max(dcm.pixel_array),
      )
    else:
      axes[row, col].imshow(
          dcm.pixel_array,
          cmap='gray',
          vmin=np.min(dcm.pixel_array),
          vmax=np.max(dcm.pixel_array),
      )
    axes[row, col].set_yticklabels([])
    axes[row, col].set_xticklabels([])
  fig.suptitle(
      title_prefix + ', Transpose: ' + str(transpose)
      + ', Figure size:' + str((fig_width, fig_height)),
      fontsize=fig_width,
  )
  fig.subplots_adjust(
      top=0.96,    # the top of the subplots of the figure
      wspace=0.1,  # the amount of width reserved for space between subplots,
                   # expressed as a fraction of the average axis width
      hspace=0.1,  # the amount of height reserved for space between subplots,
                   # expressed as a fraction of the average axis height
  )


def plot_mri_series(
    sample_mri: str, dicoms: Dict[str, pydicom.FileDataset], series_name: str, sax_sides: int,
    lax_transpose: bool, fig_width: int,
) -> None:
  """Visualize the applicable series within this DICOM.

  Args:
    sample_mri: The local or Cloud Storage path to the MRI file.
    dicoms: A dictionary of DICOMs.
    series_name: The name of the chosen series.
    sax_sides: How many sides to display for CINE_segmented_SAX_InlineVF.
    lax_transpose: Whether to transpose when plotting CINE_segmented_LAX.
    fig_width: The desired width of the figure. Note that height computed as
      the proportion of the width based on the data to be plotted.

  """
  title_prefix = f'{dicoms[series_name][0].SeriesDescription}  from MRI {os.path.basename(sample_mri)}'
  print(f'Rendering {title_prefix}.')
  if 'cine_segmented_lax' in series_name:
    plot_cardiac_short_axis(
        dicoms[series_name],
        transpose=lax_transpose,
        fig_width=fig_width,
        title_prefix=title_prefix,
    )
  elif 'cine_segmented_sax_inlinevf' in series_name:
    plot_cardiac_long_axis(
        dicoms[series_name],
        sides=sax_sides,
        fig_width=fig_width,
        title_prefix=title_prefix,
    )
  else:
    print(f'Visualization not currently implemented for {series_name}.')


def choose_mri_series(sample_mri: str) -> None:
  """Render widgets and plots for cardiac MRIs.

  Visualization is supported for CINE_segmented_SAX_InlineVF series and
  CINE_segmented_LAX series.

  Args:
    sample_mri: The local or Cloud Storage path to the MRI file.
  """
  with tempfile.TemporaryDirectory() as tmpdirname:
    local_path = os.path.join(tmpdirname, os.path.basename(sample_mri))
    try:
      tf.io.gfile.copy(src=sample_mri, dst=local_path)
      with zipfile.ZipFile(local_path, 'r') as zip_ref:
        zip_ref.extractall(tmpdirname)
    except (tf.errors.NotFoundError, tf.errors.PermissionDeniedError) as e:
      display(
          HTML(f'''<div class="alert alert-block alert-danger">
      <b>Warning:</b> Cardiac MRI not available for sample {os.path.basename(sample_mri)}:
      <hr><p><pre>{e.message}</pre></p>
      </div>'''),
      )
      return

    filtered_dicoms = collections.defaultdict(list)
    series_descriptions = []
    for dcm_file in os.listdir(tmpdirname):
      if not dcm_file.endswith('.dcm'):
        continue
      dcm = pydicom.read_file(os.path.join(tmpdirname, dcm_file))
      series_descriptions.append(dcm.SeriesDescription)
      if 'cine_segmented_lax' in dcm.SeriesDescription.lower():
        filtered_dicoms[dcm.SeriesDescription.lower()].append(dcm)
      if dcm.SeriesDescription.lower() == 'cine_segmented_sax_inlinevf':
        cur_angle = (dcm.InstanceNumber - 1) // MRI_FRAMES
        filtered_dicoms[f'{dcm.SeriesDescription.lower()}_angle_{str(cur_angle)}'].append(dcm)

    print(f'{os.path.basename(sample_mri)} contains: {str(sorted(set(series_descriptions)))}.')

    if filtered_dicoms:
      series_name_chooser = widgets.Dropdown(
          options=sorted(list(filtered_dicoms.keys())),
          value=sorted(list(filtered_dicoms.keys()))[0],
          description='Choose the MRI series to visualize:',
          style={'description_width': 'initial'},
          layout=widgets.Layout(width='800px'),
      )
      fig_width_chooser = widgets.IntSlider(
          continuous_update=False,
          value=18,
          min=8,
          description='Desired width of figure (height will be computed using input data)',
          style={'description_width': 'initial'},
          layout=widgets.Layout(width='800px'),
      )
      sax_sides_chooser = widgets.IntSlider(
          continuous_update=False,
          value=7,
          min=1,
          description='How many sides to display for CINE_segmented_SAX_InlineVF',
          style=fig_width_chooser.style,
          layout=fig_width_chooser.layout,
      )
      lax_transpose_chooser = widgets.Checkbox(
          description='Whether to transpose the images when plotting CINE_segmented_LAX',
          style=fig_width_chooser.style,
          layout=fig_width_chooser.layout,
      )
      viz_controls_ui = widgets.VBox(
          [
              widgets.HTML('<h3>Visualization controls</h3>'), series_name_chooser,
              fig_width_chooser, sax_sides_chooser, lax_transpose_chooser,
          ],
          layout=widgets.Layout(width='auto', border='solid 1px grey'),
      )
      viz_controls_output = widgets.interactive_output(
          plot_mri_series,
          {
              'sample_mri': widgets.fixed(sample_mri),
              'dicoms': widgets.fixed(filtered_dicoms),
              'series_name': series_name_chooser,
              'sax_sides': sax_sides_chooser,
              'lax_transpose': lax_transpose_chooser,
              'fig_width': fig_width_chooser,
          },
      )
      display(viz_controls_ui, viz_controls_output)
    else:
      display(
          HTML(f'''<div class="alert alert-block alert-warning">
      Neither CINE_segmented_SAX_InlineVF nor CINE_segmented_LAX available in MRI for sample {os.path.basename(sample_mri)}.
      Try a different MRI.
      </div>'''),
      )


def choose_cardiac_mri(sample_id: Union[int, str], folder: Optional[str] = None) -> None:
  """Render widget to choose the cardiac MRI to plot.

  Args:
    sample_id: The id of the ECG sample to retrieve.
    folder: The local or Cloud Storage folder under which the files reside.
  """
  if folder is None:
    folder = get_cardiac_mri_folder(sample_id)

  sample_mri_glob = str(sample_id) + '_*.zip'
  try:
    sample_mris = tf.io.gfile.glob(pattern=os.path.join(folder, sample_mri_glob))
  except (tf.errors.NotFoundError, tf.errors.PermissionDeniedError) as e:
    display(
        HTML(f'''<div class="alert alert-block alert-danger">
    <b>Warning:</b> Cardiac MRI not available for sample {sample_id} in {folder}:
    <hr><p><pre>{e.message}</pre></p>
    Use the <kbd>folder</kbd> parameter to read DICOMs from a different local directory or Cloud Storage bucket.
    </div>'''),
    )
    return

  if not sample_mris:
    display(
        HTML(f'''<div class="alert alert-block alert-danger">
    <b>Warning:</b> Cardiac MRI DICOM not available for sample {sample_id} in {folder}.<br>
    Use the <kbd>folder</kbd> parameter to read DICOMs from a different local directory or Cloud Storage bucket.
    </div>'''),
    )
    return

  mri_chooser = widgets.Dropdown(
      options=[(os.path.basename(mri), mri) for mri in sample_mris],
      value=sample_mris[0],
      description=f'Choose an MRI to visualize for sample {sample_id}:',
      style={'description_width': 'initial'},
      layout=widgets.Layout(width='800px'),
  )
  file_controls_ui = widgets.VBox(
      [widgets.HTML('<h3>File controls</h3>'), mri_chooser],
      layout=widgets.Layout(width='auto', border='solid 1px grey'),
  )
  file_controls_output = widgets.interactive_output(
      choose_mri_series, {'sample_mri': mri_chooser},
  )
  display(file_controls_ui, file_controls_output)
