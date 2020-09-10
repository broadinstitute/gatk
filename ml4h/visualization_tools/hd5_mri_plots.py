"""Methods for integration of plots of mri data processed to 3D tensors from within notebooks."""

from enum import Enum, auto
import os
import tempfile

import h5py
from IPython.display import display
from IPython.display import HTML
import ipywidgets as widgets
import matplotlib.pyplot as plt
from ml4h.runtime_data_defines import get_mri_hd5_folder
from ml4h.tensor_maps_by_hand import TMAPS
from ml4h.TensorMap import Interpretation
import numpy as np
import tensorflow as tf

# Discover applicable TMAPS.
CARDIAC_MRI_TMAP_NAMES = [k for k in TMAPS.keys() if ('_lax_' in k or '_sax_' in k) and TMAPS[k].axes() == 3]
CARDIAC_MRI_TMAP_NAMES.extend(
    [k for k in TMAPS.keys() if TMAPS[k].path_prefix == 'ukb_cardiac_mri' and TMAPS[k].axes() == 3],
)
LIVER_MRI_TMAP_NAMES = [k for k in TMAPS.keys() if TMAPS[k].path_prefix == 'ukb_liver_mri' and TMAPS[k].axes() == 3]
BRAIN_MRI_TMAP_NAMES = [k for k in TMAPS.keys() if TMAPS[k].path_prefix == 'ukb_brain_mri' and TMAPS[k].axes() == 3]
# This includes more than just MRI TMAPS, it is a best effort.
BEST_EFFORT_MRI_TMAP_NAMES = [k for k in TMAPS.keys() if TMAPS[k].interpretation == Interpretation.CONTINUOUS and TMAPS[k].axes() == 3]

MIN_IMAGE_WIDTH = 8
DEFAULT_IMAGE_WIDTH = 12
MAX_IMAGE_WIDTH = 24

MIN_COLOR_RANGE = 0
MAX_COLOR_RANGE = 6000


class PlotType(Enum):
  INTERACTIVE = auto()
  PANEL = auto()


class TensorMapCache:
  """Cache the tensor to display for reuse when re-plotting the same TMAP with different plot parameters."""

  def __init__(self, hd5, tmap_name):
    self.hd5 = hd5
    self.tmap_name = None
    self.tensor = None
    _ = self.get(tmap_name)

  def get(self, tmap_name):
    if self.tmap_name != tmap_name:
      self.tensor = TMAPS[tmap_name].tensor_from_file(TMAPS[tmap_name], self.hd5)
      self.tmap_name = tmap_name
    return self.tensor


def choose_cardiac_mri_tmap(sample_id, folder=None, tmap_name='cine_lax_4ch_192', default_tmap_names=CARDIAC_MRI_TMAP_NAMES):
  choose_mri_tmap(sample_id, folder, tmap_name, default_tmap_names)


def choose_brain_mri_tmap(sample_id, folder=None, tmap_name='t2_flair_sag_p2_1mm_fs_ellip_pf78_1', default_tmap_names=BRAIN_MRI_TMAP_NAMES):
  choose_mri_tmap(sample_id, folder, tmap_name, default_tmap_names)


def choose_liver_mri_tmap(sample_id, folder=None, tmap_name='liver_shmolli_segmented', default_tmap_names=LIVER_MRI_TMAP_NAMES):
  choose_mri_tmap(sample_id, folder, tmap_name, default_tmap_names)


def choose_mri_tmap(sample_id, folder=None, tmap_name=None, default_tmap_names=BEST_EFFORT_MRI_TMAP_NAMES):
  """Render widgets and plots for MRI tensors.

  Args:
    sample_id: The id of the sample to retrieve.
    folder: The local or Cloud Storage folder under which the files reside.
    tmap_name: The TMAP name for the 3D MRI tensor to visualize.
    default_tmap_names: Other TMAP names to offer for visualization, if present in the hd5.

  Returns:
    ipywidget or HTML upon error.
  """
  if folder is None:
    folder = get_mri_hd5_folder(sample_id)

  with tempfile.TemporaryDirectory() as tmpdirname:
    sample_hd5 = str(sample_id) + '.hd5'
    local_path = os.path.join(tmpdirname, sample_hd5)
    try:
      tf.io.gfile.copy(src=os.path.join(folder, sample_hd5), dst=local_path)
      hd5 = h5py.File(local_path, mode='r')
    except (tf.errors.NotFoundError, tf.errors.PermissionDeniedError) as e:
      return HTML(f'''
      <div class="alert alert-block alert-danger">
      <b>Warning:</b> MRI HD5 file not available for sample {sample_id} in folder {folder}:
      <hr><p><pre>{e.message}</pre></p>
      Use the <kbd>folder</kbd> parameter to read HD5s from a different local directory or Cloud Storage bucket.
      </div>''')

    sample_tmap_names = []
    # Add the passed tmap_name parameter, if it is present in this hd5.
    if tmap_name:
      if TMAPS[tmap_name].hd5_key_guess() in hd5:
        if len(TMAPS[tmap_name].shape) == 3:
          sample_tmap_names.append(tmap_name)
        else:
          print(f'{tmap_name} is not a 3D tensor, skipping it')
      else:
        print(f'{tmap_name} is not available in {sample_id}')
    # Also discover applicable TMAPS for this particular sample's HD5 file.
    sample_tmap_names.extend(
        sorted(set([k for k in default_tmap_names if TMAPS[k].hd5_key_guess() in hd5])),
    )

    if not sample_tmap_names:
      return HTML(f'''<div class="alert alert-block alert-danger">
      Neither {tmap_name} nor any of {default_tmap_names} are present in this HD5 for sample {sample_id} in {folder}.
      Use the tmap_name parameter to try a different TMAP or the folder parameter to try a different hd5 for the sample.
      </div>''')

    default_tmap_name_value = sample_tmap_names[0]
    # Display the middle instance by default in the interactive view.
    default_instance_value, max_instance_value = compute_instance_range(default_tmap_name_value)
    default_vmin_value, default_vmax_value = compute_color_range(hd5, default_tmap_name_value)

    tmap_name_chooser = widgets.Dropdown(
        options=sample_tmap_names,
        value=default_tmap_name_value,
        description='Choose the MRI tensor TMAP name to visualize:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='900px'),
    )
    fig_width_chooser = widgets.IntSlider(
        continuous_update=False,
        value=DEFAULT_IMAGE_WIDTH,
        min=MIN_IMAGE_WIDTH,
        max=MAX_IMAGE_WIDTH,
        description='Desired width of figure (height computed using input data)',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='600px'),
    )
    transpose_chooser = widgets.Checkbox(
        description='Whether to transpose the images.',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='300px'),
    )
    flip_chooser = widgets.Checkbox(
        description='Whether to flip the images.',
        style={'description_width': 'initial'},
        layout=transpose_chooser.layout,
    )
    plot_type_chooser = widgets.RadioButtons(
        options={'interactive animation': PlotType.INTERACTIVE, 'panel grid': PlotType.PANEL},
        description='Plot type',
        style={'description_width': 'initial'},
        layout=transpose_chooser.layout,
    )
    instance_chooser = widgets.IntSlider(
        continuous_update=True,
        value=default_instance_value,
        min=1,
        max=max_instance_value,
        description='Instance to display (click on slider, then use left/right arrows):',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='900px'),
    )
    color_range_chooser = widgets.IntRangeSlider(
        continuous_update=False,
        value=[default_vmin_value, default_vmax_value],
        min=MIN_COLOR_RANGE,
        max=MAX_COLOR_RANGE,
        description='Color range:',
        style={'description_width': 'initial'},
        layout=fig_width_chooser.layout,
    )
    viz_controls_ui = widgets.VBox(
        [
            widgets.HTML('<h3>Visualization controls</h3>'),
            tmap_name_chooser,
            widgets.HBox([transpose_chooser, fig_width_chooser]),
            widgets.HBox([flip_chooser, color_range_chooser]),
            widgets.HBox([plot_type_chooser, instance_chooser]),
        ],
        layout=widgets.Layout(width='auto', border='solid 1px grey'),
    )
    tmap_cache = TensorMapCache(hd5=hd5, tmap_name=tmap_name_chooser.value)
    viz_controls_output = widgets.interactive_output(
        plot_mri_tmap,
        {
            'sample_id': widgets.fixed(sample_id),
            'tmap_cache': widgets.fixed(tmap_cache),
            'tmap_name': tmap_name_chooser,
            'plot_type': plot_type_chooser,
            'instance': instance_chooser,
            'color_range': color_range_chooser,
            'transpose': transpose_chooser,
            'flip': flip_chooser,
            'fig_width': fig_width_chooser,
        },
    )

    def on_tmap_value_change(change):
      """When the TMAP name changes, update the widgets to the proper ranges for the new TMAP name."""
      color_range_chooser.value = compute_color_range(hd5, change['new'])
      instance_chooser.value, instance_chooser.max = compute_instance_range(change['new'])

    def on_plot_type_change(change):
      """Only show controls applicable to the plot type."""
      if change['new'] == PlotType.INTERACTIVE:
        instance_chooser.layout.visibility = 'visible'
      else:
        instance_chooser.layout.visibility = 'hidden'

    tmap_name_chooser.observe(on_tmap_value_change, names='value')
    plot_type_chooser.observe(on_plot_type_change, names='value')
    display(viz_controls_ui, viz_controls_output)


def compute_color_range(hd5, tmap_name):
  """Compute the mean values for the color ranges of instances in the MRI series."""
  mri_tensor = TMAPS[tmap_name].tensor_from_file(TMAPS[tmap_name], hd5)
  vmin = np.mean([np.min(mri_tensor[:, :, i]) for i in range(0, mri_tensor.shape[2])])
  vmax = np.mean([np.max(mri_tensor[:, :, i]) for i in range(0, mri_tensor.shape[2])])
  return[vmin, vmax]


def compute_instance_range(tmap_name):
  """Compute middle and max instances."""
  middle_instance = int(TMAPS[tmap_name].shape[2] / 2)
  max_instance = TMAPS[tmap_name].shape[2]
  return(middle_instance, max_instance)


def plot_mri_tmap(sample_id, tmap_cache, tmap_name, plot_type, instance, color_range, transpose, flip, fig_width):
  """Visualize the applicable MRI series within this HD5 file.

  Args:
    sample_id: The local or Cloud Storage path to the MRI file.
    tmap_cache: The cache from which to retrieve the tensor to be plotted.
    tmap_name: The name of the chosen TMAP for the MRI series.
    plot_type: Whether to display instances interactively or in a panel view.
    instance: The particular instance to display, if interactive.
    color_range: Array of minimum and maximum value for the color range.
    transpose: Whether to transpose the images.
    flip: Whether to flip the image on its vertical axis
    fig_width: The desired width of the figure. Note that height computed as
      the proportion of the width based on the data to be plotted.

  Returns:
    The plot or a notebook-friendly error message.
  """
  title_prefix = f'{tmap_name} from MRI {sample_id}'
  mri_tensor = tmap_cache.get(tmap_name)
  if plot_type == PlotType.INTERACTIVE:
    plot_mri_tensor_as_animation(
        mri_tensor=mri_tensor,
        instance=instance,
        vmin=color_range[0],
        vmax=color_range[1],
        transpose=transpose,
        flip=flip,
        fig_width=fig_width,
        title_prefix=title_prefix,
    )
  elif plot_type == PlotType.PANEL:
    # Note: this print statement causes the current image to go away while the new one is rendering.
    # Do this intensionally for the panel plots, because they are a bit slower to render. It serves
    # a purpose similar to a progress bar.
    print(f'Rendering {title_prefix} . . .')
    plot_mri_tensor_as_panels(
        mri_tensor=mri_tensor,
        vmin=color_range[0],
        vmax=color_range[1],
        transpose=transpose,
        flip=flip,
        fig_width=fig_width,
        title_prefix=title_prefix,
    )
  else:
    return HTML(f'''<div class="alert alert-block alert-danger">Invalid plot type: {plot_type}</div>''')


def plot_mri_tensor_as_panels(mri_tensor, vmin, vmax, transpose=False, flip=False, fig_width=DEFAULT_IMAGE_WIDTH, title_prefix=''):
  """Visualize an MRI series from a 3D tensor as a panel of static plots.

  Args:
    mri_tensor: The MRI 3D tensor.
    vmin: Minimum value for the color range.
    vmax: Maximum value for the color range.
    transpose: Whether to transpose the image.
    flip: Whether to flip the image on its vertical axis.
    fig_width: The desired width of the figure, note that height computed as
      the proportion of the width based on the data to be plotted.
    title_prefix: Text to display as the initial portion of the plot title.
  """
  cols = 5
  rows = int(np.ceil(mri_tensor.shape[2] / 5.0))
  if transpose:
    height = mri_tensor.shape[1]
    width = mri_tensor.shape[0]
  else:
    height = mri_tensor.shape[0]
    width = mri_tensor.shape[1]

  fig_height = int(np.ceil(fig_width * ((rows * height)/(cols * width))))
  fig, axes = plt.subplots(rows, cols, figsize=(fig_width, fig_height), facecolor='beige')
  for i in range(0, mri_tensor.shape[2]):
    col = i % cols
    row = i // cols
    pixels = mri_tensor[:, :, i]
    if transpose:
      pixels = pixels.T
    if flip:
      pixels = np.flip(pixels, 1)
    axes[row, col].imshow(pixels, cmap='gray', vmin=vmin, vmax=vmax)
    axes[row, col].set_yticklabels([])
    axes[row, col].set_xticklabels([])
  fig.suptitle(
      f'{title_prefix}\nColor range: {vmin}-{vmax}, Transpose: {transpose}, Flip: {flip}, Figure size:{fig_width}x{fig_height}',
      fontsize=fig_width,
  )
  fig.subplots_adjust(
      top=0.96,    # the top of the subplots of the figure
      wspace=0.1,  # the amount of width reserved for space between subplots,
                   # expressed as a fraction of the average axis width
      hspace=0.1,  # the amount of height reserved for space between subplots,
                   # expressed as a fraction of the average axis height
  )


def plot_mri_tensor_as_animation(mri_tensor, instance, vmin, vmax, transpose=False, flip=False, fig_width=DEFAULT_IMAGE_WIDTH, title_prefix=''):
  """Visualize an MRI series from a 3D tensor as an animation rendered one panel at a time.

  Args:
    mri_tensor: The MRI 3D tensor.
    instance: The particular instance to display.
    vmin: Minimum value for the color range.
    vmax: Maximum value for the color range.
    transpose: Whether to transpose the image.
    flip: Whether to flip the image on its vertical axis.
    fig_width: The desired width of the figure, note that height computed as
      the proportion of the width based on the data to be plotted.
    title_prefix: Text to display as the initial portion of the plot title.
  """
  if mri_tensor.shape[2] < instance:
    pixels = mri_tensor[:, :, -1]
    print(f'Instance {str(instance)} not available for {title_prefix}, using final instance instead.')
  else:
    pixels = mri_tensor[:, :, instance - 1]

  if transpose:
    pixels = pixels.T
  if flip:
    pixels = np.flip(pixels, 1)

  height = pixels.shape[0]
  width = pixels.shape[1]
  fig_height = int(np.ceil(fig_width * (height/width)))

  _, ax = plt.subplots(figsize=(fig_width, fig_height), facecolor='beige')
  ax.imshow(pixels, cmap='gray', vmin=vmin, vmax=vmax)
  ax.set_title(
      f'{title_prefix}, Instance: {instance}\nColor range: {vmin}-{vmax}, Transpose: {transpose}, Flip: {flip}, Figure size:{fig_width}x{fig_height}',
      fontsize=fig_width,
  )
  ax.set_yticklabels([])
  ax.set_xticklabels([])
