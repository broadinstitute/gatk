"""Methods for integration of plots of mri data processed to 3D tensors from within notebooks."""

from collections import OrderedDict
from enum import Enum, auto
import os
import tempfile
from typing import Any, Dict, List, Optional, Tuple, Union

from IPython.display import display
from IPython.display import HTML
import numpy as np
import h5py
import ipywidgets as widgets
import matplotlib.pyplot as plt
from ml4h.runtime_data_defines import get_mri_hd5_folder
import ml4h.tensormap.ukb.mri as ukb_mri
import ml4h.tensormap.ukb.mri_vtk as ukb_mri_vtk
from ml4h.TensorMap import Interpretation, TensorMap
import tensorflow as tf

# Discover applicable TensorMaps.
MRI_TMAPS = {
    key: value for key, value in ukb_mri.__dict__.items() if isinstance(value, TensorMap)
    and value.interpretation == Interpretation.CONTINUOUS and value.axes() == 3
}
MRI_TMAPS.update(
    {
        key: value for key, value in ukb_mri_vtk.__dict__.items()
        if isinstance(value, TensorMap) and value.interpretation == Interpretation.CONTINUOUS and value.axes() == 3
    },
)

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

  def __init__(self, hd5: Dict[str, Any], tmap: TensorMap):
    self.hd5 = hd5
    self.tmap: Optional[TensorMap] = None
    self.tensor = None
    _ = self.get(tmap)

  def get(self, tmap: TensorMap) -> np.array:
    if self.tmap != tmap:
      self.tensor = tmap.tensor_from_file(tmap, self.hd5)
      self.tmap = tmap
    return self.tensor


def choose_mri_tmap(
    sample_id: Union[int, str], folder: Optional[str] = None, tmap: Optional[TensorMap] = None,
    default_tmaps: Dict[str, TensorMap] = MRI_TMAPS,
) -> None:
  """Render widgets and plots for MRI tensors.

  Args:
    sample_id: The id of the sample to retrieve.
    folder: The local or Cloud Storage folder under which the files reside.
    tmap: The TensorMap for the 3D MRI tensor to visualize.
    default_tmaps: Other TensorMaps to offer for visualization, if present in the hd5.
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
      display(
          HTML(f'''<div class="alert alert-block alert-danger">
      <b>Warning:</b> MRI HD5 file not available for sample {sample_id} in folder {folder}:
      <hr><p><pre>{e.message}</pre></p>
      Use the <kbd>folder</kbd> parameter to read HD5s from a different local directory or Cloud Storage bucket.
      </div>'''),
      )
      return

    sample_tmaps = OrderedDict()
    # Add the passed tmap parameter, if it is present in this hd5.
    if tmap:
      if tmap.hd5_key_guess() in hd5:
        if len(tmap.shape) == 3:
          sample_tmaps[tmap.name] = tmap
        else:
          print(f'{tmap} is not a 3D tensor, skipping it')
      else:
        print(f'{tmap} is not available in {sample_id}')
    # Also discover applicable TensorMaps for this particular sample's HD5 file.
    sample_tmaps.update({n: t for n, t in sorted(default_tmaps.items(), key=lambda t: t[0]) if t.hd5_key_guess() in hd5})

    if not sample_tmaps:
      display(
          HTML(f'''<div class="alert alert-block alert-danger">
      Neither {tmap.name} nor any of {default_tmaps.keys()} are present in this HD5 for sample {sample_id} in {folder}.
      Use the tmap parameter to try a different TensorMap or the folder parameter to try a different hd5 for the sample.
      </div>'''),
      )
      return

    default_tmap_value = next(iter(sample_tmaps.values()))
    # Display the middle instance by default in the interactive view.
    default_instance_value, max_instance_value = compute_instance_range(default_tmap_value)
    default_vmin_value, default_vmax_value = compute_color_range(hd5, default_tmap_value)

    tmap_chooser = widgets.Dropdown(
        options=sample_tmaps,
        value=default_tmap_value,
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
            tmap_chooser,
            widgets.HBox([transpose_chooser, fig_width_chooser]),
            widgets.HBox([flip_chooser, color_range_chooser]),
            widgets.HBox([plot_type_chooser, instance_chooser]),
        ],
        layout=widgets.Layout(width='auto', border='solid 1px grey'),
    )
    tmap_cache = TensorMapCache(hd5=hd5, tmap=tmap_chooser.value)
    viz_controls_output = widgets.interactive_output(
        plot_mri_tmap,
        {
            'sample_id': widgets.fixed(sample_id),
            'tmap_cache': widgets.fixed(tmap_cache),
            'tmap': tmap_chooser,
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

    tmap_chooser.observe(on_tmap_value_change, names='value')
    plot_type_chooser.observe(on_plot_type_change, names='value')
    display(viz_controls_ui, viz_controls_output)


def compute_color_range(hd5: Dict[str, Any], tmap: TensorMap) -> List[int]:
  """Compute the mean values for the color ranges of instances in the MRI series."""
  mri_tensor = tmap.tensor_from_file(tmap, hd5)
  vmin = np.mean([np.min(mri_tensor[:, :, i]) for i in range(0, mri_tensor.shape[2])])
  vmax = np.mean([np.max(mri_tensor[:, :, i]) for i in range(0, mri_tensor.shape[2])])
  return [vmin, vmax]


def compute_instance_range(tmap: TensorMap) -> Tuple[int, int]:
  """Compute middle and max instances."""
  middle_instance = int(tmap.shape[2] / 2)
  max_instance = tmap.shape[2]
  return (middle_instance, max_instance)


def plot_mri_tmap(
    sample_id: Union[int, str], tmap_cache: TensorMapCache, tmap: TensorMap, plot_type: PlotType,
    instance: int, color_range: Tuple[int, int], transpose: bool, flip: bool, fig_width: int,
) -> None:
  """Visualize the applicable MRI series within this HD5 file.

  Args:
    sample_id: The local or Cloud Storage path to the MRI file.
    tmap_cache: The cache from which to retrieve the tensor to be plotted.
    tmap: The chosen TensorMap for the MRI series.
    plot_type: Whether to display instances interactively or in a panel view.
    instance: The particular instance to display, if interactive.
    color_range: Array of minimum and maximum value for the color range.
    transpose: Whether to transpose the images.
    flip: Whether to flip the image on its vertical axis
    fig_width: The desired width of the figure. Note that height computed as
      the proportion of the width based on the data to be plotted.
  """
  title_prefix = f'{tmap.name} from MRI {sample_id}'
  mri_tensor = tmap_cache.get(tmap)
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
    HTML(f'''<div class="alert alert-block alert-danger">Invalid plot type: {plot_type}</div>''')


def plot_mri_tensor_as_panels(
    mri_tensor: np.array, vmin: int, vmax: int, transpose: bool = False, flip: bool = False,
    fig_width: int = DEFAULT_IMAGE_WIDTH, title_prefix: str = '',
) -> None:
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
      f'{title_prefix}\nColor range: {vmin}-{vmax}, Transpose: {transpose}, Flip: {flip}, Figure size:{fig_width}x{fig_height}',  # pylint: disable=line-too-long
      fontsize=fig_width,
  )
  fig.subplots_adjust(
      top=0.96,    # the top of the subplots of the figure
      wspace=0.1,  # the amount of width reserved for space between subplots,
                   # expressed as a fraction of the average axis width
      hspace=0.1,  # the amount of height reserved for space between subplots,
                   # expressed as a fraction of the average axis height
  )


def plot_mri_tensor_as_animation(
    mri_tensor: np.array, instance: int, vmin: int, vmax: int,
    transpose: bool = False, flip: bool = False,
    fig_width: int = DEFAULT_IMAGE_WIDTH, title_prefix: str = '',
) -> None:
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
      f'{title_prefix}, Instance: {instance}\nColor range: {vmin}-{vmax}, Transpose: {transpose}, Flip: {flip}, Figure size:{fig_width}x{fig_height}',  # pylint: disable=line-too-long
      fontsize=fig_width,
  )
  ax.set_yticklabels([])
  ax.set_xticklabels([])
