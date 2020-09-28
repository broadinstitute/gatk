"""Methods for Altair-based interactive plots for ECG data for use within notebooks."""

import os
import tempfile
from typing import Optional, Union

from IPython.display import HTML
import altair as alt  # Interactive data visualization for plots.
from ml4h.TensorMap import TensorMap
from ml4h.visualization_tools.ecg_reshape import DEFAULT_RESTING_ECG_SIGNAL_TMAP
from ml4h.visualization_tools.ecg_reshape import reshape_exercise_ecg_to_tidy
from ml4h.visualization_tools.ecg_reshape import reshape_resting_ecg_to_tidy

# Configure for large-data altair plotting.
# https://altair-viz.github.io/user_guide/faq.html#why-does-altair-lead-to-such-extremely-large-notebooks
alt.data_transformers.disable_max_rows()

# Define this at the global level so that the tempfile exists for the duration of the Jupyter kernel.
RESTING_ECG_DATA_FILE = tempfile.NamedTemporaryFile(
    dir=os.getcwd(),
    prefix='tidy_resting_ecg_',
    suffix='.json',
)
EXERCISE_ECG_TREND_DATA_FILE = tempfile.NamedTemporaryFile(
    dir=os.getcwd(),
    prefix='tidy_exercise_ecg_trend_',
    suffix='.json',
)
EXERCISE_ECG_SIGNAL_DATA_FILE = tempfile.NamedTemporaryFile(
    dir=os.getcwd(),
    prefix='tidy_exercise_ecg_signal_',
    suffix='.json',
)


def resting_ecg_interactive_plot(
    sample_id: Union[int, str], folder: Optional[str] = None,
    tmap: TensorMap = DEFAULT_RESTING_ECG_SIGNAL_TMAP,
) -> Union[HTML, alt.Chart]:
  """Wrangle resting ECG data to tidy and present it as an interactive plot.

  Args:
    sample_id: The id of the ECG sample to retrieve.
    folder: The local or Cloud Storage folder under which the files reside.
    tmap: The TensorMap to use for ECG input.

  Returns:
    An Altair plot or a notebook-friendly error.
  """
  tidy_resting_ecg_signal = reshape_resting_ecg_to_tidy(sample_id, folder, tmap)
  if tidy_resting_ecg_signal.shape[0] == 0:
    return HTML(f'''
      <div class="alert alert-block alert-danger">
      <b>Warning:</b> Resting ECG not available for sample {sample_id}.<br>
      Use the <kbd>folder</kbd> parameter to read HD5s from a different local directory or Cloud Storage bucket.
      </div>''')

  data_file = os.path.basename(RESTING_ECG_DATA_FILE.name)
  tidy_resting_ecg_signal.query("filtering in ['raw_mV']").to_json(data_file, orient='records')

  # Define the plot components.
  brush = alt.selection(type='interval', encodings=['x'])

  lead_dropdown = alt.binding_select(options=list(tidy_resting_ecg_signal.lead.unique()))
  lead_select = alt.selection_single(
      fields=['lead'], bind=lead_dropdown,
      name='Choose just one to view',
      init={'lead': tidy_resting_ecg_signal.lead.unique()[0]},
  )

  base = alt.Chart(data_file).mark_line().encode(
      x='ts_reference:Q',
      y='signal_mV:Q',
      color=alt.Color(
          'lead:N', legend=alt.Legend(orient='top'),
          title='Lead(s) currently displayed',
      ),
  ).properties(
      width=900, height=250, title=f'Resting ECG for {sample_id}',
  ).add_selection(
      lead_select,
  ).transform_filter(lead_select)

  upper = base.encode(x=alt.X('ts_reference:Q', scale=alt.Scale(domain=brush)))

  lower = base.properties(
      height=50, title='Brush over this subplot to select a time interval.',
  ).add_selection(brush)

  return upper & lower


def exercise_ecg_interactive_plot(
    sample_id: Union[int, str], folder: Optional[str] = None, time_interval_seconds: int = 10,
) -> Union[HTML, alt.Chart]:
  """Wrangle exercise ECG data to tidy and present it as an interactive plot.

  Args:
    sample_id: The id of the ECG sample to retrieve.
    folder: The local or Cloud Storage folder under which the files reside.
    time_interval_seconds: the width of the time interval (in seconds) to display of signal data

  Returns:
    An Altair plot or a notebook-friendly error.
  """
  (exercise_ecg_trend, exercise_ecg_signal) = reshape_exercise_ecg_to_tidy(sample_id=sample_id, folder=folder)
  if(exercise_ecg_trend.shape[0] == 0 or exercise_ecg_signal.shape[0] == 0):
    return HTML(f'''
      <div class="alert alert-block alert-danger">
      <b>Warning:</b> Exercise ECG not available for sample {sample_id}.<br>
      Use the <kbd>folder</kbd> parameter to read HD5s from a different local directory or Cloud Storage bucket.
      </div>''')

  trend_data_file = os.path.basename(EXERCISE_ECG_TREND_DATA_FILE.name)
  exercise_ecg_trend.to_json(trend_data_file, orient='records')
  signal_data_file = os.path.basename(EXERCISE_ECG_SIGNAL_DATA_FILE.name)
  exercise_ecg_signal.to_json(signal_data_file, orient='records')

  brush = alt.selection_single(on='mouseover', nearest=True, fields=['time'], init={'time': 200.0})

  lead_dropdown = alt.binding_select(options=list(exercise_ecg_signal.lead.unique()))
  lead_select = alt.selection_single(
      fields=['lead'], bind=lead_dropdown, name='Choose just one to view',
      init={'lead': exercise_ecg_signal.lead.unique()[0]},
  )

  trend = alt.Chart(trend_data_file).mark_point(opacity=0.8, filled=True, size=100).encode(
      x='time:Q',
      color=alt.Color('phasename:N', legend=alt.Legend(orient='top'), title='Phase names'),
      tooltip=[
          'artifact:Q', 'grade:Q', 'heartrate:Q', 'load:Q', 'mets:Q', 'pacecount:Q',
          'phasename:N', 'phasetime:Q', 'time:Q', 'vecount:Q',
      ],
  ).properties(
      width=900, height=100, title=f'Click on a point to select a {time_interval_seconds} second time interval.',
  ).add_selection(brush)

  signal = alt.Chart(signal_data_file).mark_line().encode(
      alt.X('time:Q', axis=alt.Axis(labelAngle=15)),
      y='raw_mV:Q',
      color=alt.Color('lead:N', legend=alt.Legend(orient='top'), title='Lead names'),
  ).properties(
      width=900, height=300, title='Exercise ECG signal for {}'.format(sample_id),
  ).add_selection(
      lead_select,
  ).transform_filter(
      lead_select,
  ).transform_filter(
      # https://github.com/altair-viz/altair/issues/1960
      f'''((toNumber({brush.name}.time) - {time_interval_seconds/2.0}) < datum.time)
           && (datum.time < toNumber({brush.name}.time) + {time_interval_seconds/2.0})''',
  )

  return trend.encode(y='heartrate:Q') & trend.encode(y='load:Q') & signal
