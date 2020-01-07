"""Methods for working with raw ECG signal data."""
import os
import tempfile

from biosppy.signals.tools import filter_signal
import h5py
import numpy as np
import pandas as pd
import tensorflow as tf

DEFAULT_RESTING_ECG_HD5_FOLDERS = {
    'fake': 'gs://ml4cvd/projects/fake_ecgs/',
    'ukb': 'gs://ml4cvd/rest-ecg-hd5s/2019-11-19/'
}
DEFAULT_RAW_EXERCISE_DATA_PATH = 'gs://uk-biobank-sek-data-ttl-one-week/fc-bccc59a1-eb3e-456c-99f7-19593d55c953/ecg-rest-and-bike-with-trend-tensors/'

RAW_SCALE = 0.005  # Convert to mV.
SAMPLING_RATE = 500.0
RESTING_SIGNAL_LENGTH = 5000
EXERCISE_SIGNAL_LENGTH = 30000


def reshape_resting_ecg_to_tidy(sample_id, gcs_folder=None):
  """Wrangle raw resting ECG data to tidy.

  TODO:
    * Refactor this to reduce duplicate code from elsewhere in this repository.

  Args:
    sample_id: The id of the ECG sample to retrieve.
    gcs_folder: The local or Cloud Storage folder under which the files reside.

  Returns:
    A pandas dataframe in tidy format or a notebook-friendly error.
  """
  if gcs_folder is None:
    if 'fake' in str(sample_id):
      gcs_folder = DEFAULT_RESTING_ECG_HD5_FOLDERS['fake']
    else:
      gcs_folder = DEFAULT_RESTING_ECG_HD5_FOLDERS['ukb']

  data = {'lead': [], 'raw': [], 'ts_reference': [],
          'filtered': [], 'filtered_1': [], 'filtered_2': []}

  with tempfile.TemporaryDirectory() as tmpdirname:
    sample_hd5 = str(sample_id) + '.hd5'
    local_path = os.path.join(tmpdirname, sample_hd5)
    try:
      tf.io.gfile.copy(src=os.path.join(gcs_folder, sample_hd5),
                       dst=local_path)
    except (tf.errors.NotFoundError, tf.errors.PermissionDeniedError) as e:
      print('Warning: Resting ECG raw signal not available for sample ',
            sample_id,
            '\n\n',
            e.message)
      return pd.DataFrame(data)

    with h5py.File(local_path, mode='r') as hd5:
      if 'ecg_rest' not in hd5:
        return None

      # Loop over all ecg_rest fields. If they are the raw waveforms (not
      # medians), use biosppy package to apply band-pass filter and store
      # additional data.
      for field in list(hd5['ecg_rest'].keys()):
        signal = hd5['ecg_rest'][field][:]
        signal_length = len(signal)
        if signal_length == RESTING_SIGNAL_LENGTH:  # If 5000 steps long, this is raw.
          data['raw'].extend(signal)
          data['lead'].extend([field] * signal_length)
          data['ts_reference'].extend(np.array(
              [i*1./(SAMPLING_RATE+1.) for i in range(0, signal_length)]))
          filtered, _, _ = filter_signal(signal=signal,
                                         ftype='FIR',
                                         band='bandpass',
                                         order=int(0.3 * SAMPLING_RATE),
                                         frequency=[.9, 50],
                                         sampling_rate=SAMPLING_RATE)
          data['filtered'].extend(filtered)
          filtered_1, _, _ = filter_signal(signal=signal,
                                           ftype='FIR',
                                           band='bandpass',
                                           order=int(0.3 * SAMPLING_RATE),
                                           frequency=[.9, 20],
                                           sampling_rate=SAMPLING_RATE)
          data['filtered_1'].extend(filtered_1)
          filtered_2, _, _ = filter_signal(signal=signal,
                                           ftype='FIR',
                                           band='bandpass',
                                           order=int(0.3 * SAMPLING_RATE),
                                           frequency=[.9, 30],
                                           sampling_rate=SAMPLING_RATE)
          data['filtered_2'].extend(filtered_2)

  signal_df = pd.DataFrame(data)
  # Convert the raw signal to mV.
  signal_df['raw_mV'] = signal_df['raw'] * RAW_SCALE
  signal_df['filtered_mV'] = signal_df['filtered'] * RAW_SCALE
  signal_df['filtered_1_mV'] = signal_df['filtered_1'] * RAW_SCALE
  signal_df['filtered_2_mV'] = signal_df['filtered_2'] * RAW_SCALE
  # Reshape to tidy (long format).
  tidy_signal_df = signal_df.melt(id_vars=['lead', 'ts_reference'],
                                  value_vars=['raw_mV', 'filtered_mV',
                                              'filtered_1_mV', 'filtered_2_mV'],
                                  var_name='filtering', value_name='signal_mV')

  # The leads have a meaningful order, apply the order to this column.
  lead_factor_type = pd.api.types.CategoricalDtype(
      categories=['strip_I', 'strip_aVR', 'strip_V1', 'strip_V4',
                  'strip_II', 'strip_aVL', 'strip_V2', 'strip_V5',
                  'strip_III', 'strip_aVF', 'strip_V3', 'strip_V6'],
      ordered=True)
  tidy_signal_df['lead'] = tidy_signal_df.lead.astype(lead_factor_type)

  return tidy_signal_df


def reshape_exercise_ecg_to_tidy(sample_id,
                                 gcs_folder=DEFAULT_RAW_EXERCISE_DATA_PATH):
  """Wrangle raw exercise ECG data to tidy.

  TODO:
    * Per Puneet, this code is not correctly reading the data from the HD5 file.
    * Refactor this to reduce duplicate code from above and elsewhere in this
      repository.

  Args:
    sample_id: The id of the ECG sample to retrieve.
    gcs_folder: The local or Cloud Storage folder under which the files reside.

  Returns:
    A pandas dataframe in tidy format or a notebook-friendly error.
  """
  data = {}
  full_data = {}
  signal_data = {}

  with tempfile.TemporaryDirectory() as tmpdirname:
    sample_hd5 = str(sample_id) + '.hd5'
    local_path = os.path.join(tmpdirname, sample_hd5)
    try:
      tf.io.gfile.copy(src=os.path.join(gcs_folder, sample_hd5),
                       dst=local_path)
    except (tf.errors.NotFoundError, tf.errors.PermissionDeniedError) as e:
      print('Error: Exercise ECG raw signal not available for sample ',
            sample_id,
            '\n\n',
            e.message)
      return (pd.DataFrame(data), pd.DataFrame(data))

    with h5py.File(local_path, mode='r') as hd5:
      if 'ecg_bike_recovery' not in hd5:
        print('Warning: Exercise ECG raw signal does not contain ',
              'ecg_bike_recovery for sample ',
              sample_id,
              '\n\n',
              e.message)
        return (pd.DataFrame(data), pd.DataFrame(data))

      trends = hd5['ecg_bike_trend']
      rest_start_idx = list(trends['PhaseName']).index(2)

      # shape (57, )
      full_data['full_times'] = np.array(list(trends['time']))
      full_data['full_hrs'] = np.array(list(trends['HeartRate']))
      full_data['full_loads'] = np.array(list(trends['Load']))
      full_data['artifacts'] = np.array(list(trends['Artifact']))

      # shape (6, )
      data['times'] = np.array(list(trends['PhaseTime'])[rest_start_idx:])
      data['hrs'] = np.array(list(trends['HeartRate'])[rest_start_idx:])
      data['loads'] = np.array(list(trends['Load'])[rest_start_idx:])

      # shape (30000,)
      signal_data['ts_reference'] = np.array(
          [i * 1. / (SAMPLING_RATE + 1.) for i in range(0, EXERCISE_SIGNAL_LENGTH)])
      for lead in list(hd5['ecg_bike_recovery'].keys()):
        raw = np.array(hd5['ecg_bike_recovery'][lead])
        filtered, _, _ = filter_signal(signal=np.array(raw),
                                       ftype='FIR',
                                       band='bandpass',
                                       order=int(0.3 * SAMPLING_RATE),
                                       frequency=[.9, 50],
                                       sampling_rate=SAMPLING_RATE)
        signal_data['raw_mV_' + lead] = raw * RAW_SCALE
        signal_data['filtered_mV_' + lead] = filtered * RAW_SCALE
      signal_data['raw_mV_average'] = (
          1/3.0 * (signal_data['raw_mV_lead_I']
                   + signal_data['raw_mV_lead_2']
                   + signal_data['raw_mV_lead_3']))
      signal_data['filtered_mV_average'] = (
          1/3.0 * (signal_data['filtered_mV_lead_I']
                   + signal_data['filtered_mV_lead_2']
                   + signal_data['filtered_mV_lead_3']))

      # shape (1,)
      data['max_hr'] = list(hd5['continuous']['bike_max_hr'])[0]
      data['resting_hr'] = list(hd5['continuous']['bike_resting_hr'])[0]
      data['max_pred_hr'] = list(hd5['continuous']['bike_max_pred_hr'])[0]

  # shape (1,)
  data['hr_recovery_60s'] = data['max_hr'] - data['hrs'][-1]
  data['hr_increase'] = data['max_hr'] - data['resting_hr']
  data['max_hr_achieved_ratio'] = data['max_hr'] / data['max_pred_hr']
  data['pred_peak_exercise'] = data['max_pred_hr'] - data['resting_hr']
  data['peak_hr_reserve'] = data['max_hr'] / data['pred_peak_exercise']

  full_df = pd.DataFrame(full_data)
  signal_df = pd.DataFrame(signal_data)
  tidy_signal_df = pd.wide_to_long(signal_df,
                                   stubnames=['raw_mV', 'filtered_mV'],
                                   i='ts_reference',
                                   j='lead',
                                   sep='_',
                                   suffix='.*')
  tidy_signal_df.reset_index(inplace=True)  # Turn pd multiindex into columns.

  # The leads have a meaningful order, apply the order to this column.
  lead_factor_type = pd.api.types.CategoricalDtype(
      categories=['average', 'lead_I', 'lead_2', 'lead_3'],
      ordered=True)
  tidy_signal_df['lead'] = tidy_signal_df.lead.astype(lead_factor_type)

  return (full_df, tidy_signal_df)

