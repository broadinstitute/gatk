import h5py
import numpy as np
from ml4cvd.TensorMap import TensorMap
from typing import List, Dict
from ml4cvd.defines import DataSetType
from ml4cvd.tensor_writer_ukbb import tensor_path, path_date_to_datetime


"""
For now, all we will map `group` in TensorMap to `source` in tensor_path and `name` to `name`
"""


def _all_dates(hd5: h5py.File, source: str, dtype: DataSetType, name: str) -> List[str]:
    """
    Gets the dates in the hd5 with source, dtype, name.
    """
    # TODO: This ideally would be implemented to not depend on the order of name, date, dtype, source in the hd5s
    # Unfortunately, that's hard to do efficiently
    return hd5[source][str(dtype)][name]


def _fail_nan(tensor):
    if np.isnan(tensor).any():
        raise ValueError('Tensor contains nans.')
    return tensor


def _nan_to_mean(tensor, max_allowed_nan_fraction=.2):
    tensor_isnan = np.isnan(tensor)
    if np.count_nonzero(tensor_isnan) / tensor.size > max_allowed_nan_fraction:
        raise ValueError('Tensor contains too many nans.')
    tensor[tensor_isnan] = np.nanmean(tensor)
    return tensor


def _get_tensor_at_first_date(hd5: h5py.File, source: str, dtype: DataSetType, name: str, handle_nan=_fail_nan):
    """
    Gets the numpy array at the first date of source, dtype, name.
    """
    dates = _all_dates(hd5, source, dtype, name)
    if not dates:
        raise ValueError(f'No {name} values values available.')
    # TODO: weird to convert date from string to datetime, because it just gets converted back.
    first_date = path_date_to_datetime(min(dates))  # Date format is sortable. 
    first_date_path = tensor_path(source=source, dtype=dtype, name=name, date=first_date)
    tensor = np.array(hd5[first_date_path])
    tensor = handle_nan(tensor)
    return tensor


def _pad_array_to_shape(tm: TensorMap, original: np.ndarray):
    if tm.shape == original.shape:
        return original
    padded = np.zeros(tm.shape)
    padded[:original.shape[0]] = original.reshape((original.shape[0],) + tm.shape[1:])
    return padded


def normalized_first_date(tm: TensorMap, hd5: h5py.File, dependents=None):
    tensor = _get_tensor_at_first_date(hd5, tm.group, tm.dtype, tm.name)
    if tm.dtype == DataSetType.CONTINUOUS:
        return tm.normalize(tensor)
    if tm.dtype == DataSetType.FLOAT_ARRAY:
        tensor = tm.zero_mean_std1(tensor)
        return _pad_array_to_shape(tm, tensor)
    raise ValueError(f'normalize_first_date not implemented for {tm.dtype}')


TMAPS: Dict[str, TensorMap] = {}


# BIKE ECG
def _check_phase_full_len(hd5: h5py.File, phase: str):
    phase_len = _get_tensor_at_first_date(hd5, 'ecg_bike', DataSetType.CONTINUOUS, f'{phase}_duration')
    valid = True
    if phase == 'pretest':
        valid &= phase_len == 15
    elif phase == 'exercise':
        valid &= phase_len == 360
    elif phase == 'rest':
        valid &= phase_len == 60
    else:
        raise ValueError(f'Phase {phase} is not a valid phase.')
    if not valid:
        raise ValueError(f'{phase} phase is not full length')


def _first_date_bike_recovery(tm: TensorMap, hd5: h5py.File, dependents=None):
    _check_phase_full_len(hd5, 'rest')
    original = _get_tensor_at_first_date(hd5, tm.group, DataSetType.FLOAT_ARRAY, tm.name)
    flat = original.mean(axis=1)  # all leads are basically the same
    recovery = flat[-tm.shape[0]:]
    return tm.zero_mean_std1(recovery).reshape(tm.shape)


def _first_date_bike_pretest(tm: TensorMap, hd5: h5py.File, dependents=None):
    _check_phase_full_len(hd5, 'pretest')
    original = _get_tensor_at_first_date(hd5, tm.group, DataSetType.FLOAT_ARRAY, tm.name)
    flat = original.mean(axis=1)  # all leads are basically the same
    recovery = flat[:tm.shape[0]]
    return tm.normalize(recovery).reshape(tm.shape)


def _first_date_hrr(tm: TensorMap, hd5: h5py.File, dependents=None):
    _check_phase_full_len(hd5, 'rest')
    last_hr = _get_tensor_at_first_date(hd5, 'ecg_bike', DataSetType.FLOAT_ARRAY, 'trend_heartrate')[-1]
    max_hr = _get_tensor_at_first_date(hd5, 'ecg_bike', DataSetType.CONTINUOUS, 'max_hr')
    return tm.normalize(max_hr - last_hr)


def _healthy_check(hd5):
    for phase in ('pretest', 'exercise', 'rest'):
        _check_phase_full_len(hd5, phase)
    max_load = max(_get_tensor_at_first_date(hd5, 'ecg_bike', DataSetType.FLOAT_ARRAY, 'trend_load'))
    if max_load < 60:
        raise ValueError('Max load not high enough')


def _healthy_bike(tm: TensorMap, hd5: h5py.File, dependents=None):
    _healthy_check(hd5)
    return normalized_first_date(tm, hd5)


def _healthy_hrr(tm: TensorMap, hd5: h5py.File, dependents=None):
    _healthy_check(hd5)
    return _first_date_hrr(tm, hd5)


def _median_pretest(tm: TensorMap, hd5: h5py.File, dependents=None):
    _healthy_check(hd5)
    times = _get_tensor_at_first_date(hd5, 'ecg_bike', DataSetType.FLOAT_ARRAY, 'trend_time')
    tensor = np.abs(_get_tensor_at_first_date(hd5, tm.group, DataSetType.FLOAT_ARRAY, tm.name))
    return tm.normalize(np.median(tensor[times <= 15]))


TMAPS['ecg-bike-hrr'] = TensorMap('hrr', group='ecg_bike', loss='logcosh', metrics=['mape'], shape=(1,),
                                  normalization={'mean': 30.55, 'std': 12.81},
                                  tensor_from_file=_first_date_hrr, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-healthy-max-hr'] = TensorMap('max_hr', group='ecg_bike', loss='logcosh', metrics=['mape'],
                                             normalization={'mean': 113.7, 'std': 13.3}, shape=(1,),
                                             tensor_from_file=_healthy_bike, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-healthy-hrr'] = TensorMap('hrr', group='ecg_bike', loss='logcosh', metrics=['mape'], shape=(1,),
                                          normalization={'mean': 30.47, 'std': 11.76},
                                          tensor_from_file=_healthy_hrr, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-healthy-resting'] = TensorMap('resting_hr', group='ecg_bike', loss='logcosh', metrics=['mape'], shape=(1,),
                                              normalization={'mean': 70.0, 'std': 11.62},
                                              tensor_from_file=_healthy_bike, dtype=DataSetType.CONTINUOUS)

TMAPS['ecg-bike-med-pretest-hr'] = TensorMap('trend_heartrate', group='ecg_bike', loss='logcosh', metrics=['mape'], shape=(1,),
                                             normalization={'mean': 70., 'std': 11.},
                                             tensor_from_file=_median_pretest, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-med-pretest-stamp'] = TensorMap('trend_stamplitude', group='ecg_bike', loss='logcosh', metrics=['mape'], shape=(1,),
                                                normalization={'mean': .03, 'std': .03},
                                                tensor_from_file=_median_pretest, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-med-pretest-jpoint'] = TensorMap('trend_jpointamplitude', group='ecg_bike', loss='logcosh', metrics=['mape'], shape=(1,),
                                                 normalization={'mean': .032, 'std': .46},
                                                 tensor_from_file=_median_pretest, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-med-pretest-stamp20'] = TensorMap('trend_stamplitude20ms', group='ecg_bike', loss='logcosh', metrics=['mape'], shape=(1,),
                                                  normalization={'mean': .03, 'std': .03},
                                                  tensor_from_file=_median_pretest, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-recovery'] = TensorMap('full', shape=(30000, 1), group='ecg_bike',
                                       tensor_from_file=_first_date_bike_recovery, dtype=DataSetType.FLOAT_ARRAY)
TMAPS['ecg-bike-pretest'] = TensorMap('full', shape=(500 * 15, 1), group='ecg_bike',
                                      normalization={'mean': 8.3, 'std': 33.},
                                      tensor_from_file=_first_date_bike_pretest, dtype=DataSetType.FLOAT_ARRAY)
