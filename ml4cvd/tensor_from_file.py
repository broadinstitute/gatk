import datetime
from typing import List, Dict, Tuple, Callable

import os
import csv
import h5py
import logging
import numpy as np
from keras.utils import to_categorical

from ml4cvd.metrics import weighted_crossentropy
from ml4cvd.tensor_writer_ukbb import tensor_path, path_date_to_datetime
from ml4cvd.TensorMap import TensorMap, no_nans, str2date, make_range_validator
from ml4cvd.defines import ECG_REST_LEADS, ECG_REST_MEDIAN_LEADS, ECG_REST_AMP_LEADS
from ml4cvd.defines import DataSetType, EPS, MRI_TO_SEGMENT, MRI_SEGMENTED, MRI_SEGMENTED_CHANNEL_MAP


"""
For now, all we will map `group` in TensorMap to `source` in tensor_path and `name` to `name`
"""


def normalized_first_date(tm: TensorMap, hd5: h5py.File, dependents=None):
    tensor = _get_tensor_at_first_date(hd5, tm.group, tm.dtype, tm.name)
    if tm.dtype == DataSetType.CONTINUOUS:
        return tm.normalize_and_validate(tensor)
    if tm.dtype == DataSetType.FLOAT_ARRAY:
        tensor = tm.normalize_and_validate(tensor)
        return _pad_or_crop_array_to_shape(tm.shape, tensor)
    if tm.dtype == DataSetType.CATEGORICAL:
        tensor = tm.normalize_and_validate(tensor)
        return _pad_or_crop_array_to_shape(tm.shape, tensor)
    raise ValueError(f'normalize_first_date not implemented for {tm.dtype}')


def _random_slice_tensor(tensor_key, dependent_key=None):
    def _random_slice_tensor_from_file(tm: TensorMap, hd5: h5py.File, dependents=None):
        big_tensor = _get_tensor_at_first_date(hd5, tm.group, tm.dtype, tensor_key)
        cur_slice = np.random.choice(range(big_tensor.shape[-1]))
        tensor = np.zeros(tm.shape, dtype=np.float32)
        tensor[..., 0] = big_tensor[..., cur_slice]
        if dependent_key is not None:
            dependents[tm.dependent_map] = np.zeros(tm.dependent_map.shape, dtype=np.float32)
            label_tensor = np.array(hd5[dependent_key][..., cur_slice], dtype=np.float32)
            dependents[tm.dependent_map][:, :, :] = to_categorical(label_tensor, tm.dependent_map.shape[-1])
        return tm.normalize_and_validate(tensor)
    return _random_slice_tensor_from_file


def _slice_subset_tensor(tensor_key, start, stop, step=1, dependent_key=None, pad_shape=None, dtype_override=None, allow_channels=True):
    def _slice_subset_tensor_from_file(tm: TensorMap, hd5: h5py.File, dependents=None):
        if dtype_override is not None:
            big_tensor = _get_tensor_at_first_date(hd5, tm.group, dtype_override, tensor_key)
        else:
            big_tensor = _get_tensor_at_first_date(hd5, tm.group, tm.dtype, tensor_key)

        if pad_shape is not None:
            big_tensor = _pad_or_crop_array_to_shape(pad_shape, big_tensor)

        if allow_channels and tm.shape[-1] < (stop-start) // step:
            tensor = big_tensor[..., np.arange(start, stop, step), :]
        else:
            tensor = big_tensor[..., np.arange(start, stop, step)]

        if dependent_key is not None:
            label_tensor = np.array(hd5[dependent_key][..., start:stop], dtype=np.float32)
            dependents[tm.dependent_map] = to_categorical(label_tensor, tm.dependent_map.shape[-1])
        return tm.normalize_and_validate(tensor)
    return _slice_subset_tensor_from_file


def _build_tensor_from_file(file_name: str, target_column: str, normalization: bool = False, delimiter: str = '\t'):
    """
    Build a tensor_from_file function from a column in a file.
    Only works for continuous values.
    When normalization is True values will be normalized according to the mean and std of all of the values in the column.
    """
    error = None
    try:
        with open(file_name, 'r') as f:
            reader = csv.reader(f, delimiter=delimiter)
            header = next(reader)
            index = header.index(target_column)
            table = {row[0]: np.array([float(row[index])]) for row in reader}
            if normalization:
                value_array = np.array([sub_array[0] for sub_array in table.values()])
                mean = value_array.mean()
                std = value_array.std()
                logging.info(f'Normalizing TensorMap from file {file_name}, column {target_column} with mean: '
                             f'{mean:.2f}, std: {std:.2f}')
    except FileNotFoundError as e:
        error = e

    def tensor_from_file(tm: TensorMap, hd5: h5py.File, dependents=None):
        if error:
            raise error
        try:
            t = table[os.path.basename(hd5.filename).replace('.hd5', '')]
            if normalization:
                tm.normalization = {'mean': mean, 'std': std}
            tn = tm.normalize_and_validate(t)
            return tn
        except KeyError:
            raise KeyError(f'User id not in file {file_name}.')
    return tensor_from_file


def _survival_tensor(start_date_key, day_window):
    def _survival_tensor_from_file(tm: TensorMap, hd5: h5py.File, dependents=None):
        assess_date = str2date(str(hd5[start_date_key][0]))
        has_disease = 0   # Assume no disease if the tensor does not have the dataset
        if tm.name in hd5['categorical']:
            has_disease = int(hd5['categorical'][tm.name][0])

        if tm.name + '_date' in hd5['dates']:
            censor_date = str2date(str(hd5['dates'][tm.name + '_date'][0]))
        elif 'phenotype_censor' in hd5['dates']:
            censor_date = str2date(str(hd5['dates/phenotype_censor']))
        else:
            raise ValueError(f'No date found for survival {tm.name}')

        intervals = int(tm.shape[0] / 2)
        days_per_interval = day_window / intervals
        survival_then_censor = np.zeros(tm.shape, dtype=np.float32)
        for i, day_delta in enumerate(np.arange(0, day_window, days_per_interval)):
            cur_date = assess_date + datetime.timedelta(days=day_delta)
            survival_then_censor[i] = float(cur_date < censor_date)
            survival_then_censor[intervals+i] = has_disease * float(censor_date <= cur_date < censor_date + datetime.timedelta(days=days_per_interval))
            if i == 0 and censor_date <= cur_date:  # Handle prevalent diseases
                survival_then_censor[intervals] = has_disease
        return survival_then_censor

    return _survival_tensor_from_file


def _age_in_years_tensor(date_key, birth_key='continuous/34_Year-of-birth_0_0'):
    def age_at_tensor_from_file(tm: TensorMap, hd5: h5py.File, dependents=None):
        assess_date = str2date(str(hd5[date_key][0]))
        birth_year = hd5[birth_key][0]
        return tm.normalize_and_validate(np.array([assess_date.year-birth_year]))
    return age_at_tensor_from_file


def _all_dates(hd5: h5py.File, source: str, dtype: DataSetType, name: str) -> List[str]:
    """
    Gets the dates in the hd5 with source, dtype, name.
    """
    # TODO: This ideally would be implemented to not depend on the order of name, date, dtype, source in the hd5s
    # Unfortunately, that's hard to do efficiently
    return hd5[source][str(dtype)][name]


def _pass_nan(tensor):
    return tensor


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
    tensor = np.array(hd5[first_date_path], dtype=np.float32)
    tensor = handle_nan(tensor)
    return tensor


def _pad_or_crop_array_to_shape(new_shape: Tuple, original: np.ndarray):
    if new_shape == original.shape:
        return original
    
    result = np.zeros(new_shape)
    slices = tuple(slice(min(original.shape[i], new_shape[i])) for i in range(len(original.shape)))

    # Allow expanding one dimension eg (256, 256) can become (256, 256, 1)
    if len(new_shape) - len(original.shape) == 1:
        padded = result[..., 0]
    else:
        padded = result

    padded[slices] = original[slices]
    return result


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
    recovery = original[-tm.shape[0]:]
    return tm.normalize_and_validate(recovery).reshape(tm.shape)


def _first_date_bike_pretest(tm: TensorMap, hd5: h5py.File, dependents=None):
    _check_phase_full_len(hd5, 'pretest')
    original = _get_tensor_at_first_date(hd5, tm.group, DataSetType.FLOAT_ARRAY, tm.name)
    pretest = original[:tm.shape[0]]
    return tm.normalize_and_validate(pretest).reshape(tm.shape)


def _first_date_hrr(tm: TensorMap, hd5: h5py.File, dependents=None):
    _check_phase_full_len(hd5, 'rest')
    last_hr = _get_tensor_at_first_date(hd5, 'ecg_bike', DataSetType.FLOAT_ARRAY, 'trend_heartrate')[-1]
    max_hr = _get_tensor_at_first_date(hd5, 'ecg_bike', DataSetType.CONTINUOUS, 'max_hr')
    return tm.normalize_and_validate(max_hr - last_hr)


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
    return tm.normalize_and_validate(np.median(tensor[times <= 15]))


def _new_hrr(tm: TensorMap, hd5: h5py.File, dependents=None):
    _check_phase_full_len(hd5, 'rest')
    hrs = _get_tensor_at_first_date(hd5, 'ecg_bike', DataSetType.FLOAT_ARRAY, 'trend_heartrate')
    phases = _get_tensor_at_first_date(hd5, 'ecg_bike', DataSetType.FLOAT_ARRAY, 'trend_phasename')
    min_hr = hrs[phases == 2].min()
    max_hr = _get_tensor_at_first_date(hd5, 'ecg_bike', DataSetType.CONTINUOUS, 'max_hr')
    max_pred = _get_tensor_at_first_date(hd5, 'ecg_bike', DataSetType.CONTINUOUS, 'max_pred_hr')
    hrr = max_hr - min_hr
    if max_hr / max_pred > 150:
        raise ValueError('Max hr / max pred hr too high.')
    if hrr > 80:
        raise ValueError('HRR too high.')
    return tm.normalize_and_validate(hrr)


_HRR_SENTINEL = -1000


def _sentinel_hrr(tm: TensorMap, hd5: h5py.File, dependents=None):
    try:
        return _new_hrr(tm, hd5)
    except ValueError:
        return _HRR_SENTINEL


def _hr_achieved(tm: TensorMap, hd5: h5py.File, dependents=None):
    _check_phase_full_len(hd5, 'rest')
    max_hr = _get_tensor_at_first_date(hd5, 'ecg_bike', DataSetType.CONTINUOUS, 'max_hr')
    max_pred = _get_tensor_at_first_date(hd5, 'ecg_bike', DataSetType.CONTINUOUS, 'max_pred_hr')
    return tm.normalize_and_validate(max_hr / max_pred)


TMAPS: Dict[str, TensorMap] = dict()


TMAPS['ecg-bike-hrr'] = TensorMap('hrr', group='ecg_bike', loss='logcosh', metrics=['mae'], shape=(1,),
                                  normalization={'mean': 30.55, 'std': 12.81},
                                  tensor_from_file=_first_date_hrr, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-healthy-max-hr'] = TensorMap('max_hr', group='ecg_bike', loss='logcosh', metrics=['mae'],
                                             normalization={'mean': 113.7, 'std': 13.3}, shape=(1,),
                                             tensor_from_file=_healthy_bike, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-healthy-hrr'] = TensorMap('hrr', group='ecg_bike', loss='logcosh', metrics=['mae'], shape=(1,),
                                          normalization={'mean': 30.47, 'std': 11.76},
                                          tensor_from_file=_healthy_hrr, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-healthy-resting'] = TensorMap('resting_hr', group='ecg_bike', loss='logcosh', metrics=['mae'], shape=(1,),
                                              normalization={'mean': 70.0, 'std': 11.62},
                                              tensor_from_file=_healthy_bike, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-med-pretest-hr'] = TensorMap('trend_heartrate', group='ecg_bike', loss='logcosh', metrics=['mae'], shape=(1,),
                                             normalization={'mean': 70., 'std': 11.},
                                             tensor_from_file=_median_pretest, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-med-pretest-stamp'] = TensorMap('trend_stamplitude', group='ecg_bike', loss='logcosh', metrics=['mae'], shape=(1,),
                                                normalization={'mean': .03, 'std': .03},
                                                tensor_from_file=_median_pretest, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-med-pretest-jpoint'] = TensorMap('trend_jpointamplitude', group='ecg_bike', loss='logcosh', metrics=['mae'], shape=(1,),
                                                 normalization={'mean': .032, 'std': .46},
                                                 tensor_from_file=_median_pretest, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-med-pretest-stamp20'] = TensorMap('trend_stamplitude20ms', group='ecg_bike', loss='logcosh', metrics=['mae'], shape=(1,),
                                                  normalization={'mean': .03, 'std': .03},
                                                  tensor_from_file=_median_pretest, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-recovery'] = TensorMap('full', shape=(30000, 1), group='ecg_bike', validator=no_nans,
                                       tensor_from_file=_first_date_bike_recovery, dtype=DataSetType.FLOAT_ARRAY)
TMAPS['ecg-bike-pretest'] = TensorMap('full', shape=(500 * 15 - 4, 3), group='ecg_bike', validator=no_nans,
                                      normalization={'mean': np.array([7, -7, 3.5])[np.newaxis], 'std': np.array([31, 30, 16])[np.newaxis]},
                                      tensor_from_file=_first_date_bike_pretest, dtype=DataSetType.FLOAT_ARRAY)
TMAPS['ecg-bike-pretest-5k'] = TensorMap('full', shape=(5000, 3), group='ecg_bike', validator=no_nans,
                                      normalization={'mean': np.array([7, -7, 3.5])[np.newaxis], 'std': np.array([31, 30, 16])[np.newaxis]},
                                      tensor_from_file=_first_date_bike_pretest, dtype=DataSetType.FLOAT_ARRAY)
TMAPS['ecg-bike-new-hrr'] = TensorMap('hrr', group='ecg_bike', loss='logcosh', metrics=['mae'], shape=(1,),
                                      normalization={'mean': 31, 'std': 12},
                                      tensor_from_file=_new_hrr, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-hrr-sentinel'] = TensorMap('hrr', group='ecg_bike', metrics=['mae'], shape=(1,),
                                           normalization={'mean': 31, 'std': 12}, sentinel=_HRR_SENTINEL,
                                           tensor_from_file=_sentinel_hrr, dtype=DataSetType.CONTINUOUS)
TMAPS['ecg-bike-hrr-student'] = TensorMap('hrr', group='ecg_bike', metrics=['mae'], shape=(1,),
                                          normalization={'mean': 31, 'std': 12}, sentinel=_HRR_SENTINEL, dtype=DataSetType.CONTINUOUS,
                                          tensor_from_file=_build_tensor_from_file('inference.tsv', 'ecg-bike-hrr-sentinel_prediction'))
TMAPS['ecg-bike-hr-achieved'] = TensorMap('hr_achieved', group='ecg_bike', loss='logcosh', metrics=['mae'], shape=(1,),
                                          normalization={'mean': .68, 'std': .1},
                                          tensor_from_file=_hr_achieved, dtype=DataSetType.CONTINUOUS)

TMAPS['ecg_rest_afib_hazard'] = TensorMap('atrial_fibrillation_or_flutter', group='proportional_hazard', shape=(100,),
                                          tensor_from_file=_survival_tensor('ecg_rest_date', 365 * 5), dtype=DataSetType.SERIES)
TMAPS['ecg_rest_cad_hazard'] = TensorMap('coronary_artery_disease', group='proportional_hazard', shape=(100,),
                                         tensor_from_file=_survival_tensor('ecg_rest_date', 365 * 5), dtype=DataSetType.SERIES)
TMAPS['ecg_rest_hyp_hazard'] = TensorMap('hypertension', group='proportional_hazard', shape=(100,),
                                         tensor_from_file=_survival_tensor('ecg_rest_date', 365 * 5), dtype=DataSetType.SERIES)
TMAPS['ecg_rest_cad_hazard'] = TensorMap('coronary_artery_disease', group='proportional_hazard', shape=(100,),
                                         tensor_from_file=_survival_tensor('ecg_rest_date', 365 * 5), dtype=DataSetType.SERIES)
TMAPS['enroll_cad_hazard'] = TensorMap('coronary_artery_disease', group='proportional_hazard', shape=(100,),
                                       tensor_from_file=_survival_tensor('dates/enroll_date', 365 * 10), dtype=DataSetType.SERIES)
TMAPS['enroll_hyp_hazard'] = TensorMap('hypertension', group='proportional_hazard', shape=(100,),
                                       tensor_from_file=_survival_tensor('dates/enroll_date', 365 * 10), dtype=DataSetType.SERIES)
TMAPS['enroll_afib_hazard'] = TensorMap('atrial_fibrillation_or_flutter', group='proportional_hazard', shape=(100,),
                                        tensor_from_file=_survival_tensor('dates/enroll_date', 365 * 10), dtype=DataSetType.SERIES)
TMAPS['enroll_chol_hazard'] = TensorMap('hypercholesterolemia', group='proportional_hazard', shape=(100,),
                                        tensor_from_file=_survival_tensor('dates/enroll_date', 365 * 10), dtype=DataSetType.SERIES)
TMAPS['enroll_diabetes2_hazard'] = TensorMap('diabetes_type_2', group='proportional_hazard', shape=(100,),
                                             tensor_from_file=_survival_tensor('dates/enroll_date', 365 * 10), dtype=DataSetType.SERIES)


def _make_ecg_rest(population_normalize: float = None):
    def ecg_rest_from_file(tm, hd5, dependents={}):
        tensor = np.zeros(tm.shape, dtype=np.float32)
        if tm.dependent_map is not None:
            dependents[tm.dependent_map] = np.zeros(tm.dependent_map.shape, dtype=np.float32)
            key_choices = [k for k in hd5[tm.group] if tm.name in k]
            lead_idx = np.random.choice(key_choices)
            tensor = np.reshape(hd5[tm.group][lead_idx][: tensor.shape[0] * tensor.shape[1]], tensor.shape, order='F')
            dependents[tm.dependent_map][:, 0] = np.array(hd5[tm.group][lead_idx.replace(tm.name, tm.dependent_map.name)])
            dependents[tm.dependent_map] = tm.zero_mean_std1(dependents[tm.dependent_map])
        else:
            for k in hd5[tm.group]:
                if k in tm.channel_map:
                    if len(tensor.shape) == 3:  # Grab the stacked tensor maps
                        window_size = tensor.shape[0]
                        channels = tensor.shape[2]
                        new_shape = (window_size, channels)
                        new_total = window_size * channels
                        tensor[:, tm.channel_map[k], :] = np.reshape(hd5[tm.group][k][:new_total], new_shape, order='F')
                    elif tm.name == 'ecg_rest_fft':
                        tensor[:, tm.channel_map[k]] = np.log(np.abs(np.fft.fft(hd5[tm.group][k])) + EPS)
                    else:
                        tensor[:, tm.channel_map[k]] = hd5[tm.group][k]
        if population_normalize is None:
            tensor = tm.zero_mean_std1(tensor)
        else:
            tensor /= population_normalize
        return tensor
    return ecg_rest_from_file


TMAPS['ecg_rest_raw'] = TensorMap('ecg_rest_raw', shape=(5000, 12), group='ecg_rest', tensor_from_file=_make_ecg_rest(population_normalize=2000.0),
                                  channel_map=ECG_REST_LEADS)
TMAPS['ecg_rest_raw_100'] = TensorMap('ecg_rest_raw_100', shape=(5000, 12), group='ecg_rest', tensor_from_file=_make_ecg_rest(population_normalize=100.0),
                                      channel_map=ECG_REST_LEADS)

TMAPS['ecg_rest'] = TensorMap('strip', shape=(5000, 12), group='ecg_rest', tensor_from_file=_make_ecg_rest(),
                              channel_map=ECG_REST_LEADS)


TMAPS['ecg_rest_fft'] = TensorMap('ecg_rest_fft', shape=(5000, 12), group='ecg_rest', tensor_from_file=_make_ecg_rest(),
                                  channel_map=ECG_REST_LEADS)

TMAPS['ecg_rest_stack'] = TensorMap('strip', shape=(600, 12, 8), group='ecg_rest', tensor_from_file=_make_ecg_rest(),
                                    channel_map=ECG_REST_LEADS)

TMAPS['ecg_rest_median_raw'] = TensorMap('median', group='ecg_rest', shape=(600, 12), loss='logcosh', activation='linear', tensor_from_file=_make_ecg_rest(population_normalize=2000.0),
                                     metrics=['mse', 'mae', 'logcosh'], channel_map=ECG_REST_MEDIAN_LEADS)

TMAPS['ecg_rest_median'] = TensorMap('median', group='ecg_rest', shape=(600, 12), loss='logcosh', activation='linear', tensor_from_file=_make_ecg_rest(),
                                     metrics=['mse', 'mae', 'logcosh'], channel_map=ECG_REST_MEDIAN_LEADS)

TMAPS['ecg_rest_median_stack'] = TensorMap('median', group='ecg_rest', shape=(600, 12, 1), activation='linear', tensor_from_file=_make_ecg_rest(),
                                           metrics=['mse', 'mae', 'logcosh'], loss='logcosh', loss_weight=1.0,
                                           channel_map=ECG_REST_MEDIAN_LEADS)

TMAPS['ecg_median_1lead'] = TensorMap('median', group='ecg_rest', shape=(600, 1), loss='logcosh', loss_weight=10.0, tensor_from_file=_make_ecg_rest(),
                                      activation='linear', metrics=['mse', 'mae', 'logcosh'], channel_map={'lead': 0})

TMAPS['ecg_rest_1lead'] = TensorMap('strip', shape=(600, 8), group='ecg_rest', channel_map={'lead': 0}, tensor_from_file=_make_ecg_rest(),
                                    dependent_map=TMAPS['ecg_median_1lead'])


def _get_lead_cm(length):
    lead_cm = {}
    lead_weights = []
    for i in range(length):
        wave_val = i - (length//2)
        lead_cm['w'+str(wave_val).replace('-', '_')] = i
        lead_weights.append((np.abs(wave_val+1)/(length/2)) + 1.0)
    return lead_cm, lead_weights


TMAPS['ecg_median_1lead_categorical'] = TensorMap('median', group='categorical', shape=(600, 32), activation='softmax', tensor_from_file=_make_ecg_rest(),
                                                  channel_map=_get_lead_cm(32)[0],
                                                  loss=weighted_crossentropy(_get_lead_cm(32)[1], 'ecg_median_categorical'))

TMAPS['ecg_rest_1lead_categorical'] = TensorMap('strip', shape=(600, 8), group='ecg_rest', tensor_from_file=_make_ecg_rest(),
                                                channel_map={'window0': 0, 'window1': 1, 'window2': 2, 'window3': 3,
                                                             'window4': 4, 'window5': 5, 'window6': 6, 'window7': 7},
                                                dependent_map=TMAPS['ecg_median_1lead_categorical'])


def _make_rhythm_tensor(skip_poor=True):
    def rhythm_tensor_from_file(tm, hd5, dependents={}):
        categorical_data = np.zeros(tm.shape, dtype=np.float32)
        if skip_poor and 'poor_data_quality' in hd5['categorical']:
            raise ValueError(f'Poor data quality skipped by {tm.name}.')
        ecg_interpretation = str(hd5['ecg_rest_text'][0])
        for channel in tm.channel_map:
            if channel in hd5['categorical']:
                categorical_data[tm.channel_map[channel]] = 1.0
                return categorical_data
        for afib in ['Atrial fibrillation']:
            if afib in ecg_interpretation:
                categorical_data[tm.channel_map['Atrial_fibrillation']] = 1.0
                return categorical_data
        for rhythm in ['sinus', 'Sinus']:
            if rhythm in ecg_interpretation:
                categorical_data[tm.channel_map['Other_sinus_rhythm']] = 1.0
                return categorical_data
        categorical_data[tm.channel_map['Other_rhythm']] = 1.0
        return categorical_data
    return rhythm_tensor_from_file


TMAPS['ecg_rhythm'] = TensorMap('ecg_rhythm', group='categorical', tensor_from_file=_make_rhythm_tensor(),
                                loss=weighted_crossentropy([1.0, 2.0, 3.0, 3.0, 20.0, 20.0], 'ecg_rhythm'),
                                channel_map={'Normal_sinus_rhythm': 0, 'Sinus_bradycardia': 1, 'Marked_sinus_bradycardia': 2, 'Other_sinus_rhythm': 3, 'Atrial_fibrillation': 4, 'Other_rhythm': 5})
TMAPS['ecg_rhythm_poor'] = TensorMap('ecg_rhythm', group='categorical', tensor_from_file=_make_rhythm_tensor(False),
                                loss=weighted_crossentropy([1.0, 2.0, 3.0, 3.0, 20.0, 20.0], 'ecg_rhythm'),
                                channel_map={'Normal_sinus_rhythm': 0, 'Sinus_bradycardia': 1, 'Marked_sinus_bradycardia': 2, 'Other_sinus_rhythm': 3, 'Atrial_fibrillation': 4, 'Other_rhythm': 5})

TMAPS['ecg_rest_age'] = TensorMap('ecg_rest_age', group='continuous', tensor_from_file=_age_in_years_tensor('ecg_rest_date'), loss='logcosh',
                                  channel_map={'ecg_rest_age': 0}, validator=make_range_validator(0, 110), normalization = {'mean': 65, 'std': 7.7})


# Extract RAmplitude and SAmplitude for LVH criteria
def _make_ukb_ecg_rest(population_normalize: float = None):
    def ukb_ecg_rest_from_file(tm, hd5, dependents={}):
        if 'ukb_ecg_rest' not in hd5:
            raise ValueError('Group with R and S amplitudes not present in hd5')
        tensor = _get_tensor_at_first_date(hd5, tm.group, DataSetType.FLOAT_ARRAY, tm.name, _pass_nan)        
        try:            
            if population_normalize is None:
                tensor = tm.zero_mean_std1(tensor)
            else:
                tensor /= population_normalize
        except:
            ValueError(f'Cannot normalize {tm.name}')
        return tensor
    return ukb_ecg_rest_from_file


TMAPS['ecg_rest_ramplitude_raw'] = TensorMap('ramplitude', group='ukb_ecg_rest', shape=(12,), tensor_from_file=_make_ukb_ecg_rest(1.0),
                            loss='logcosh', metrics=['mse', 'mape', 'mae'], loss_weight=1.0)

TMAPS['ecg_rest_samplitude_raw'] = TensorMap('samplitude', group='ukb_ecg_rest', shape=(12,), tensor_from_file=_make_ukb_ecg_rest(1.0),
                            loss='logcosh', metrics=['mse', 'mape', 'mae'], loss_weight=1.0)

TMAPS['ecg_rest_ramplitude'] = TensorMap('ramplitude', group='ukb_ecg_rest', shape=(12,), tensor_from_file=_make_ukb_ecg_rest(),
                            loss='logcosh', metrics=['mse', 'mape', 'mae'], loss_weight=1.0)

TMAPS['ecg_rest_samplitude'] = TensorMap('samplitude', group='ukb_ecg_rest', shape=(12,), tensor_from_file=_make_ukb_ecg_rest(),
                            loss='logcosh', metrics=['mse', 'mape', 'mae'], loss_weight=1.0)


def _make_ukb_ecg_rest_lvh():
    def ukb_ecg_rest_lvh_from_file(tm, hd5, dependents={}):
        # Lead order seems constant and standard throughout, but we could eventually tensorize it from XML
        lead_order = ECG_REST_AMP_LEADS
        avl_min = 1100.0
        sl_min = 3500.0
        cornell_female_min = 2000.0
        cornell_male_min = 2800.0
        if 'ukb_ecg_rest' not in hd5:
            raise ValueError('Group with R and S amplitudes not present in hd5')        
        tensor_ramp = _get_tensor_at_first_date(hd5, tm.group, DataSetType.FLOAT_ARRAY, 'ramplitude', _pass_nan)
        tensor_samp = _get_tensor_at_first_date(hd5, tm.group, DataSetType.FLOAT_ARRAY, 'samplitude', _pass_nan)
        criteria_sleads = [lead_order[l] for l in ['V1', 'V3']]
        criteria_rleads = [lead_order[l] for l in ['aVL', 'V5', 'V6']]
        if np.any(np.isnan(np.union1d(tensor_ramp[criteria_rleads], tensor_samp[criteria_sleads]))):
            raise ValueError('Missing some of the R and S amplitude readings needed to evaluate LVH criteria')        
        is_female = 'Genetic-sex_Female_0_0' in hd5['categorical']
        is_male   = 'Genetic-sex_Male_0_0' in hd5['categorical']
        # If genetic sex not available, try phenotypic
        if not(is_female or is_male):
            is_female = 'Sex_Female_0_0' in hd5['categorical']
            is_male   = 'Sex_Male_0_0' in hd5['categorical']
        # If neither available, raise error
        if not(is_female or is_male):
            raise ValueError('Sex info required to evaluate LVH criteria')        
        if tm.name == 'avl_lvh':
            is_lvh = tensor_ramp[lead_order['aVL']] > avl_min
        elif tm.name == 'sokolow_lyon_lvh':
            is_lvh = tensor_samp[lead_order['V1']] +\
                     np.maximum(tensor_ramp[lead_order['V5']], tensor_ramp[lead_order['V6']]) > sl_min
        elif tm.name == 'cornell_lvh':            
            is_lvh = tensor_ramp[lead_order['aVL']] + tensor_samp[lead_order['V3']]     
            if is_female:
                is_lvh = is_lvh > cornell_female_min
            if is_male:
                is_lvh = is_lvh > cornell_male_min
        else:
            raise ValueError(f'{tm.name} criterion for LVH is not accounted for')
        # Following convention from categorical TMAPS, positive has cmap index 1
        tensor = np.zeros(tm.shape, dtype=np.float32)
        index = 0    
        if is_lvh:
            index = 1
        tensor[index] = 1.0
        return tensor
    return ukb_ecg_rest_lvh_from_file
        

TMAPS['ecg_rest_lvh_avl'] = TensorMap('avl_lvh', group='ukb_ecg_rest', tensor_from_file=_make_ukb_ecg_rest_lvh(),
                            channel_map={'no_avl_lvh': 0, 'aVL LVH': 1},
                            loss=weighted_crossentropy([0.006, 1.0], 'avl_lvh'))

TMAPS['ecg_rest_lvh_sokolow_lyon'] = TensorMap('sokolow_lyon_lvh', group='ukb_ecg_rest', tensor_from_file=_make_ukb_ecg_rest_lvh(),
                            channel_map={'no_sokolow_lyon_lvh': 0, 'Sokolow Lyon LVH': 1},
                            loss=weighted_crossentropy([0.005, 1.0], 'sokolov_lyon_lvh'))

TMAPS['ecg_rest_lvh_cornell'] = TensorMap('cornell_lvh', group='ukb_ecg_rest', tensor_from_file=_make_ukb_ecg_rest_lvh(),
                            channel_map={'no_cornell_lvh': 0, 'Cornell LVH': 1},
                            loss=weighted_crossentropy([0.003, 1.0], 'cornell_lvh'))
    

TMAPS['t2_flair_sag_p2_1mm_fs_ellip_pf78_1'] = TensorMap('t2_flair_sag_p2_1mm_fs_ellip_pf78_1', shape=(256, 256, 192), group='ukb_brain_mri',
                                                         tensor_from_file=normalized_first_date, dtype=DataSetType.FLOAT_ARRAY,
                                                         normalization={'zero_mean_std1': True})
TMAPS['t2_flair_sag_p2_1mm_fs_ellip_pf78_2'] = TensorMap('t2_flair_sag_p2_1mm_fs_ellip_pf78_2', shape=(256, 256, 192), group='ukb_brain_mri',
                                                         tensor_from_file=normalized_first_date, dtype=DataSetType.FLOAT_ARRAY,
                                                         normalization={'zero_mean_std1': True})
TMAPS['t2_flair_slice_1'] = TensorMap('t2_flair_slice_1', shape=(256, 256, 1), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                      tensor_from_file=_random_slice_tensor('t2_flair_sag_p2_1mm_fs_ellip_pf78_1'), normalization={'zero_mean_std1': True})
TMAPS['t2_flair_slice_2'] = TensorMap('t2_flair_slice_2', shape=(256, 256, 1), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                      tensor_from_file=_random_slice_tensor('t2_flair_sag_p2_1mm_fs_ellip_pf78_2'), normalization={'zero_mean_std1': True})
TMAPS['t1_p2_1mm_fov256_sag_ti_880_1'] = TensorMap('t1_p2_1mm_fov256_sag_ti_880_1', shape=(256, 256, 208), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                                   normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)
TMAPS['t1_p2_1mm_fov256_sag_ti_880_2'] = TensorMap('t1_p2_1mm_fov256_sag_ti_880_2', shape=(256, 256, 208), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                                   normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)
TMAPS['t1_slice_1'] = TensorMap('t1_slice_1', shape=(256, 256, 1), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY, normalization={'zero_mean_std1': True},
                                tensor_from_file=_random_slice_tensor('t1_p2_1mm_fov256_sag_ti_880_1'))
TMAPS['t1_slice_2'] = TensorMap('t1_slice_2', shape=(256, 256, 1), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY, normalization={'zero_mean_std1': True},
                                tensor_from_file=_random_slice_tensor('t1_p2_1mm_fov256_sag_ti_880_2'))
TMAPS['t1_20_slices_1'] = TensorMap('t1_20_slices_1', shape=(256, 256, 20), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                    normalization={'zero_mean_std1': True},
                                    tensor_from_file=_slice_subset_tensor('t1_p2_1mm_fov256_sag_ti_880_1', 94, 114))
TMAPS['t1_20_slices_2'] = TensorMap('t1_20_slices_2', shape=(256, 256, 20), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                    normalization={'zero_mean_std1': True},
                                    tensor_from_file=_slice_subset_tensor('t1_p2_1mm_fov256_sag_ti_880_2', 94, 114))
TMAPS['t2_20_slices_1'] = TensorMap('t2_20_slices_1', shape=(256, 256, 20), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                    normalization={'zero_mean_std1': True},
                                    tensor_from_file=_slice_subset_tensor('t2_flair_sag_p2_1mm_fs_ellip_pf78_1', 86, 106))
TMAPS['t2_20_slices_2'] = TensorMap('t2_20_slices_2', shape=(256, 256, 20), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                    normalization={'zero_mean_std1': True},
                                    tensor_from_file=_slice_subset_tensor('t2_flair_sag_p2_1mm_fs_ellip_pf78_2', 86, 106))
TMAPS['t1_40_slices_1'] = TensorMap('t1_40_slices_1', shape=(256, 256, 40), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                    normalization={'zero_mean_std1': True},
                                    tensor_from_file=_slice_subset_tensor('t1_p2_1mm_fov256_sag_ti_880_1', 64, 144, 2))
TMAPS['t2_40_slices_1'] = TensorMap('t2_40_slices_1', shape=(256, 256, 40), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                    normalization={'zero_mean_std1': True},
                                    tensor_from_file=_slice_subset_tensor('t2_flair_sag_p2_1mm_fs_ellip_pf78_1', 56, 136, 2))
TMAPS['sos_te1'] = TensorMap('SOS_TE1', shape=(256, 288, 48), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                             normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)
TMAPS['sos_te2'] = TensorMap('SOS_TE2', shape=(256, 288, 48), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                             normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)
TMAPS['swi'] = TensorMap('SWI', shape=(256, 288, 48), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                             normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)
TMAPS['swi_total_mag'] = TensorMap('SWI_TOTAL_MAG', shape=(256, 288, 48), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                             normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)
TMAPS['swi_total_mag_te2_orig'] = TensorMap('SWI_TOTAL_MAG_TE2_orig', shape=(256, 288, 48), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                             normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)
TMAPS['swi_total_mag_orig'] = TensorMap('SWI_TOTAL_MAG_orig', shape=(256, 288, 48), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                             normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)
TMAPS['t2star'] = TensorMap('T2star', shape=(256, 288, 48), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                             normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)
TMAPS['brain_mask_normed'] = TensorMap('brain_mask_normed', shape=(256, 288, 48), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)

TMAPS['filtered_phase'] = TensorMap('filtered_phase', shape=(256, 288, 48), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                    normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)
TMAPS['swi_to_t1_40_slices'] = TensorMap('swi_to_t1_40_slices', shape=(173, 231, 40), group='ukb_brain_mri',
                                         dtype=DataSetType.FLOAT_ARRAY, normalization={'zero_mean_std1': True},
                                         tensor_from_file=_slice_subset_tensor('SWI_TOTAL_MAG_to_T1', 60, 140, 2))
TMAPS['t2star_to_t1_40_slices'] = TensorMap('t2star_to_t1_40_slices', shape=(173, 231, 40), group='ukb_brain_mri',
                                            dtype=DataSetType.FLOAT_ARRAY, normalization={'zero_mean_std1': True},
                                            tensor_from_file=_slice_subset_tensor('T2star_to_T1', 60, 140, 2))

TMAPS['t1'] = TensorMap('T1', shape=(192, 256, 256, 1), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY, normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)
TMAPS['t1_brain'] = TensorMap('T1_brain', shape=(192, 256, 256, 1), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY, normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)
TMAPS['t1_brain_30_slices'] = TensorMap('t1_brain_30_slices', shape=(192, 256, 30), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                        normalization={'zero_mean_std1': True}, tensor_from_file=_slice_subset_tensor('T1_brain', 66, 126, 2, pad_shape=(192, 256, 256)))
TMAPS['t1_30_slices'] = TensorMap('t1_30_slices', shape=(192, 256, 30), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                  normalization={'zero_mean_std1': True}, tensor_from_file=_slice_subset_tensor('T1', 90, 150, 2, pad_shape=(192, 256, 256)))

TMAPS['t1_brain_to_mni'] = TensorMap('T1_brain_to_MNI', shape=(192, 256, 256, 1), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY, normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)
TMAPS['t1_fast_t1_brain_bias'] = TensorMap('T1_fast_T1_brain_bias', shape=(192, 256, 256, 1), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY, normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)

TMAPS['t2_flair'] = TensorMap('T2_FLAIR', shape=(192, 256, 256, 1), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY, normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)
TMAPS['t2_flair_brain'] = TensorMap('T2_FLAIR_brain', shape=(192, 256, 256, 1), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY, normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)
TMAPS['t2_flair_brain_30_slices'] = TensorMap('t2_flair_brain_30_slices', shape=(192, 256, 30), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                              normalization={'zero_mean_std1': True},
                                              tensor_from_file=_slice_subset_tensor('T2_FLAIR_brain', 66, 126, 2, pad_shape=(192, 256, 256)))
TMAPS['t2_flair_30_slices'] = TensorMap('t2_flair_30_slices', shape=(192, 256, 30), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                        normalization={'zero_mean_std1': True},
                                        tensor_from_file=_slice_subset_tensor('T2_FLAIR', 90, 150, 2, pad_shape=(192, 256, 256)))
TMAPS['t2_flair_30_slices_4d'] = TensorMap('t2_flair_30_slices_4d', shape=(192, 256, 30, 1), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                           tensor_from_file=_slice_subset_tensor('T2_FLAIR', 90, 150, 2, pad_shape=(192, 256, 256, 1)),
                                           normalization={'zero_mean_std1': True})
TMAPS['t2_flair_unbiased_brain'] = TensorMap('T2_FLAIR_unbiased_brain', shape=(192, 256, 256, 1), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY, normalization={'zero_mean_std1': True}, tensor_from_file=normalized_first_date)


def _mask_from_file(tm: TensorMap, hd5: h5py.File, dependents=None):
    original = _get_tensor_at_first_date(hd5, tm.group, DataSetType.FLOAT_ARRAY, tm.name)
    reshaped = _pad_or_crop_array_to_shape(tm.shape, original)
    tensor = to_categorical(reshaped[..., 0], tm.shape[-1])
    return tm.normalize_and_validate(tensor)


def _mask_subset_tensor(tensor_key, start, stop, step=1, pad_shape=None):
    slice_subset_tensor_from_file = _slice_subset_tensor(tensor_key, start, stop, step=step, pad_shape=pad_shape, dtype_override=DataSetType.FLOAT_ARRAY)

    def mask_subset_from_file(tm: TensorMap, hd5: h5py.File, dependents=None):
        original = slice_subset_tensor_from_file(tm, hd5, dependents)
        tensor = to_categorical(original[..., 0], tm.shape[-1])
        return tm.normalize_and_validate(tensor)
    return mask_subset_from_file


TMAPS['swi_brain_mask'] = TensorMap('SWI_brain_mask', shape=(256, 288, 48, 2), group='ukb_brain_mri', dtype=DataSetType.CATEGORICAL,
                                    tensor_from_file=_mask_from_file, channel_map={'not_brain': 0, 'brain': 1})
TMAPS['t1_brain_mask'] = TensorMap('T1_brain_mask', shape=(192, 256, 256, 2), group='ukb_brain_mri', dtype=DataSetType.CATEGORICAL,
                                   tensor_from_file=_mask_from_file, channel_map={'not_brain': 0, 'brain': 1})
TMAPS['t1_seg'] = TensorMap('T1_fast_T1_brain_seg', shape=(192, 256, 256, 4), group='ukb_brain_mri', dtype=DataSetType.CATEGORICAL,
                            tensor_from_file=_mask_from_file, channel_map={'not_brain_tissue': 0, 'csf': 1, 'grey': 2, 'white': 3})
TMAPS['t1_seg_30_slices'] = TensorMap('T1_fast_T1_brain_seg_30_slices', shape=(192, 256, 30, 4), group='ukb_brain_mri', dtype=DataSetType.CATEGORICAL,
                                      tensor_from_file=_mask_subset_tensor('T1_fast_T1_brain_seg', 90, 150, 2, pad_shape=(192, 256, 256, 1)),
                                      channel_map={'not_brain_tissue': 0, 'csf': 1, 'grey': 2, 'white': 3})
TMAPS['t1_brain_mask_30_slices'] = TensorMap('T1_brain_mask_30_slices', shape=(192, 256, 30, 2), group='ukb_brain_mri', dtype=DataSetType.CATEGORICAL,
                                             tensor_from_file=_mask_subset_tensor('T1_brain_mask', 90, 150, 2, pad_shape=(192, 256, 256, 1)),
                                             channel_map={'not_brain': 0, 'brain': 1})
TMAPS['lesions'] = TensorMap('lesions_final_mask', shape=(192, 256, 256, 2), group='ukb_brain_mri', dtype=DataSetType.CATEGORICAL,
                             tensor_from_file=_mask_from_file, channel_map={'not_lesion': 0, 'lesion': 1}, loss=weighted_crossentropy([0.01, 10.0], 'lesion'))


def _combined_subset_tensor(tensor_keys, start, stop, step=1, pad_shape=None):
    slice_subsets = [_slice_subset_tensor(k, start, stop, step=step, pad_shape=pad_shape, allow_channels=False) for k in tensor_keys]

    def mask_subset_from_file(tm: TensorMap, hd5: h5py.File, dependents=None):
        tensor = np.zeros(tm.shape, dtype=np.float32)
        for i, slice_subset_tensor_from_file in enumerate(slice_subsets):
            tensor[..., i] = slice_subset_tensor_from_file(tm, hd5, dependents)
        return tm.normalize_and_validate(tensor)
    return mask_subset_from_file


TMAPS['t1_and_t2_flair_30_slices'] = TensorMap('t1_and_t2_flair_30_slices', shape=(192, 256, 30, 2), group='ukb_brain_mri', dtype=DataSetType.FLOAT_ARRAY,
                                               tensor_from_file=_combined_subset_tensor(['T1', 'T2_FLAIR'], 90, 150, 2, pad_shape=(192, 256, 256)),
                                               normalization={'zero_mean_std1': True})


def _ttn_tensor_from_file(tm, hd5, dependents={}):
    index = 0
    categorical_data = np.zeros(tm.shape, dtype=np.float32)
    if 'has_exome' not in hd5['categorical']:
        raise ValueError('Skipping people without exome sequencing.')
    if tm.name in hd5['categorical'] and int(hd5['categorical'][tm.name][0]) != 0:
        index = 1
    categorical_data[index] = 1.0
    return categorical_data


TMAPS['ttntv'] = TensorMap('has_ttntv', group='categorical_flag', channel_map={'no_TTN_tv': 0, 'TTN_tv': 1}, tensor_from_file=_ttn_tensor_from_file)
TMAPS['ttntv_10x'] = TensorMap('has_ttntv', group='categorical_flag', channel_map={'no_TTN_tv': 0, 'TTN_tv': 1}, loss_weight=10.0, tensor_from_file=_ttn_tensor_from_file)


def _make_index_tensor_from_file(index_map_name):
    def indexed_lvmass_tensor_from_file(tm, hd5, dependents={}):
        tensor = np.zeros(tm.shape, dtype=np.float32)
        for k in tm.channel_map:
            tensor = np.array(hd5[tm.group][k], dtype=np.float32)
        index = np.array(hd5[tm.group][index_map_name], dtype=np.float32)
        return tm.normalize_and_validate(tensor / index)
    return indexed_lvmass_tensor_from_file


TMAPS['lv_mass_dubois_index'] = TensorMap('lv_mass_dubois_index', group='continuous', activation='linear', loss='logcosh', loss_weight=1.0,
                                          tensor_from_file=_make_index_tensor_from_file('bsa_dubois'),
                                          channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8})
TMAPS['lv_mass_mosteller_index'] = TensorMap('lv_mass_mosteller_index', group='continuous', activation='linear', loss='logcosh', loss_weight=1.0,
                                             tensor_from_file=_make_index_tensor_from_file('bsa_mosteller'),
                                             channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8})
TMAPS['lv_mass_dubois_index_sentinel'] = TensorMap('lv_mass_dubois_index', group='continuous', activation='linear', sentinel=0, loss_weight=1.0,
                                                   tensor_from_file=_make_index_tensor_from_file('bsa_dubois'),
                                                   channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8})
TMAPS['lv_mass_mosteller_index_sentinel'] = TensorMap('lv_mass_mosteller_index', group='continuous', activation='linear', sentinel=0, loss_weight=1.0,
                                                      tensor_from_file=_make_index_tensor_from_file('bsa_mosteller'),
                                                      channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8})
TMAPS['lv_mass_dubois_indexp'] = TensorMap('lv_mass_dubois_index', group='continuous', activation='linear', loss='logcosh', loss_weight=1.0,
                                           parents=['output_mri_systole_diastole_8_segmented_categorical'],
                                           tensor_from_file=_make_index_tensor_from_file('bsa_dubois'),
                                           channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8})
TMAPS['lv_mass_mosteller_indexp'] = TensorMap('lv_mass_mosteller_index', group='continuous', activation='linear', loss='logcosh', loss_weight=1.0,
                                              parents=['output_mri_systole_diastole_8_segmented_categorical'],
                                              tensor_from_file=_make_index_tensor_from_file('bsa_mosteller'),
                                              channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8})
TMAPS['lvm_dubois_index'] = TensorMap('lvm_dubois_index', group='continuous', activation='linear', loss='logcosh', loss_weight=1.0,
                                      tensor_from_file=_make_index_tensor_from_file('bsa_dubois'),
                                      channel_map={'LVM': 0}, normalization={'mean': 89.7, 'std': 24.8})
TMAPS['lvm_mosteller_index'] = TensorMap('lvm_mosteller_index', group='continuous', activation='linear', loss='logcosh', loss_weight=1.0,
                                         tensor_from_file=_make_index_tensor_from_file('bsa_mosteller'),
                                         channel_map={'LVM': 0}, normalization={'mean': 89.7, 'std': 24.8})
TMAPS['lvm_dubois_index_w4'] = TensorMap('lvm_dubois_index', group='continuous', activation='linear', loss='logcosh', loss_weight=4.0,
                                          tensor_from_file=_make_index_tensor_from_file('bsa_dubois'),
                                          channel_map={'LVM': 0}, normalization={'mean': 89.7, 'std': 24.8})
TMAPS['lvm_mosteller_index_w4'] = TensorMap('lvm_mosteller_index', group='continuous', activation='linear', loss='logcosh', loss_weight=4.0,
                                             tensor_from_file=_make_index_tensor_from_file('bsa_mosteller'),
                                             channel_map={'LVM': 0}, normalization={'mean': 89.7, 'std': 24.8})
TMAPS['lvm_dubois_index_sentinel'] = TensorMap('lvm_dubois_index', group='continuous', activation='linear', sentinel=0, loss_weight=1.0,
                                               tensor_from_file=_make_index_tensor_from_file('bsa_dubois'),
                                               channel_map={'LVM': 0}, normalization={'mean': 89.7, 'std': 24.8})
TMAPS['lvm_mosteller_index_sentinel'] = TensorMap('lvm_mosteller_index', group='continuous', activation='linear', sentinel=0, loss_weight=1.0,
                                                  tensor_from_file=_make_index_tensor_from_file('bsa_mosteller'),
                                                  channel_map={'LVM': 0}, normalization={'mean': 89.7, 'std': 24.8})


def _select_tensor_from_file(selection_predicate: Callable):
    def selected_tensor_from_file(tm, hd5, dependents={}):
        if not selection_predicate(hd5):
            raise ValueError(f'Tensor did not meet selection criteria:{selection_predicate.__name__} with Tensor Map:{tm.name}')
        tensor = np.zeros(tm.shape, dtype=np.float32)
        for k in tm.channel_map:
            tensor = np.array(hd5[tm.group][k], dtype=np.float32)
        return tm.normalize_and_validate(tensor)
    return selected_tensor_from_file


def _is_genetic_man(hd5):
    return 'Genetic-sex_Male_0_0' in hd5['categorical']


def _is_genetic_woman(hd5):
    return 'Genetic-sex_Female_0_0' in hd5['categorical']


TMAPS['myocardial_mass_noheritable_men_only'] = TensorMap('inferred_myocardial_mass_noheritable', group='continuous', activation='linear', loss='logcosh',
                                                          tensor_from_file=_select_tensor_from_file(_is_genetic_man),
                                                          channel_map={'inferred_myocardial_mass_noheritable': 0}, normalization={'mean': 100.0, 'std': 18.0})
TMAPS['myocardial_mass_noheritable_women_only'] = TensorMap('inferred_myocardial_mass_noheritable', group='continuous', activation='linear', loss='logcosh',
                                                            tensor_from_file=_select_tensor_from_file(_is_genetic_woman),
                                                            channel_map={'inferred_myocardial_mass_noheritable': 0}, normalization={'mean': 78.0, 'std': 16.0})


def _make_lvh_from_lvm_tensor_from_file(bsa_key, lvm_key, group='continuous', male_lvh_threshold=72, female_lvh_threshold=55):
    def lvh_from_lvm_tensor_from_file(tm, hd5, dependents={}):
        tensor = np.zeros(tm.shape, dtype=np.float32)
        lvm = float(hd5[group][lvm_key])
        bsa = float(hd5[group][bsa_key])
        lvm_indexed = lvm / bsa
        index = 0
        if _is_genetic_man(hd5) and lvm_indexed > male_lvh_threshold:
            index = 1
        elif _is_genetic_woman(hd5) and lvm_indexed > female_lvh_threshold:
            index = 1
        tensor[index] = 1
        return tensor
    return lvh_from_lvm_tensor_from_file


TMAPS['lvh_from_lvm_actual'] = TensorMap('lvh_from_lvm_actual', group='categorical', channel_map={'no_lvh': 0, 'lvh': 1},
                                         tensor_from_file=_make_lvh_from_lvm_tensor_from_file('bsa_dubois', 'myocardial_mass_noheritable_sentinel_actual'))
TMAPS['lvh_from_lvm_predict'] = TensorMap('lvh_from_lvm_predict', group='categorical', channel_map={'no_lvh': 0, 'lvh': 1},
                                          tensor_from_file=_make_lvh_from_lvm_tensor_from_file('bsa_dubois', 'myocardial_mass_noheritable_sentinel_prediction'))


def _mri_slice_blackout_tensor_from_file(tm, hd5, dependents={}):
    cur_slice = np.random.choice(list(hd5[MRI_TO_SEGMENT].keys()))
    tensor = np.zeros(tm.shape, dtype=np.float32)
    dependents[tm.dependent_map] = np.zeros(tm.dependent_map.shape, dtype=np.float32)
    tensor[:, :, 0] = np.array(hd5[MRI_TO_SEGMENT][cur_slice], dtype=np.float32)
    label_tensor = np.array(hd5[MRI_SEGMENTED][cur_slice], dtype=np.float32)
    dependents[tm.dependent_map][:, :, :] = to_categorical(label_tensor, tm.dependent_map.shape[-1])
    tensor[:, :, 0] *= np.not_equal(label_tensor, 0, dtype=np.float32)
    return tm.zero_mean_std1(tensor)


TMAPS['mri_slice_blackout_segmented_weighted'] = TensorMap('mri_slice_segmented', (256, 256, 3), group='categorical', channel_map=MRI_SEGMENTED_CHANNEL_MAP,
                                                           loss=weighted_crossentropy([0.1, 25.0, 25.0], 'mri_slice_blackout_segmented'))
TMAPS['mri_slice_blackout'] = TensorMap('mri_slice_blackout', (256, 256, 1), tensor_from_file=_mri_slice_blackout_tensor_from_file,
                                        dependent_map=TMAPS['mri_slice_blackout_segmented_weighted'])


def _slice_tensor(tensor_key, slice_index):
    def _slice_tensor_from_file(tm, hd5, dependents={}):
        tensor = np.zeros(tm.shape, dtype=np.float32)
        tensor[..., 0] = np.array(hd5[tensor_key][slice_index], dtype=np.float32)
        return tm.normalize_and_validate(tensor)
    return _slice_tensor_from_file


TMAPS['lax_4ch_diastole_slice'] = TensorMap('lax_4ch_diastole_slice', (256, 256, 1), group='root_array', loss='logcosh',
                                            tensor_from_file=_slice_tensor('cine_segmented_lax_4ch', 0),
                                            normalization={'zero_mean_std1': True})


def _make_fallback_tensor_from_file(tensor_keys):
    def fallback_tensor_from_file(tm, hd5, dependents={}):
        for k in tensor_keys:
            if k in hd5:
                tensor = np.array(hd5[k], dtype=np.float32)
                return tm.normalize_and_validate(tensor)
        raise ValueError(f'No fallback tensor found from keys: {tensor_keys}')
    return fallback_tensor_from_file


TMAPS['shmolli_192i_both'] = TensorMap('shmolli_192i', (288, 384, 7), group='root_array',
                                       tensor_from_file=_make_fallback_tensor_from_file(['shmolli_192i', 'shmolli_192i_liver']))
TMAPS['shmolli_192i_liver_only'] = TensorMap('shmolli_192i', (288, 384, 7), group='root_array',
                                             tensor_from_file=_make_fallback_tensor_from_file(['shmolli_192i_liver']))
