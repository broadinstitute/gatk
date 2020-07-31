import os
import csv
import copy
import h5py
import logging
import datetime
import numpy as np
import pandas as pd
from itertools import product
from collections import defaultdict
from typing import Callable, Dict, List, Tuple, Union

from ml4cvd.tensor_maps_by_hand import TMAPS
from ml4cvd.defines import ECG_REST_AMP_LEADS, PARTNERS_DATE_FORMAT, STOP_CHAR, PARTNERS_DATETIME_FORMAT, CARDIAC_SURGERY_DATE_FORMAT
from ml4cvd.TensorMap import TensorMap, str2date, Interpretation, make_range_validator, decompress_data, TimeSeriesOrder
from ml4cvd.normalizer import Standardize, ZeroMeanStd1


YEAR_DAYS = 365.26
INCIDENCE_CSV = '/media/erisone_snf13/lc_outcomes.csv'
CARDIAC_SURGERY_OUTCOMES_CSV = '/data/sts-data/mgh-preop-ecg-outcome-labels.csv'
PARTNERS_PREFIX = 'partners_ecg_rest'


def _hd5_filename_to_mrn_int(filename: str) -> int:
    return int(os.path.basename(filename).split('.')[0])


def _get_ecg_dates(tm, hd5):
    if not hasattr(_get_ecg_dates, 'mrn_lookup'):
        _get_ecg_dates.mrn_lookup = dict()
    mrn = _hd5_filename_to_mrn_int(hd5.filename)
    if mrn in _get_ecg_dates.mrn_lookup:
        return _get_ecg_dates.mrn_lookup[mrn]

    dates = list(hd5[tm.path_prefix])
    if tm.time_series_lookup is not None:
        start, end = tm.time_series_lookup[mrn]
        dates = [date for date in dates if start < date < end]
    if tm.time_series_order == TimeSeriesOrder.NEWEST:
        dates.sort()
    elif tm.time_series_order == TimeSeriesOrder.OLDEST:
        dates.sort(reverse=True)
    elif tm.time_series_order == TimeSeriesOrder.RANDOM:
        np.random.shuffle(dates)
    else:
        raise NotImplementedError(f'Unknown option "{tm.time_series_order}" passed for which tensors to use in multi tensor HD5')
    start_idx = tm.time_series_limit if tm.time_series_limit is not None else 1
    dates = dates[-start_idx:]  # If num_tensors is 0, get all tensors
    dates.sort(reverse=True)
    _get_ecg_dates.mrn_lookup[mrn] = dates
    return dates


def validator_no_empty(tm: TensorMap, tensor: np.ndarray, hd5: h5py.File):
    if any(tensor == ''):
        raise ValueError(f'TensorMap {tm.name} failed empty string check.')


def validator_no_negative(tm: TensorMap, tensor: np.ndarray, hd5: h5py.File):
    if any(tensor < 0):
        raise ValueError(f'TensorMap {tm.name} failed non-negative check')


def validator_not_all_zero(tm: TensorMap, tensor: np.ndarray, hd5: h5py.File):
    if np.count_nonzero(tensor) == 0:
        raise ValueError(f'TensorMap {tm.name} failed all-zero check')


def _is_dynamic_shape(tm: TensorMap, num_ecgs: int) -> Tuple[bool, Tuple[int, ...]]:
    if tm.shape[0] is None:
        return True, (num_ecgs,) + tm.shape[1:]
    return False, tm.shape


def _make_hd5_path(tm, ecg_date, value_key):
    return f'{tm.path_prefix}/{ecg_date}/{value_key}'


def _resample_voltage(voltage, desired_samples):
    if len(voltage) == desired_samples:
        return voltage
    elif len(voltage) == 2500 and desired_samples == 5000:
        x = np.arange(2500)
        x_interp = np.linspace(0, 2500, 5000)
        return np.interp(x_interp, x, voltage)
    elif len(voltage) == 5000 and desired_samples == 2500:
        return voltage[::2]
    else:
        raise ValueError(f'Voltage length {len(voltage)} is not desired {desired_samples} and re-sampling method is unknown.')


def _resample_voltage_with_rate(voltage, desired_samples, rate, desired_rate):
    if len(voltage) == desired_samples and rate == desired_rate:
        return voltage
    elif desired_samples / len(voltage) == 2 and desired_rate / rate == 2:
        x = np.arange(len(voltage))
        x_interp = np.linspace(0, len(voltage), desired_samples)
        return np.interp(x_interp, x, voltage)
    elif desired_samples / len(voltage) == 0.5 and desired_rate / rate == 0.5:
        return voltage[::2]
    elif desired_samples / len(voltage) == 2 and desired_rate == rate:
        return np.pad(voltage, (0, len(voltage)))
    elif desired_samples / len(voltage) == 0.5 and desired_rate == rate:
        return voltage[:len(voltage)//2]
    else:
        raise ValueError(f'Voltage length {len(voltage)} is not desired {desired_samples} with desired rate {desired_rate} and rate {rate}.')


def make_voltage(exact_length = False):
    def get_voltage_from_file(tm, hd5, dependents={}):
        ecg_dates = _get_ecg_dates(tm, hd5)
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        voltage_length = shape[1] if dynamic else shape[0]
        tensor = np.zeros(shape, dtype=np.float32)
        for i, ecg_date in enumerate(ecg_dates):
            for cm in tm.channel_map:
                try:
                    path = _make_hd5_path(tm, ecg_date, cm)
                    voltage = decompress_data(data_compressed=hd5[path][()], dtype=hd5[path].attrs['dtype'])
                    if exact_length:
                        assert len(voltage) == voltage_length
                    voltage = _resample_voltage(voltage, voltage_length)
                    slices = (i, ..., tm.channel_map[cm]) if dynamic else (..., tm.channel_map[cm])
                    tensor[slices] = voltage
                except (KeyError, AssertionError, ValueError):
                    logging.debug(f'Could not get voltage for lead {cm} with {voltage_length} samples in {hd5.filename}')
        return tensor
    return get_voltage_from_file

# Creates 12 TMaps:
# partners_ecg_2500      partners_ecg_2500_exact      partners_ecg_5000      partners_ecg_5000_exact
# partners_ecg_2500_std  partners_ecg_2500_std_exact  partners_ecg_5000_std  partners_ecg_5000_std_exact
# partners_ecg_2500_raw  partners_ecg_2500_raw_exact  partners_ecg_5000_raw  partners_ecg_5000_raw_exact
#
# default normalizes with ZeroMeanStd1 and resamples
# _std normalizes with Standardize mean = 0, std = 2000
# _raw does not normalize
# _exact does not resample
length_options = [2500, 5000]
exact_options = [True, False]
normalize_options = [ZeroMeanStd1(), Standardize(mean=0, std=2000), None]
for length, exact_length, normalization in product(length_options, exact_options, normalize_options):
    norm = '' if isinstance(normalization, ZeroMeanStd1) else '_std' if isinstance(normalization, Standardize) else '_raw'
    exact = '_exact' if exact_length else ''
    name = f'partners_ecg_{length}{norm}{exact}'
    TMAPS[name] = TensorMap(
        name,
        shape=(None, length, 12),
        path_prefix=PARTNERS_PREFIX,
        tensor_from_file=make_voltage(exact_length),
        normalization=normalization,
        channel_map=ECG_REST_AMP_LEADS,
        time_series_limit=0,
        validator=validator_not_all_zero,
    )


def voltage_stat(tm, hd5, dependents={}):
    ecg_dates = _get_ecg_dates(tm, hd5)
    dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
    tensor = np.zeros(shape, dtype=np.float32)
    for i, ecg_date in enumerate(ecg_dates):
        try:
            slices = lambda stat: (i, tm.channel_map[stat]) if dynamic else (tm.channel_map[stat],)
            path = lambda lead: _make_hd5_path(tm, ecg_date, lead)
            voltages = np.array([decompress_data(data_compressed=hd5[path(lead)][()], dtype='int16') for lead in ECG_REST_AMP_LEADS])
            tensor[slices('mean')] = np.mean(voltages)
            tensor[slices('std')] = np.std(voltages)
            tensor[slices('min')] = np.min(voltages)
            tensor[slices('max')] = np.max(voltages)
            tensor[slices('median')] = np.median(voltages)
        except KeyError:
            logging.warning(f'Could not get voltage stats for ECG at {hd5.filename}')
    return tensor


TMAPS['partners_ecg_voltage_stats'] = TensorMap(
    'partners_ecg_voltage_stats',
    shape=(None, 5),
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=voltage_stat,
    channel_map={'mean': 0, 'std': 1, 'min': 2, 'max': 3, 'median': 4},
    time_series_limit=0,
)


def make_voltage_attr(volt_attr: str = ""):
    def get_voltage_attr_from_file(tm, hd5, dependents={}):
        ecg_dates = _get_ecg_dates(tm, hd5)
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        tensor = np.zeros(shape, dtype=np.float32)
        for i, ecg_date in enumerate(ecg_dates):
            for cm in tm.channel_map:
                try:
                    path = _make_hd5_path(tm, ecg_date, cm)
                    slices = (i, tm.channel_map[cm]) if dynamic else (tm.channel_map[cm],)
                    tensor[slices] = hd5[path].attrs[volt_attr]
                except KeyError:
                    pass
        return tensor
    return get_voltage_attr_from_file


TMAPS["voltage_len"] = TensorMap(
    "voltage_len",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_voltage_attr(volt_attr="len"),
    shape=(None, 12),
    channel_map=ECG_REST_AMP_LEADS,
    time_series_limit=0,
)


def make_partners_ecg_label(keys: Union[str, List[str]] = "read_md_clean", dict_of_list: Dict = dict(), not_found_key: str = "unspecified"):
    if type(keys) == str:
        keys = [keys]

    def get_partners_ecg_label(tm, hd5, dependents={}):
        ecg_dates = _get_ecg_dates(tm, hd5)
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        label_array = np.zeros(shape, dtype=np.float32)
        for i, ecg_date in enumerate(ecg_dates):
            found = False
            for channel, idx in sorted(tm.channel_map.items(), key=lambda cm: cm[1]):
                if channel not in dict_of_list:
                    continue
                for key in keys:
                    path = _make_hd5_path(tm, ecg_date, key)
                    if path not in hd5:
                        continue
                    read = decompress_data(data_compressed=hd5[path][()], dtype=hd5[path].attrs['dtype'])
                    for string in dict_of_list[channel]:
                        if string not in read:
                            continue
                        slices = (i, idx) if dynamic else (idx,)
                        label_array[slices] = 1
                        found = True
                        break
                    if found:
                        break
                if found:
                    break
            if not found:
                slices = (i, tm.channel_map[not_found_key]) if dynamic else (tm.channel_map[not_found_key],)
                label_array[slices] = 1
        return label_array
    return get_partners_ecg_label


def partners_ecg_datetime(tm, hd5, dependents={}):
    ecg_dates = _get_ecg_dates(tm, hd5)
    dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
    tensor = np.full(shape , '', dtype=f'<U19')
    for i, ecg_date in enumerate(ecg_dates):
        tensor[i] = ecg_date
    return tensor


task = "partners_ecg_datetime"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=partners_ecg_datetime,
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
)


def make_voltage_len_categorical_tmap(lead, channel_prefix = '_', channel_unknown = 'other'):
    def _tensor_from_file(tm, hd5, dependents = {}):
        ecg_dates = _get_ecg_dates(tm, hd5)
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        tensor = np.zeros(shape, dtype=float)
        for i, ecg_date in enumerate(ecg_dates):
            path = _make_hd5_path(tm, ecg_date, lead)
            try:
                lead_len = hd5[path].attrs['len']
                lead_len = f'{channel_prefix}{lead_len}'
                matched = False
                for cm in tm.channel_map:
                    if lead_len.lower() == cm.lower():
                        slices = (i, tm.channel_map[cm]) if dynamic else (tm.channel_map[cm],)
                        tensor[slices] = 1.0
                        matched = True
                        break
                if not matched:
                    slices = (i, tm.channel_map[channel_unknown]) if dynamic else (tm.channel_map[channel_unknown],)
                    tensor[slices] = 1.0
            except KeyError:
                logging.debug(f'Could not get voltage length for lead {lead} from ECG on {ecg_date} in {hd5.filename}')
        return tensor
    return _tensor_from_file


for lead in ECG_REST_AMP_LEADS:
    tmap_name = f'lead_{lead}_len'
    TMAPS[tmap_name] = TensorMap(
        tmap_name,
        interpretation=Interpretation.CATEGORICAL,
        path_prefix=PARTNERS_PREFIX,
        tensor_from_file=make_voltage_len_categorical_tmap(lead=lead),
        channel_map={'_2500': 0, '_5000': 1, 'other': 2},
        time_series_limit=0,
        validator=validator_not_all_zero,
    )


def make_partners_ecg_tensor(key: str, fill: float = 0, channel_prefix: str = '', channel_unknown: str = 'other'):
    def get_partners_ecg_tensor(tm, hd5, dependents={}):
        ecg_dates = _get_ecg_dates(tm, hd5)
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        if tm.interpretation == Interpretation.LANGUAGE:
            tensor = np.full(shape, '', dtype=object)
        elif tm.interpretation == Interpretation.CONTINUOUS:
            tensor = np.zeros(shape, dtype=np.float32) if fill == 0 else np.full(shape, fill, dtype=np.float32)
        elif tm.interpretation == Interpretation.CATEGORICAL:
            tensor = np.zeros(shape, dtype=float)
        else:
            raise NotImplementedError(f'unsupported interpretation for partners tmaps: {tm.interpretation}')

        for i, ecg_date in enumerate(ecg_dates):
            path = _make_hd5_path(tm, ecg_date, key)
            try:
                data = decompress_data(data_compressed=hd5[path][()], dtype=hd5[path].attrs['dtype'])
                if tm.interpretation == Interpretation.CATEGORICAL:
                    matched = False
                    data = f'{channel_prefix}{data}'
                    for cm in tm.channel_map:
                        if data.lower() == cm.lower():
                            slices = (i, tm.channel_map[cm]) if dynamic else (tm.channel_map[cm],)
                            tensor[slices] = 1.0
                            matched = True
                            break
                    if not matched:
                        slices = (i, tm.channel_map[channel_unknown]) if dynamic else (tm.channel_map[channel_unknown],)
                        tensor[slices] = 1.0
                else:
                    tensor[i] = data
            except (KeyError, ValueError):
                logging.debug(f'Could not obtain tensor {tm.name} from ECG on {ecg_date} in {hd5.filename}')

        if tm.interpretation == Interpretation.LANGUAGE:
            tensor = tensor.astype(str)
        return tensor
    return get_partners_ecg_tensor


def make_partners_language_tensor(key: str):
    def language_tensor(tm, hd5, dependents={}):
        words = str(decompress_data(data_compressed=hd5[key][()], dtype=hd5[key].attrs['dtype']))
        tensor = np.zeros(tm.shape, dtype=np.float32)
        for i, c in enumerate(words):
            if i >= tm.shape[0]:
                logging.debug(f'Text {words} is longer than {tm.name} can store in shape:{tm.shape}, truncating...')
                break
            tensor[i, tm.channel_map[c]] = 1.0
        tensor[min(tm.shape[0]-1, i+1), tm.channel_map[STOP_CHAR]] = 1.0
        return tensor
    return language_tensor


task = "partners_ecg_read_md"
TMAPS[task] = TensorMap(
    task,
    #annotation_units=128,
    #channel_map=PARTNERS_CHAR_2_IDX,
    interpretation=Interpretation.LANGUAGE,
    #shape=(512, len(PARTNERS_CHAR_2_IDX)),
    path_prefix=PARTNERS_PREFIX,
    #tensor_from_file=make_partners_language_tensor(key="read_md_clean"),
    tensor_from_file=make_partners_ecg_tensor(key="read_md_clean"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
)


task = "partners_ecg_read_pc"
TMAPS[task] = TensorMap(
    task,
    #annotation_units=128,
    #channel_map=PARTNERS_CHAR_2_IDX,
    interpretation=Interpretation.LANGUAGE,
    #tensor_from_file=make_partners_language_tensor(key="read_pc_clean"),
    #shape=(512, len(PARTNERS_CHAR_2_IDX)),
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="read_pc_clean"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
)


task = "partners_ecg_patientid"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="patientid"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
)


def validator_clean_mrn(tm: TensorMap, tensor: np.ndarray, hd5: h5py.File):
    int(tensor)


task = "partners_ecg_patientid_clean"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="patientid_clean"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_clean_mrn,
)


task = "partners_ecg_firstname"
TMAPS[task] = TensorMap(
    task,
    #channel_map=PARTNERS_CHAR_2_IDX,
    interpretation=Interpretation.LANGUAGE,
    #tensor_from_file=make_partners_language_tensor(key="patientfirstname"),
    #shape=(512, len(PARTNERS_CHAR_2_IDX)),
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="patientfirstname"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
)


task = "partners_ecg_lastname"
TMAPS[task] = TensorMap(
    task,
    #channel_map=PARTNERS_CHAR_2_IDX,
    interpretation=Interpretation.LANGUAGE,
    #tensor_from_file=make_partners_language_tensor(key="patientlastname"),
    #shape=(512, len(PARTNERS_CHAR_2_IDX)),
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="patientlastname"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
)


task = "partners_ecg_sex"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="gender"),
    channel_map={'female': 0, 'male': 1},
    time_series_limit=0,
    validator=validator_not_all_zero,
)

task = "partners_ecg_date"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="acquisitiondate"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
)


task = "partners_ecg_time"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="acquisitiontime"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
)


task = "partners_ecg_sitename"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="sitename"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
)


task = "partners_ecg_location"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="location"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
)


task = "partners_ecg_dob"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="dateofbirth"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
)


def make_sampling_frequency_from_file(lead: str = "I", duration: int = 10, channel_prefix: str = "_", channel_unknown: str = "other", fill: int = -1):
    def sampling_frequency_from_file(tm: TensorMap, hd5: h5py.File, dependents: Dict = {}):
        ecg_dates = _get_ecg_dates(tm, hd5)
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        if tm.interpretation == Interpretation.CATEGORICAL:
            tensor = np.zeros(shape, dtype=np.float32)
        else:
            tensor = np.full(shape, fill, dtype=np.float32)
        for i, ecg_date in enumerate(ecg_dates):
            path = _make_hd5_path(tm, ecg_date, lead)
            lead_length = hd5[path].attrs["len"]
            sampling_frequency = lead_length / duration
            try:
                if tm.interpretation == Interpretation.CATEGORICAL:
                    matched = False
                    sampling_frequency = f'{channel_prefix}{sampling_frequency}'
                    for cm in tm.channel_map:
                        if sampling_frequency.lower() == cm.lower():
                            slices = (i, tm.channel_map[cm]) if dynamic else (tm.channel_map[cm],)
                            tensor[slices] = 1.0
                            matched = True
                            break
                    if not matched:
                        slices = (i, tm.channel_map[channel_unknown]) if dynamic else (tm.channel_map[channel_unknown],)
                        tensor[slices] = 1.0
                else:
                    tensor[i] = sampling_frequency
            except (KeyError, ValueError):
                logging.debug(f'Could not calculate sampling frequency from ECG on {ecg_date} in {hd5.filename}')
        return tensor
    return sampling_frequency_from_file


task = "partners_ecg_sampling_frequency"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_sampling_frequency_from_file(),
    channel_map={'_250': 0, '_500': 1, 'other': 2},
    time_series_limit=0,
    validator=validator_not_all_zero,
)


task = "partners_ecg_sampling_frequency_pc"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="ecgsamplebase_pc", channel_prefix='_'),
    channel_map={'_0': 0, '_250': 1, '_500': 2, 'other': 3},
    time_series_limit=0,
    validator=validator_not_all_zero,
)


task = "partners_ecg_sampling_frequency_md"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="ecgsamplebase_md", channel_prefix='_'),
    channel_map={'_0': 0, '_250': 1, '_500': 2, 'other': 3},
    time_series_limit=0,
    validator=validator_not_all_zero,
)


task = "partners_ecg_sampling_frequency_wv"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="waveform_samplebase", channel_prefix='_'),
    channel_map={'_0': 0, '_240': 1, '_250': 2, '_500': 3, 'other': 4},
    time_series_limit=0,
    validator=validator_not_all_zero,
)


task = "partners_ecg_sampling_frequency_continuous"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_sampling_frequency_from_file(),
    time_series_limit=0,
    shape=(None, 1),
    validator=validator_no_negative,
)


task = "partners_ecg_sampling_frequency_pc_continuous"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="ecgsamplebase_pc", fill=-1),
    time_series_limit=0,
    shape=(None, 1),
    validator=validator_no_negative,
)


task = "partners_ecg_sampling_frequency_md_continuous"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="ecgsamplebase_md", fill=-1),
    time_series_limit=0,
    shape=(None, 1),
    validator=validator_no_negative,
)


task = "partners_ecg_sampling_frequency_wv_continuous"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="waveform_samplebase", fill=-1),
    time_series_limit=0,
    shape=(None, 1),
    validator=validator_no_negative,
)


task = "partners_ecg_time_resolution"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="intervalmeasurementtimeresolution", channel_prefix='_'),
    channel_map={'_25': 0, '_50': 1, '_100': 2, 'other': 3},
    time_series_limit=0,
    validator=validator_not_all_zero,
)


task = "partners_ecg_amplitude_resolution"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="intervalmeasurementamplituderesolution", channel_prefix='_'),
    channel_map={'_10': 0, '_20': 1, '_40': 2, 'other': 3},
    time_series_limit=0,
    validator=validator_not_all_zero,
)


task = "partners_ecg_measurement_filter"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="intervalmeasurementfilter", channel_prefix='_'),
    time_series_limit=0,
    channel_map={'_None': 0, '_40': 1, '_80': 2, 'other': 3},
    validator=validator_not_all_zero,
)


task = "partners_ecg_high_pass_filter"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="waveform_highpassfilter", fill=-1),
    time_series_limit=0,
    shape=(None, 1),
    validator=validator_no_negative,
)


task = "partners_ecg_low_pass_filter"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="waveform_lowpassfilter", fill=-1),
    time_series_limit=0,
    shape=(None, 1),
    validator=validator_no_negative,
)


task = "partners_ecg_ac_filter"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="waveform_acfilter", channel_prefix='_'),
    time_series_limit=0,
    channel_map={'_None': 0, '_50': 1, '_60': 2, 'other': 3},
    validator=validator_not_all_zero,
)


task = "partners_ecg_rate_pc"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="ventricularrate_pc"),
    shape=(None, 1),
    normalization=Standardize(mean=59.3, std=10.6),
    time_series_limit=0,
    validator=make_range_validator(10, 200),
)


task = "partners_ecg_rate_md"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="ventricularrate_md"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(10, 200),
)


task = "partners_ecg_qrs_pc"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="qrsduration_pc"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(20, 400),
)


task = "partners_ecg_qrs_md"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="qrsduration_md"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(20, 400),
)


task = "partners_ecg_pr_pc"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="printerval_pc"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(50, 500),
)


task = "partners_ecg_pr_md"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="printerval_md"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(50, 500),
)


task = "partners_ecg_qt_pc"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="qtinterval_pc"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(100, 800),
)


task = "partners_ecg_qt_md"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="qtinterval_md"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(100, 800),
)


task = "partners_ecg_qtc_pc"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="qtcorrected_pc"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(100, 800),
)


task = "partners_ecg_qtc_md"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="qtcorrected_md"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(100, 800),
)


task = "partners_ecg_paxis_pc"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="paxis_pc", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_paxis_md"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="paxis_md", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_raxis_pc"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="raxis_pc", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_raxis_md"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="raxis_md", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_taxis_pc"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="taxis_pc", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_taxis_md"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="taxis_md", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_measuredamplitudepeak_r"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="measuredamplitudepeak_IE_R", fill=np.nan),
    shape=(None, 12),
    time_series_limit=0,
)


task = "partners_ecg_acquisitiondevice"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="acquisitiondevice", fill=999),
    shape=(None, 1),
    time_series_limit=0,
)


task = "partners_ecg_weight_lbs"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="weightlbs"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(100, 800),
)


def partners_ecg_age(tm, hd5, dependents={}):
    ecg_dates = _get_ecg_dates(tm, hd5)
    dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
    tensor = np.zeros(shape, dtype=float)
    for i, ecg_date in enumerate(ecg_dates):
        if i >= shape[0]:
            break
        path = lambda key: _make_hd5_path(tm, ecg_date, key)
        try:
            birthday = decompress_data(data_compressed=hd5[path('dateofbirth')][()], dtype='str')
            acquisition = decompress_data(data_compressed=hd5[path('acquisitiondate')][()], dtype='str')
            delta = _partners_str2date(acquisition) - _partners_str2date(birthday)
            years = delta.days / YEAR_DAYS
            tensor[i] = years
        except KeyError:
            try:
                tensor[i] = decompress_data(data_compressed=hd5[path('patientage')][()], dtype='str')
            except KeyError:
                logging.debug(f'Could not get patient date of birth or age from ECG on {ecg_date} in {hd5.filename}')
    return tensor


TMAPS['partners_ecg_age'] = TensorMap('partners_ecg_age', path_prefix=PARTNERS_PREFIX, loss='logcosh', tensor_from_file=partners_ecg_age, shape=(None, 1), time_series_limit=0)


def partners_ecg_acquisition_year(tm, hd5, dependents={}):
    ecg_dates = _get_ecg_dates(tm, hd5)
    dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
    tensor = np.zeros(shape, dtype=int)
    for i, ecg_date in enumerate(ecg_dates):
        path = _make_hd5_path(tm, ecg_date, 'acquisitiondate')
        try:
            acquisition = decompress_data(data_compressed=hd5[path][()], dtype='str')
            tensor[i] = _partners_str2date(acquisition).year
        except KeyError:
            pass
    return tensor


TMAPS['partners_ecg_acquisition_year'] = TensorMap('partners_ecg_acquisition_year', path_prefix=PARTNERS_PREFIX, loss='logcosh',  tensor_from_file=partners_ecg_acquisition_year, shape=(None, 1), time_series_limit=0)


def partners_bmi(tm, hd5, dependents={}):
    ecg_dates = _get_ecg_dates(tm, hd5)
    dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
    tensor = np.zeros(shape, dtype=float)
    for i, ecg_date in enumerate(ecg_dates):
        path = lambda key: _make_hd5_path(tm, ecg_date, key)
        try:
            weight_lbs = decompress_data(data_compressed=hd5[path('weightlbs')][()], dtype='str')
            weight_kg = 0.453592 * float(weight_lbs)
            height_in = decompress_data(data_compressed=hd5[path('heightin')][()], dtype='str')
            height_m = 0.0254 * float(height_in)
            bmi = weight_kg / (height_m*height_m)
            logging.info(f' Height was {height_in} weight: {weight_lbs} bmi is {bmi}')
            tensor[i] = bmi
        except KeyError:
            pass
    return tensor


TMAPS['partners_ecg_bmi'] = TensorMap('partners_ecg_bmi', path_prefix=PARTNERS_PREFIX, channel_map={'bmi': 0}, tensor_from_file=partners_bmi, time_series_limit=0)


def partners_channel_string(hd5_key, race_synonyms={}, unspecified_key=None):
    def tensor_from_string(tm, hd5, dependents={}):
        ecg_dates = _get_ecg_dates(tm, hd5)
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        tensor = np.zeros(shape, dtype=np.float32)
        for i, ecg_date in enumerate(ecg_dates):
            path = _make_hd5_path(tm, ecg_date,hd5_key)
            found = False
            try:
                hd5_string = decompress_data(data_compressed=hd5[path][()], dtype='str')
                for key in tm.channel_map:
                    slices = (i, tm.channel_map[key]) if dynamic else (tm.channel_map[key],)
                    if hd5_string.lower() == key.lower():
                        tensor[slices] = 1.0
                        found = True
                        break
                    if key in race_synonyms:
                        for synonym in race_synonyms[key]:
                            if hd5_string.lower() == synonym.lower():
                                tensor[slices] = 1.0
                                found = True
                            if found: break
                        if found: break
            except KeyError:
                pass
            if not found:
                if unspecified_key is None:
                    # TODO Do we want to try to continue to get tensors for other ECGs in HD5?
                    raise ValueError(f'No channel keys found in {hd5_string} for {tm.name} with channel map {tm.channel_map}.')
                slices = (i, tm.channel_map[unspecified_key]) if dynamic else (tm.channel_map[unspecified_key],)
                tensor[slices] = 1.0
        return tensor
    return tensor_from_string


race_synonyms = {'asian': ['oriental'], 'hispanic': ['latino'], 'white': ['caucasian']}
TMAPS['partners_ecg_race'] = TensorMap(
    'partners_ecg_race', interpretation=Interpretation.CATEGORICAL, path_prefix=PARTNERS_PREFIX, channel_map={'asian': 0, 'black': 1, 'hispanic': 2, 'white': 3, 'unknown': 4},
    tensor_from_file=partners_channel_string('race', race_synonyms), time_series_limit=0,
)


def _partners_adult(hd5_key, minimum_age=18):
    def tensor_from_string(tm, hd5, dependents={}):
        ecg_dates = _get_ecg_dates(tm, hd5)
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        tensor = np.zeros(shape, dtype=np.float32)
        for i, ecg_date in enumerate(ecg_dates):
            path = lambda key: _make_hd5_path(tm, ecg_date, key)
            birthday = decompress_data(data_compressed=hd5[path('dateofbirth')][()], dtype='str')
            acquisition = decompress_data(data_compressed=hd5[path('acquisitiondate')][()], dtype='str')
            delta = _partners_str2date(acquisition) - _partners_str2date(birthday)
            years = delta.days / YEAR_DAYS
            if years < minimum_age:
                raise ValueError(f'ECG taken on patient below age cutoff.')
            hd5_string = decompress_data(data_compressed=hd5[path(hd5_key)][()], dtype=hd5[path(hd5_key)].attrs['dtype'])
            found = False
            for key in tm.channel_map:
                if hd5_string.lower() == key.lower():
                    slices = (i, tm.channel_map[key]) if dynamic else (tm.channel_map[key],)
                    tensor[slices] = 1.0
                    found = True
                    break
            if not found:
                # TODO Do we want to try to continue to get tensors for other ECGs in HD5?
                raise ValueError(f'No channel keys found in {hd5_string} for {tm.name} with channel map {tm.channel_map}.')
        return tensor
    return tensor_from_string


TMAPS['partners_adult_sex'] = TensorMap(
    'adult_sex', interpretation=Interpretation.CATEGORICAL, path_prefix=PARTNERS_PREFIX, channel_map={'female': 0, 'male': 1},
    tensor_from_file=_partners_adult('gender'), time_series_limit=0,
)


def voltage_zeros(tm, hd5, dependents={}):
    ecg_dates = _get_ecg_dates(tm, hd5)
    dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
    tensor = np.zeros(shape, dtype=np.float32)
    for i, ecg_date in enumerate(ecg_dates):
        for cm in tm.channel_map:
            path = _make_hd5_path(tm, ecg_date, cm)
            voltage = decompress_data(data_compressed=hd5[path][()], dtype=hd5[path].attrs['dtype'])
            slices = (i, tm.channel_map[cm]) if dynamic else (tm.channel_map[cm],)
            tensor[slices] = np.count_nonzero(voltage == 0)
    return tensor


TMAPS["voltage_zeros"] = TensorMap(
    "voltage_zeros",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=voltage_zeros,
    shape=(None, 12),
    channel_map=ECG_REST_AMP_LEADS,
    time_series_limit=0,
)


def v6_zeros_validator(tm: TensorMap, tensor: np.ndarray, hd5: h5py.File):
    voltage = decompress_data(data_compressed=hd5['V6'][()], dtype=hd5['V6'].attrs['dtype'])
    if np.count_nonzero(voltage == 0) > 10:
        raise ValueError(f'TensorMap {tm.name} has too many zeros in V6.')


def build_partners_time_series_tensor_maps(
        needed_tensor_maps: List[str],
        time_series_limit: int = 1,
) -> Dict[str, TensorMap]:
    name2tensormap: Dict[str:TensorMap] = {}

    for needed_name in needed_tensor_maps:
        if needed_name.endswith('_newest'):
            base_split = '_newest'
            time_series_order = TimeSeriesOrder.NEWEST
        elif needed_name.endswith('_oldest'):
            base_split = '_oldest'
            time_series_order = TimeSeriesOrder.OLDEST
        elif needed_name.endswith('_random'):
            base_split = '_random'
            time_series_order = TimeSeriesOrder.RANDOM
        else:
            continue

        base_name = needed_name.split(base_split)[0]
        if base_name not in TMAPS:
            continue

        time_tmap = copy.deepcopy(TMAPS[base_name])
        time_tmap.name = needed_name
        time_tmap.shape = time_tmap.shape[1:]
        time_tmap.time_series_limit = time_series_limit
        time_tmap.time_series_order = time_series_order
        time_tmap.metrics = None
        time_tmap.infer_metrics()

        name2tensormap[needed_name] = time_tmap
    return name2tensormap


# Date formatting
def _partners_str2date(d) -> datetime.date:
    return datetime.datetime.strptime(d, PARTNERS_DATE_FORMAT).date()


def _loyalty_str2date(date_string: str) -> datetime.date:
    return str2date(date_string.split(' ')[0])


def _cardiac_surgery_str2date(input_date: str, date_format: str = CARDIAC_SURGERY_DATE_FORMAT) -> datetime.datetime:
    return datetime.datetime.strptime(input_date, date_format)


def build_incidence_tensor_from_file(
    file_name: str, patient_column: str = 'Mrn', birth_column: str = 'birth_date',
    diagnosis_column: str = 'first_stroke', start_column: str = 'start_fu',
    delimiter: str = ',', incidence_only: bool = False, check_birthday: bool = True,
) -> Callable:
    """Build a tensor_from_file function for future (and prior) diagnoses given a TSV of patients and diagnosis dates.

    The tensor_from_file function returned here should be used
    with CATEGORICAL TensorMaps to classify patients by disease state.

    :param file_name: CSV or TSV file with header of patient IDs (MRNs) dates of enrollment and dates of diagnosis
    :param patient_column: The header name of the column of patient ids
    :param diagnosis_column: The header name of the column of disease diagnosis dates
    :param start_column: The header name of the column of enrollment dates
    :param delimiter: The delimiter separating columns of the TSV or CSV
    :param incidence_only: Flag to skip patients whose diagnosis date is prior to acquisition date of input data
    :return: The tensor_from_file function to provide to TensorMap constructors
    """
    error = None
    try:
        with open(file_name, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter=delimiter)
            header = next(reader)
            patient_index = header.index(patient_column)
            birth_index = header.index(birth_column)
            start_index = header.index(start_column)
            date_index = header.index(diagnosis_column)
            date_table = {}
            birth_table = {}
            patient_table = {}
            for row in reader:
                try:
                    patient_key = int(row[patient_index])
                    patient_table[patient_key] = _loyalty_str2date(row[start_index])
                    birth_table[patient_key] = _loyalty_str2date(row[birth_index])
                    if row[date_index] == '' or row[date_index] == 'NULL':
                        continue
                    date_table[patient_key] = _loyalty_str2date(row[date_index])
                    if len(patient_table) % 2000 == 0:
                        logging.debug(f'Processed: {len(patient_table)} patient rows.')
                except ValueError as e:
                    logging.warning(f'val err {e}')
            logging.info(f'Done processing {diagnosis_column} Got {len(patient_table)} patient rows and {len(date_table)} events.')
    except FileNotFoundError as e:
        error = e

    def tensor_from_file(tm: TensorMap, hd5: h5py.File, dependents=None):
        if error:
            raise error

        ecg_dates = _get_ecg_dates(tm, hd5)
        if len(ecg_dates) > 1:
            raise NotImplementedError('Diagnosis models for multiple ECGs are not implemented.')
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        categorical_data = np.zeros(shape, dtype=np.float32)
        for i, ecg_date in enumerate(ecg_dates):
            path = lambda key: _make_hd5_path(tm, ecg_date, key)
            mrn = _hd5_filename_to_mrn_int(hd5.filename)
            mrn_int = int(mrn)

            if mrn_int not in patient_table:
                raise KeyError(f'{tm.name} mrn not in incidence csv')

            if check_birthday:
                birth_date = _partners_str2date(decompress_data(data_compressed=hd5[path('dateofbirth')][()], dtype=hd5[path('dateofbirth')].attrs['dtype']))
                if birth_date != birth_table[mrn_int]:
                    raise ValueError(f'Birth dates do not match! CSV had {birth_table[patient_key]} but HD5 has {birth_date}')

            assess_date = _partners_str2date(decompress_data(data_compressed=hd5[path('acquisitiondate')][()], dtype=hd5[path('acquisitiondate')].attrs['dtype']))
            if assess_date < patient_table[mrn_int]:
                raise ValueError(f'{tm.name} Assessed earlier than enrollment')
            if mrn_int not in date_table:
                index = 0
            else:
                disease_date = date_table[mrn_int]

                if incidence_only and disease_date < assess_date:
                    raise ValueError(f'{tm.name} is skipping prevalent cases.')
                elif incidence_only and disease_date >= assess_date:
                    index = 1
                else:
                    index = 1 if disease_date < assess_date else 2
                logging.debug(f'mrn: {mrn_int}  Got disease_date: {disease_date} assess  {assess_date} index  {index}.')
            slices = (i, index) if dynamic else (index,)
            categorical_data[slices] = 1.0
        return categorical_data
    return tensor_from_file


def _diagnosis_channels(disease: str, incidence_only: bool = False):
    if incidence_only:
        return {f'no_{disease}': 0,  f'future_{disease}': 1}
    return {f'no_{disease}': 0, f'prior_{disease}': 1, f'future_{disease}': 2}


def _outcome_channels(outcome: str):
    return {f'no_{outcome}': 0,  f'{outcome}': 1}


def loyalty_time_to_event(
    file_name: str, incidence_only: bool = False, patient_column: str = 'Mrn',
    follow_up_start_column: str = 'start_fu', follow_up_total_column: str = 'total_fu',
    diagnosis_column: str = 'first_stroke', delimiter: str = ',',
):
    """Build a tensor_from_file function for modeling relative time to event of diagnoses given a TSV of patients and dates.

    The tensor_from_file function returned here should be used
    with TIME_TO_EVENT TensorMaps to model relative time free from a diagnosis for a given disease.

    :param file_name: CSV or TSV file with header of patient IDs (MRNs) dates of enrollment and dates of diagnosis
    :param incidence_only: Flag to skip patients whose diagnosis date is prior to acquisition date of input data
    :param patient_column: The header name of the column of patient ids
    :param follow_up_start_column: The header name of the column of enrollment dates
    :param follow_up_total_column: The header name of the column with total enrollment time (in years)
    :param diagnosis_column: The header name of the column of disease diagnosis dates
    :param delimiter: The delimiter separating columns of the TSV or CSV
    :return: The tensor_from_file function to provide to TensorMap constructors
    """
    error = None
    disease_dicts = defaultdict(dict)
    try:
        with open(file_name, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter=delimiter)
            header = next(reader)
            follow_up_start_index = header.index(follow_up_start_column)
            follow_up_total_index = header.index(follow_up_total_column)
            patient_index = header.index(patient_column)
            date_index = header.index(diagnosis_column)
            for row in reader:
                try:
                    patient_key = int(row[patient_index])
                    disease_dicts['follow_up_start'][patient_key] = _loyalty_str2date(row[follow_up_start_index])
                    disease_dicts['follow_up_total'][patient_key] = float(row[follow_up_total_index])
                    if row[date_index] == '' or row[date_index] == 'NULL':
                        continue
                    disease_dicts['diagnosis_dates'][patient_key] = _loyalty_str2date(row[date_index])
                    if len(disease_dicts['follow_up_start']) % 2000 == 0:
                        logging.debug(f"Processed: {len(disease_dicts['follow_up_start'])} patient rows.")
                except ValueError as e:
                    logging.warning(f'val err {e}')
            logging.info(f"Done processing {diagnosis_column} Got {len(disease_dicts['follow_up_start'])} patient rows and {len(disease_dicts['diagnosis_dates'])} events.")
    except FileNotFoundError as e:
        error = e

    def _cox_tensor_from_file(tm: TensorMap, hd5: h5py.File, dependents=None):
        if error:
            raise error

        ecg_dates = _get_ecg_dates(tm, hd5)
        if len(ecg_dates) > 1:
            raise NotImplementedError('Cox hazard models for multiple ECGs are not implemented.')
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        tensor = np.zeros(tm.shape, dtype=np.float32)
        for i, ecg_date in enumerate(ecg_dates):
            patient_key_from_ecg = _hd5_filename_to_mrn_int(hd5.filename)
            if patient_key_from_ecg not in disease_dicts['follow_up_start']:
                raise KeyError(f'{tm.name} mrn not in incidence csv')

            path = _make_hd5_path(tm, ecg_date, 'acquisitiondate')
            assess_date = _partners_str2date(decompress_data(data_compressed=hd5[path][()], dtype=hd5[path].attrs['dtype']))
            if assess_date < disease_dicts['follow_up_start'][patient_key_from_ecg]:
                raise ValueError(f'Assessed earlier than enrollment.')

            if patient_key_from_ecg not in disease_dicts['diagnosis_dates']:
                has_disease = 0
                censor_date = disease_dicts['follow_up_start'][patient_key_from_ecg] + datetime.timedelta(
                    days=YEAR_DAYS * disease_dicts['follow_up_total'][patient_key_from_ecg],
                )
            else:
                has_disease = 1
                censor_date = disease_dicts['diagnosis_dates'][patient_key_from_ecg]

            if incidence_only and censor_date <= assess_date and has_disease:
                raise ValueError(f'{tm.name} only considers incident diagnoses')

            tensor[(i, 0) if dynamic else 0] = has_disease
            tensor[(i, 1) if dynamic else 1] = (censor_date - assess_date).days
        return tensor
    return _cox_tensor_from_file


def _survival_from_file(
    day_window: int, file_name: str, incidence_only: bool = False, patient_column: str = 'Mrn',
    follow_up_start_column: str = 'start_fu', follow_up_total_column: str = 'total_fu',
    diagnosis_column: str = 'first_stroke', delimiter: str = ',',
) -> Callable:
    """Build a tensor_from_file function for modeling survival curves of diagnoses given a TSV of patients and dates.

    The tensor_from_file function returned here should be used
    with SURVIVAL_CURVE TensorMaps to model survival curves of patients for a given disease.

    :param day_window: Total number of days of follow up the length of the survival curves to learn.
    :param file_name: CSV or TSV file with header of patient IDs (MRNs) dates of enrollment and dates of diagnosis
    :param incidence_only: Flag to skip patients whose diagnosis date is prior to acquisition date of input data
    :param patient_column: The header name of the column of patient ids
    :param follow_up_start_column: The header name of the column of enrollment dates
    :param follow_up_total_column: The header name of the column with total enrollment time (in years)
    :param diagnosis_column: The header name of the column of disease diagnosis dates
    :param delimiter: The delimiter separating columns of the TSV or CSV
    :return: The tensor_from_file function to provide to TensorMap constructors
    """
    error = None
    disease_dicts = defaultdict(dict)
    try:
        with open(file_name, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter=delimiter)
            header = next(reader)
            follow_up_start_index = header.index(follow_up_start_column)
            follow_up_total_index = header.index(follow_up_total_column)
            patient_index = header.index(patient_column)
            date_index = header.index(diagnosis_column)
            for row in reader:
                try:
                    patient_key = int(row[patient_index])
                    disease_dicts['follow_up_start'][patient_key] = _loyalty_str2date(row[follow_up_start_index])
                    disease_dicts['follow_up_total'][patient_key] = float(row[follow_up_total_index])
                    if row[date_index] == '' or row[date_index] == 'NULL':
                        continue
                    disease_dicts['diagnosis_dates'][patient_key] = _loyalty_str2date(row[date_index])
                    if len(disease_dicts['follow_up_start']) % 2000 == 0:
                        logging.debug(f"Processed: {len(disease_dicts['follow_up_start'])} patient rows.")
                except ValueError as e:
                    logging.warning(f'val err {e}')
            logging.info(f"Done processing {diagnosis_column} Got {len(disease_dicts['follow_up_start'])} patient rows and {len(disease_dicts['diagnosis_dates'])} events.")
    except FileNotFoundError as e:
        error = e

    def tensor_from_file(tm: TensorMap, hd5: h5py.File, dependents=None):
        if error:
            raise error

        ecg_dates = _get_ecg_dates(tm, hd5)
        if len(ecg_dates) > 1:
            raise NotImplementedError('Survival curve models for multiple ECGs are not implemented.')
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        survival_then_censor = np.zeros(shape, dtype=np.float32)
        for ed, ecg_date in enumerate(ecg_dates):
            patient_key_from_ecg = _hd5_filename_to_mrn_int(hd5.filename)
            if patient_key_from_ecg not in disease_dicts['follow_up_start']:
                raise KeyError(f'{tm.name} mrn not in incidence csv')

            path = _make_hd5_path(tm, ecg_date, 'acquisitiondate')
            assess_date = _partners_str2date(decompress_data(data_compressed=hd5[path][()], dtype=hd5[path].attrs['dtype']))
            if assess_date < disease_dicts['follow_up_start'][patient_key_from_ecg]:
                raise ValueError(f'Assessed earlier than enrollment.')

            if patient_key_from_ecg not in disease_dicts['diagnosis_dates']:
                has_disease = 0
                censor_date = disease_dicts['follow_up_start'][patient_key_from_ecg] + datetime.timedelta(days=YEAR_DAYS*disease_dicts['follow_up_total'][patient_key_from_ecg])
            else:
                has_disease = 1
                censor_date = disease_dicts['diagnosis_dates'][patient_key_from_ecg]

            intervals = int(shape[1] if dynamic else shape[0] / 2)
            days_per_interval = day_window / intervals

            for i, day_delta in enumerate(np.arange(0, day_window, days_per_interval)):
                cur_date = assess_date + datetime.timedelta(days=day_delta)
                survival_then_censor[(ed, i) if dynamic else i] = float(cur_date < censor_date)
                survival_then_censor[(ed, intervals+i) if dynamic else intervals+i] = has_disease * float(censor_date <= cur_date < censor_date + datetime.timedelta(days=days_per_interval))
                if i == 0 and censor_date <= cur_date:  # Handle prevalent diseases
                    survival_then_censor[(ed, intervals) if dynamic else intervals] = has_disease
                    if has_disease and incidence_only:
                        raise ValueError(f'{tm.name} is skipping prevalent cases.')
            logging.debug(
                f"Got survival disease {has_disease}, censor: {censor_date}, assess {assess_date}, fu start {disease_dicts['follow_up_start'][patient_key_from_ecg]} "
                f"fu total {disease_dicts['follow_up_total'][patient_key_from_ecg]} tensor:{(survival_then_censor[ed] if dynamic else survival_then_censor)[:4]} mid tense: {(survival_then_censor[ed] if dynamic else survival_then_censor)[intervals:intervals+4]} ",
            )
        return survival_then_censor
    return tensor_from_file


def build_partners_tensor_maps(needed_tensor_maps: List[str]) -> Dict[str, TensorMap]:
    name2tensormap: Dict[str, TensorMap] = {}
    diagnosis2column = {
        'atrial_fibrillation': 'first_af', 'blood_pressure_medication': 'first_bpmed',
        'coronary_artery_disease': 'first_cad', 'cardiovascular_disease': 'first_cvd',
        'death': 'death_date', 'diabetes_mellitus': 'first_dm', 'heart_failure': 'first_hf',
        'hypertension': 'first_htn', 'left_ventricular_hypertrophy': 'first_lvh',
        'myocardial_infarction': 'first_mi', 'pulmonary_artery_disease': 'first_pad',
        'stroke': 'first_stroke', 'valvular_disease': 'first_valvular_disease',
    }
    logging.info(f'needed name {needed_tensor_maps}')
    for diagnosis in diagnosis2column:
        # Build diagnosis classification TensorMaps
        name = f'diagnosis_{diagnosis}'
        if name in needed_tensor_maps:
            tensor_from_file_fxn = build_incidence_tensor_from_file(INCIDENCE_CSV, diagnosis_column=diagnosis2column[diagnosis])
            name2tensormap[name] = TensorMap(f'{name}_newest', Interpretation.CATEGORICAL, path_prefix=PARTNERS_PREFIX, channel_map=_diagnosis_channels(diagnosis), tensor_from_file=tensor_from_file_fxn)
        name = f'incident_diagnosis_{diagnosis}'
        if name in needed_tensor_maps:
            tensor_from_file_fxn = build_incidence_tensor_from_file(INCIDENCE_CSV, diagnosis_column=diagnosis2column[diagnosis], incidence_only=True)
            name2tensormap[name] = TensorMap(f'{name}_newest', Interpretation.CATEGORICAL, path_prefix=PARTNERS_PREFIX, channel_map=_diagnosis_channels(diagnosis, incidence_only=True), tensor_from_file=tensor_from_file_fxn)

        # Build time to event TensorMaps
        name = f'cox_{diagnosis}'
        if name in needed_tensor_maps:
            tff = loyalty_time_to_event(INCIDENCE_CSV, diagnosis_column=diagnosis2column[diagnosis])
            name2tensormap[name] = TensorMap(f'{name}_newest', Interpretation.TIME_TO_EVENT, path_prefix=PARTNERS_PREFIX, tensor_from_file=tff)
        name = f'incident_cox_{diagnosis}'
        if name in needed_tensor_maps:
            tff = loyalty_time_to_event(INCIDENCE_CSV, diagnosis_column=diagnosis2column[diagnosis], incidence_only=True)
            name2tensormap[name] = TensorMap(f'{name}_newest', Interpretation.TIME_TO_EVENT, path_prefix=PARTNERS_PREFIX, tensor_from_file=tff)

        # Build survival curve TensorMaps
        for needed_name in needed_tensor_maps:
            if 'survival' not in needed_name:
                continue
            potential_day_string = needed_name.split('_')[-1]
            try:
                days_window = int(potential_day_string)
            except ValueError:
                days_window = 1825  # Default to 5 years of follow up
            name = f'survival_{diagnosis}'
            if name in needed_name:
                tff = _survival_from_file(days_window, INCIDENCE_CSV, diagnosis_column=diagnosis2column[diagnosis])
                name2tensormap[needed_name] = TensorMap(f'{needed_name}_newest', Interpretation.SURVIVAL_CURVE, path_prefix=PARTNERS_PREFIX, shape=(50,), days_window=days_window, tensor_from_file=tff)
            name = f'incident_survival_{diagnosis}'
            if name in needed_name:
                tff = _survival_from_file(days_window, INCIDENCE_CSV, diagnosis_column=diagnosis2column[diagnosis], incidence_only=True)
                name2tensormap[needed_name] = TensorMap(f'{needed_name}_newest', Interpretation.SURVIVAL_CURVE, path_prefix=PARTNERS_PREFIX, shape=(50,), days_window=days_window, tensor_from_file=tff)
    logging.info(f'return names {list(name2tensormap.keys())}')
    return name2tensormap


def build_cardiac_surgery_dict(
    filename: str = CARDIAC_SURGERY_OUTCOMES_CSV,
    patient_column: str = 'medrecn',
    date_column: str = 'surgdt',
    additional_columns: List[str] = [],
) -> Dict[int, Dict[str, Union[int, str]]]:
    keys = [date_column] + additional_columns
    cardiac_surgery_dict = {}
    df = pd.read_csv(
        filename,
        low_memory=False,
        usecols=[patient_column]+keys,
    ).sort_values(by=[patient_column, date_column])
    # sort dataframe such that newest surgery per patient appears later and is used in lookup table
    for row in df.itertuples():
        patient_key = getattr(row, patient_column)
        cardiac_surgery_dict[patient_key] = {key: getattr(row, key) for key in keys}
    return cardiac_surgery_dict


def build_date_interval_lookup(
    cardiac_surgery_dict: Dict[int, Dict[str, Union[int, str]]],
    start_column: str = 'surgdt',
    start_offset: int = -30,
    end_column: str = 'surgdt',
    end_offset: int = 0,
) -> Dict[int, Tuple[str, str]]:
    date_interval_lookup = {}
    for mrn in cardiac_surgery_dict:
        start_date = (_cardiac_surgery_str2date(cardiac_surgery_dict[mrn][start_column], PARTNERS_DATETIME_FORMAT.replace('T', ' ')) + datetime.timedelta(days=start_offset)).strftime(PARTNERS_DATETIME_FORMAT)
        end_date = (_cardiac_surgery_str2date(cardiac_surgery_dict[mrn][end_column], PARTNERS_DATETIME_FORMAT.replace('T', ' ')) + datetime.timedelta(days=end_offset)).strftime(PARTNERS_DATETIME_FORMAT)
        date_interval_lookup[mrn] = (start_date, end_date)
    return date_interval_lookup


def make_cardiac_surgery_outcome_tensor_from_file(
    cardiac_surgery_dict: Dict[int, Dict[str, Union[int, str]]],
    outcome_column: str,
) -> Callable:
    def tensor_from_file(tm: TensorMap, hd5: h5py.File, dependents: Dict = {}):
        mrn = _hd5_filename_to_mrn_int(hd5.filename)
        tensor = np.zeros(tm.shape, dtype=np.float32)
        outcome = cardiac_surgery_dict[mrn][outcome_column]

        if type(outcome) is float and not outcome.is_integer():
            raise ValueError(f'Cardiac Surgery categorical outcome {tm.name} ({outcome_column}) got non-discrete value: {outcome}')

        # ensure binary outcome
        if outcome != 0 and outcome != 1:
            raise ValueError(f'Cardiac Surgery categorical outcome {tm.name} ({outcome_column}) got non-binary value: {outcome}')

        tensor[outcome] = 1
        return tensor
    return tensor_from_file


def build_cardiac_surgery_tensor_maps(
    needed_tensor_maps: List[str],
) -> Dict[str, TensorMap]:
    name2tensormap: Dict[str, TensorMap] = {}
    outcome2column = {
        "sts_death": "mtopd",
        "sts_stroke": "cnstrokp",
        "sts_renal_failure": "crenfail",
        "sts_prolonged_ventilation": "cpvntlng",
        "sts_dsw_infection": "deepsterninf",
        "sts_reoperation": "reop",
        "sts_any_morbidity": "anymorbidity",
        "sts_long_stay": "llos",
    }

    cardiac_surgery_dict = None
    date_interval_lookup = None
    for needed_name in needed_tensor_maps:
        if needed_name in outcome2column:
            if cardiac_surgery_dict is None:
                cardiac_surgery_dict = build_cardiac_surgery_dict(additional_columns=[column for outcome, column in outcome2column.items() if outcome in needed_tensor_maps])
            channel_map = _outcome_channels(needed_name)
            sts_tmap = TensorMap(
                needed_name,
                Interpretation.CATEGORICAL,
                path_prefix=PARTNERS_PREFIX,
                tensor_from_file=make_cardiac_surgery_outcome_tensor_from_file(cardiac_surgery_dict, outcome2column[needed_name]),
                channel_map=channel_map,
                validator=validator_not_all_zero,
            )
        else:
            if not needed_name.endswith('_sts'):
                continue

            base_name = needed_name.split('_sts')[0]
            if base_name not in TMAPS:
                TMAPS.update(build_partners_time_series_tensor_maps([base_name]))
                if base_name not in TMAPS:
                    continue

            if cardiac_surgery_dict is None:
                cardiac_surgery_dict = build_cardiac_surgery_dict(additional_columns=[column for outcome, column in outcome2column.items() if outcome in needed_tensor_maps])
            if date_interval_lookup is None:
                date_interval_lookup = build_date_interval_lookup(cardiac_surgery_dict)
            sts_tmap = copy.deepcopy(TMAPS[base_name])
            sts_tmap.name = needed_name
            sts_tmap.time_series_lookup = date_interval_lookup

        name2tensormap[needed_name] = sts_tmap

    return name2tensormap


# Measurement matrix TMAPS -- indices from MUSE XML dev manual, page 49 and following
measurement_matrix_leads = {
    'I': 0, 'II': 1, 'V1': 2, 'V2': 3, 'V3': 4, 'V4':5, 'V5': 6, 'V6': 7, 'III': 8, 'aVR': 9, 'aVL': 10, 'aVF': 11
}
measurement_matrix_global_measures = {
    'pon': 1,       # P-wave onset in median beat (in samples)
    'poff': 2,      # P-wave offset in median beat
    'qon': 3,       # Q-Onset in median beat
    'qoff': 4,      # Q-Offset in median beat
    'ton': 5,       # T-Onset in median beat
    'toff': 6,      # T-Offset in median beat
    'nqrs': 7,      # Number of QRS Complexes
    'qrsdur': 8,    # QRS Duration
    'qt': 9,        # QT Interval
    'qtc': 10,      # QT Corrected
    'print': 11,    # PR Interval
    'vrate': 12,    # Ventricular Rate
    'avgrr': 13,    # Average R-R Interval
}
measurement_matrix_lead_measures = {
    'pona': 1,      # P Wave amplitude at P-onset
    'pamp': 2,      # P wave amplitude
    'pdur': 3,      # P wave duration
    'bmpar': 4,     # P wave area
    'bmpi': 5,      # P wave intrinsicoid (time from P onset to peak of P)
    'ppamp': 6,     # P Prime amplitude
    'ppdur': 7,     # P Prime duration
    'bmppar': 8,    # P Prime area
    'bmppi': 9,     # P Prime intrinsicoid (time from P onset to peak of P')
    'qamp': 10,     # Q wave amplitude
    'qdur': 11,     # Q wave duration
    'bmqar': 12,    # Q wave area
    'bmqi': 13,     # Q intrinsicoid (time from Q onset to peak of Q)
    'ramp': 14,     # R amplitude
    'rdur': 15,     # R duration
    'bmrar': 16,    # R wave area
    'bmri': 17,     # R intrinsicoid (time from R onset to peak of R)
    'samp': 18,     # S amplitude
    'sdur': 19,     # S duration
    'bmsar': 20,    # S wave area
    'bmsi': 21,     # S intrinsicoid (time from Q onset to peak of S)
    'rpamp': 22,    # R Prime amplitude
    'rpdur': 23,    # R Prime duration
    'bmrpar': 24,   # R Prime wave area
    'bmrpi': 25,    # R Prime intrinsicoid (time from Q onset to peak of R Prime)
    'spamp': 26,    # S Prime Amplitude
    'spdur': 27,    # S Prime Duration
    'bmspar': 28,   # S Prime wave area
    'bmspi': 29,    # S intriniscoid (time from Q onset to peak of S prime)
    'stj': 30,      # STJ point, End of QRS Point Amplitude
    'stm': 31,      # STM point, Middle of the ST Segment Amplitude
    'ste': 32,      # STE point, End of ST Segment Amplitude
    'mxsta': 33,    # Maximum of STJ, STM, STE Amplitudes
    'mnsta': 34,    # Minimum of STJ and STM Amplitudes
    'spta': 35,     # Special T-Wave amplitude
    'qrsa': 36,     # Total QRS area
    'qrsdef': 37,   # QRS Deflection
    'maxra': 38,    # Maximum R Amplitude (R or R Prime)
    'maxsa': 39,    # Maximum S Amplitude (S or S Prime)
    'tamp': 40,     # T amplitude
    'tdur': 41,     # T duration
    'bmtar': 42,    # T wave area
    'bmti': 43,     # T intriniscoid (time from STE to peak of T)
    'tpamp': 44,    # T Prime amplitude
    'tpdur': 45,    # T Prime duration
    'bmtpar': 46,   # T Prime area
    'bmtpi': 47,    # T Prime intriniscoid (time from STE to peak of T)
    'tend': 48,     # T Amplitude at T offset
    'parea': 49,    # P wave area, includes P and P Prime
    'qrsar': 50,    # QRS area
    'tarea': 51,    # T wave area, include T and T Prime
    'qrsint': 52    # QRS intriniscoid (see following)
}


def _get_measurement_matrix_entry(matrix: np.ndarray, key: str, lead: str = None):
    # First 18 words of measurement matrix are for global measurements, then each lead has 53*2 words
    lead_start = 18
    lead_words = 53 * 2
    if lead is None:
        idx = measurement_matrix_global_measures[key]
    else:
        idx = lead_start + measurement_matrix_leads[lead] * lead_words + (measurement_matrix_lead_measures[key]-1)*2+1
    return matrix[idx]


def make_measurement_matrix_from_file(key: str, lead: str = None):
    def measurement_matrix_from_file(tm: TensorMap, hd5: h5py.File, dependents: Dict = {}):        
        ecg_dates = _get_ecg_dates(tm, hd5)
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        tensor = np.zeros(shape, dtype=float)
        for i, ecg_date in enumerate(ecg_dates):
            path = _make_hd5_path(tm, ecg_date, 'measurementmatrix')
            matrix = decompress_data(data_compressed=hd5[path][()], dtype=hd5[path].attrs['dtype'])
            tensor[i] = _get_measurement_matrix_entry(matrix, key, lead)
        return tensor
    return measurement_matrix_from_file


for measurement in measurement_matrix_global_measures:
    TMAPS[f'partners_ecg_measurementmatrix_{measurement}'] = TensorMap(
        f'partners_ecg_measurementmatrix_{measurement}',
        interpretation=Interpretation.CONTINUOUS,
        shape=(None, 1),
        path_prefix=PARTNERS_PREFIX,
        loss='logcosh',
        time_series_limit=0,
        tensor_from_file=make_measurement_matrix_from_file(measurement)
    )


for lead in measurement_matrix_leads:
    for measurement in measurement_matrix_lead_measures:
        TMAPS[f'partners_ecg_measurementmatrix_{lead}_{measurement}'] = TensorMap(
              f'partners_ecg_measurementmatrix_{lead}_{measurement}',
              interpretation=Interpretation.CONTINUOUS,
              shape=(None, 1),
              path_prefix=PARTNERS_PREFIX,
              loss='logcosh',
              time_series_limit=0,
              tensor_from_file=make_measurement_matrix_from_file(measurement, lead=lead)
        )


def ecg_lvh_from_file(tm, hd5, dependents={}):
    # Lead order seems constant and standard throughout, but we could eventually tensorize it from XML
    avl_min = 1100.0
    sl_min = 3500.0
    cornell_female_min = 2000.0
    cornell_male_min = 2800.0
    sleads = ['V1', 'V3']
    rleads = ['aVL', 'V5', 'V6']
    ecg_dates = _get_ecg_dates(tm, hd5)
    dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
    tensor = np.zeros(shape, dtype=float)

    for i, ecg_date in enumerate(ecg_dates):
        path = _make_hd5_path(tm, ecg_date, 'measurementmatrix')
        matrix = decompress_data(data_compressed=hd5[path][()], dtype=hd5[path].attrs['dtype'])
        criteria_sleads = {lead: _get_measurement_matrix_entry(matrix, 'samp', lead) for lead in sleads}
        criteria_rleads = {lead: _get_measurement_matrix_entry(matrix, 'ramp', lead) for lead in rleads}
        sex_path = _make_hd5_path(tm, ecg_date, 'gender')
        is_female = 'female' in decompress_data(data_compressed=hd5[sex_path][()], dtype=hd5[sex_path].attrs['dtype'])
        if 'avl_lvh' in tm.name:
            is_lvh = criteria_rleads['aVL'] > avl_min
        elif 'sokolow_lyon_lvh' in tm.name:
            is_lvh = criteria_sleads['V1'] + np.maximum(criteria_rleads['V5'], criteria_rleads['V6']) > sl_min
        elif 'cornell_lvh' in tm.name:
            is_lvh = criteria_rleads['aVL'] + criteria_sleads['V3']
            if is_female:
                is_lvh = is_lvh > cornell_female_min
            else:
                is_lvh = is_lvh > cornell_male_min
        else:
            raise ValueError(f'{tm.name} criterion for LVH is not accounted for')
        # Following convention from categorical TMAPS, positive has cmap index 1
        index = 1 if is_lvh else 0
        slices = (i, index) if dynamic else (index,)
        tensor[slices] = 1.0
    return tensor


for criterion in ['avl_lvh', 'sokolow_lyon_lvh', 'cornell_lvh']:
    TMAPS[f'partners_ecg_{criterion}'] = TensorMap(
        f'partners_ecg_{criterion}',
        interpretation=Interpretation.CATEGORICAL,
        path_prefix=PARTNERS_PREFIX,
        tensor_from_file=ecg_lvh_from_file,
        channel_map={f'no_{criterion}': 0, criterion: 1},
        shape=(None, 2),
        time_series_limit=0,
    )
