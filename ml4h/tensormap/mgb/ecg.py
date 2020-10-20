import os
import copy
import logging
import datetime
from itertools import product
from collections import defaultdict
from typing import Callable, Dict, List, Tuple, Union

import csv
import h5py
import numpy as np
import pandas as pd

from ml4h.metrics import weighted_crossentropy
from ml4h.normalizer import Standardize, ZeroMeanStd1
from ml4h.TensorMap import TensorMap, str2date, Interpretation, make_range_validator, decompress_data, TimeSeriesOrder
from ml4h.defines import ECG_REST_AMP_LEADS, PARTNERS_DATE_FORMAT, STOP_CHAR, PARTNERS_DATETIME_FORMAT, CARDIAC_SURGERY_DATE_FORMAT

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


# Date formatting
def _partners_str2date(d) -> datetime.date:
    return datetime.datetime.strptime(d, PARTNERS_DATE_FORMAT).date()


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


partners_ecg_voltage_stats = TensorMap(
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


voltage_len = TensorMap(
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


partners_ecg_datetime = TensorMap(
    "partners_ecg_datetime",
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=partners_ecg_datetime,
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
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
                data = decompress_data(data_compressed=hd5[path][()], dtype='str')
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


partners_ecg_read_md = TensorMap(
    "partners_ecg_read_md",
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


partners_ecg_read_pc = TensorMap(
    "partners_ecg_read_pc",
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


partners_ecg_patientid = TensorMap(
    "partners_ecg_patientid",
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="patientid"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
)


def validator_clean_mrn(tm: TensorMap, tensor: np.ndarray, hd5: h5py.File):
    int(tensor)


partners_ecg_patientid_clean = TensorMap(
    "partners_ecg_patientid_clean",
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="patientid_clean"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_clean_mrn,
)


partners_ecg_firstname = TensorMap(
    "partners_ecg_firstname",
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


partners_ecg_lastname = TensorMap(
    "partners_ecg_lastname",
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


partners_ecg_sex = TensorMap(
    "partners_ecg_sex",
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="gender"),
    channel_map={'female': 0, 'male': 1},
    time_series_limit=0,
    validator=validator_not_all_zero,
)

partners_ecg_date = TensorMap(
    "partners_ecg_date",
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="acquisitiondate"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
)


partners_ecg_time = TensorMap(
    "partners_ecg_time",
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="acquisitiontime"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
)


partners_ecg_sitename = TensorMap(
    "partners_ecg_sitename",
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="sitename"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
)


partners_ecg_location = TensorMap(
    "partners_ecg_location",
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="location"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_no_empty,
)


partners_ecg_dob = TensorMap(
    "partners_ecg_dob",
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


partners_ecg_sampling_frequency = TensorMap(
    "partners_ecg_sampling_frequency",
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_sampling_frequency_from_file(),
    channel_map={'_250': 0, '_500': 1, 'other': 2},
    time_series_limit=0,
    validator=validator_not_all_zero,
)


partners_ecg_sampling_frequency_pc = TensorMap(
    "partners_ecg_sampling_frequency_pc",
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="ecgsamplebase_pc", channel_prefix='_'),
    channel_map={'_0': 0, '_250': 1, '_500': 2, 'other': 3},
    time_series_limit=0,
    validator=validator_not_all_zero,
)


partners_ecg_sampling_frequency_md = TensorMap(
    "partners_ecg_sampling_frequency_md",
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="ecgsamplebase_md", channel_prefix='_'),
    channel_map={'_0': 0, '_250': 1, '_500': 2, 'other': 3},
    time_series_limit=0,
    validator=validator_not_all_zero,
)


partners_ecg_sampling_frequency_wv = TensorMap(
    "partners_ecg_sampling_frequency_wv",
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="waveform_samplebase", channel_prefix='_'),
    channel_map={'_0': 0, '_240': 1, '_250': 2, '_500': 3, 'other': 4},
    time_series_limit=0,
    validator=validator_not_all_zero,
)


partners_ecg_sampling_frequency_continuous = TensorMap(
    "partners_ecg_sampling_frequency_continuous",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_sampling_frequency_from_file(),
    time_series_limit=0,
    shape=(None, 1),
    validator=validator_no_negative,
)


partners_ecg_sampling_frequency_pc_continuous = TensorMap(
    "partners_ecg_sampling_frequency_pc_continuous",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="ecgsamplebase_pc", fill=-1),
    time_series_limit=0,
    shape=(None, 1),
    validator=validator_no_negative,
)


partners_ecg_sampling_frequency_md_continuous = TensorMap(
    "partners_ecg_sampling_frequency_md_continuous",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="ecgsamplebase_md", fill=-1),
    time_series_limit=0,
    shape=(None, 1),
    validator=validator_no_negative,
)


partners_ecg_sampling_frequency_wv_continuous = TensorMap(
    "partners_ecg_sampling_frequency_wv_continuous",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="waveform_samplebase", fill=-1),
    time_series_limit=0,
    shape=(None, 1),
    validator=validator_no_negative,
)


partners_ecg_time_resolution = TensorMap(
    "partners_ecg_time_resolution",
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="intervalmeasurementtimeresolution", channel_prefix='_'),
    channel_map={'_25': 0, '_50': 1, '_100': 2, 'other': 3},
    time_series_limit=0,
    validator=validator_not_all_zero,
)


partners_ecg_amplitude_resolution = TensorMap(
    "partners_ecg_amplitude_resolution",
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="intervalmeasurementamplituderesolution", channel_prefix='_'),
    channel_map={'_10': 0, '_20': 1, '_40': 2, 'other': 3},
    time_series_limit=0,
    validator=validator_not_all_zero,
)


partners_ecg_measurement_filter = TensorMap(
    "partners_ecg_measurement_filter",
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="intervalmeasurementfilter", channel_prefix='_'),
    time_series_limit=0,
    channel_map={'_None': 0, '_40': 1, '_80': 2, 'other': 3},
    validator=validator_not_all_zero,
)


partners_ecg_high_pass_filter = TensorMap(
    "partners_ecg_high_pass_filter",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="waveform_highpassfilter", fill=-1),
    time_series_limit=0,
    shape=(None, 1),
    validator=validator_no_negative,
)


partners_ecg_low_pass_filter = TensorMap(
    "partners_ecg_low_pass_filter",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="waveform_lowpassfilter", fill=-1),
    time_series_limit=0,
    shape=(None, 1),
    validator=validator_no_negative,
)


partners_ecg_ac_filter = TensorMap(
    "partners_ecg_ac_filter",
    interpretation=Interpretation.CATEGORICAL,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="waveform_acfilter", channel_prefix='_'),
    time_series_limit=0,
    channel_map={'_None': 0, '_50': 1, '_60': 2, 'other': 3},
    validator=validator_not_all_zero,
)


partners_ecg_rate_pc = TensorMap(
    "partners_ecg_rate_pc",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="ventricularrate_pc"),
    shape=(None, 1),
    normalization=Standardize(mean=59.3, std=10.6),
    time_series_limit=0,
    validator=make_range_validator(10, 200),
)


partners_ecg_rate_md = TensorMap(
    "partners_ecg_rate_md",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="ventricularrate_md"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(10, 200),
)


partners_ecg_qrs_pc = TensorMap(
    "partners_ecg_qrs_pc",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="qrsduration_pc"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(20, 400),
)


partners_ecg_qrs_md = TensorMap(
    "partners_ecg_qrs_md",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="qrsduration_md"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(20, 400),
)


partners_ecg_pr_pc = TensorMap(
    "partners_ecg_pr_pc",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="printerval_pc"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(50, 500),
)


partners_ecg_pr_md = TensorMap(
    "partners_ecg_pr_md",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="printerval_md"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(50, 500),
)


partners_ecg_qt_pc = TensorMap(
    "partners_ecg_qt_pc",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="qtinterval_pc"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(100, 800),
)


partners_ecg_qt_md = TensorMap(
    "partners_ecg_qt_md",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="qtinterval_md"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(100, 800),
)


partners_ecg_qtc_pc = TensorMap(
    "partners_ecg_qtc_pc",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="qtcorrected_pc"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(100, 800),
)


partners_ecg_qtc_md = TensorMap(
    "partners_ecg_qtc_md",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="qtcorrected_md"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(100, 800),
)


partners_ecg_paxis_pc = TensorMap(
    "partners_ecg_paxis_pc",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="paxis_pc", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


partners_ecg_paxis_md = TensorMap(
    "partners_ecg_paxis_md",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="paxis_md", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


partners_ecg_raxis_pc = TensorMap(
    "partners_ecg_raxis_pc",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="raxis_pc", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


partners_ecg_raxis_md = TensorMap(
    "partners_ecg_raxis_md",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="raxis_md", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


partners_ecg_taxis_pc = TensorMap(
    "partners_ecg_taxis_pc",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="taxis_pc", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


partners_ecg_taxis_md = TensorMap(
    "partners_ecg_taxis_md",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="taxis_md", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


partners_ecg_measuredamplitudepeak_r = TensorMap(
    "partners_ecg_measuredamplitudepeak_r",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="measuredamplitudepeak_IE_R", fill=np.nan),
    shape=(None, 12),
    time_series_limit=0,
)


partners_ecg_acquisitiondevice = TensorMap(
    "partners_ecg_acquisitiondevice",
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="acquisitiondevice", fill=999),
    shape=(None, 1),
    time_series_limit=0,
)


partners_ecg_weight_lbs = TensorMap(
    "partners_ecg_weight_lbs",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="weightlbs"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(100, 800),
)


def _partners_ecg_age_from_hd5(tm, hd5, dependents={}):
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


partners_ecg_age = TensorMap('partners_ecg_age', path_prefix=PARTNERS_PREFIX, loss='logcosh', tensor_from_file=_partners_ecg_age_from_hd5, shape=(None, 1), time_series_limit=0)


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


partners_ecg_acquisition_year = TensorMap('partners_ecg_acquisition_year', path_prefix=PARTNERS_PREFIX, loss='logcosh',  tensor_from_file=partners_ecg_acquisition_year, shape=(None, 1), time_series_limit=0)


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


partners_ecg_bmi = TensorMap('partners_ecg_bmi', path_prefix=PARTNERS_PREFIX, channel_map={'bmi': 0}, tensor_from_file=partners_bmi, time_series_limit=0)


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
partners_ecg_race = TensorMap(
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


partners_adult_sex = TensorMap(
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


voltage_zeros = TensorMap(
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


partners_ecg_read_md_clean_supranodal_rhythms = TensorMap('partners_ecg_read_md_clean_supranodal_rhythms', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest', channel_map={'retrograde_atrial_activation': 0, 'narrow_qrs_tachycardia': 1, 'unspecified': 2, 'svt': 3, 'ectopic_atrial_tachycardia': 4, 'ectopic_atrial_rhythm': 5, 'atrial_fibrillation': 6, 'atrial_flutter': 7, 'sinus_rhythm': 8}, tensor_from_file=make_partners_ecg_label(keys='read_md_clean', dict_of_list = {'atrial_fibrillation': ['atrial fib', 'afib', 'atrial fibrillation', 'atrial  fibrillation', 'afibrillation'], 'sinus_rhythm': ['marked sinus arrhythmia', 'sinus arrhythmia', 'normal ecg', 'atrial bigeminal rhythm', 'atrial bigeminy and ventricular bigeminy', 'atrial trigeminy', 'atrialbigeminy', 'normal sinus rhythm', 'normal when compared with ecg of', 'rhythm has reverted to normal', 'rhythm is now clearly sinus', 'sinus bradycardia', 'sinus rhythm', 'sinus rhythm at a rate', 'sinus tachycardia', 'tracing within normal limits', 'tracing is within normal limits', 'tracing within normal limits', 'sinoatrial block', 'sa block', 'sinoatrial block, type ii', 'sa block, type i', 'type i sinoatrial block', 'type i sa block', 'sinoatrial block, type ii', 'sa block, type i', 'type ii sinoatrial block', 'type ii sa block', 'sinus pause', 'sinus arrest', 'sinus slowing', 'with occasional native sinus beats', 'frequent native sinus beats', 'conducted sinus impulses', 'sinus mechanism has replaced', 'rhythm is normal sinus', 'rhythm remains normal sinus'], 'atrial_flutter': ['fibrillation/flutter', 'aflutter', 'atrial flutter', 'probable flutter', 'atrial flutter fixed block', 'atrial flutter variable block', 'atrial flutter unspecified block', 'tachycardia possibly flutter'], 'ectopic_atrial_rhythm': ['unifocal atrial tachycardia', 'unifocal ectopic atrial rhythm', 'unifocal ear', 'multifocal atrial tachycardia', 'multifocal ectopic atrial rhythm', 'multifocal ear', 'wandering atrial tachycardia', 'wandering ectopic atrial rhythm', 'wandering ear', 'wandering atrial pacemaker', 'atrial rhythm', 'ectopic atrial rhythm', 'ear', 'dual atrial foci ', 'ectopic atrial bradycardia', 'multiple atrial foci', 'abnormal p vector', 'nonsinus atrial mechanism', 'p wave axis suggests atrial rather than sinus mechanism', 'low atrial bradycardia', 'multifocal atrial rhythm', 'multifocal atrialrhythm', 'unusual p wave axis', 'low atrial pacer'], 'ectopic_atrial_tachycardia': ['ectopic atrial tachycardia', 'ectopic atrial tachycardia, unspecified', 'unspecified ectopic atrial tachycardia', 'ectopic atrial tachycardia, unifocal', 'unifocal ectopic atrial tachycardia', 'ectopic atrial tachycardia, multifocal', 'multifocal ectopic atrial tachycardia'], 'narrow_qrs_tachycardia': ['narrow qrs tachycardia', 'narrow qrs tachycardia', 'tachycardia narrow qrs', 'narrow complex tachycardia'], 'retrograde_atrial_activation': ['retrograde atrial activation'], 'svt': ['sinus arrhythmia accelerated atrioventricular junctional rhythm', 'supraventricular tachycardia', 'accelerated atrioventricular nodal rhythm', 'accelerated nodal rhythm', 'atrial tachycardia', 'av nodal reentry tachycardia', 'atrioventricular nodal reentry tachycardia', 'avnrt', 'atrioventricular reentrant tachycardia ', 'av reentrant tachycardia ', 'avrt'], 'unspecified': ['junctional tachycardia', 'atrial arrhythmia', 'technically poor tracing ', 'accelerated idioventricular rhythm', 'atrial activity is indistinct', 'rhythm uncertain', 'rhythm unclear', 'uncertain rhythm', 'undetermined rhythm', 'supraventricular rhythm']}))

partners_ecg_read_pc_clean_supranodal_rhythms = TensorMap('partners_ecg_read_pc_clean_supranodal_rhythms', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest', channel_map={'retrograde_atrial_activation': 0, 'narrow_qrs_tachycardia': 1, 'unspecified': 2, 'svt': 3, 'ectopic_atrial_tachycardia': 4, 'ectopic_atrial_rhythm': 5, 'atrial_fibrillation': 6, 'atrial_flutter': 7, 'sinus_rhythm': 8}, tensor_from_file=make_partners_ecg_label(keys='read_pc_clean', dict_of_list = {'atrial_fibrillation': ['atrial fib', 'afib', 'atrial fibrillation', 'atrial  fibrillation', 'afibrillation'], 'sinus_rhythm': ['marked sinus arrhythmia', 'sinus arrhythmia', 'normal ecg', 'atrial bigeminal rhythm', 'atrial bigeminy and ventricular bigeminy', 'atrial trigeminy', 'atrialbigeminy', 'normal sinus rhythm', 'normal when compared with ecg of', 'rhythm has reverted to normal', 'rhythm is now clearly sinus', 'sinus bradycardia', 'sinus rhythm', 'sinus rhythm at a rate', 'sinus tachycardia', 'tracing within normal limits', 'tracing is within normal limits', 'tracing within normal limits', 'sinoatrial block', 'sa block', 'sinoatrial block, type ii', 'sa block, type i', 'type i sinoatrial block', 'type i sa block', 'sinoatrial block, type ii', 'sa block, type i', 'type ii sinoatrial block', 'type ii sa block', 'sinus pause', 'sinus arrest', 'sinus slowing', 'with occasional native sinus beats', 'frequent native sinus beats', 'conducted sinus impulses', 'sinus mechanism has replaced', 'rhythm is normal sinus', 'rhythm remains normal sinus'], 'atrial_flutter': ['fibrillation/flutter', 'aflutter', 'atrial flutter', 'probable flutter', 'atrial flutter fixed block', 'atrial flutter variable block', 'atrial flutter unspecified block', 'tachycardia possibly flutter'], 'ectopic_atrial_rhythm': ['unifocal atrial tachycardia', 'unifocal ectopic atrial rhythm', 'unifocal ear', 'multifocal atrial tachycardia', 'multifocal ectopic atrial rhythm', 'multifocal ear', 'wandering atrial tachycardia', 'wandering ectopic atrial rhythm', 'wandering ear', 'wandering atrial pacemaker', 'atrial rhythm', 'ectopic atrial rhythm', 'ear', 'dual atrial foci ', 'ectopic atrial bradycardia', 'multiple atrial foci', 'abnormal p vector', 'nonsinus atrial mechanism', 'p wave axis suggests atrial rather than sinus mechanism', 'low atrial bradycardia', 'multifocal atrial rhythm', 'multifocal atrialrhythm', 'unusual p wave axis', 'low atrial pacer'], 'ectopic_atrial_tachycardia': ['ectopic atrial tachycardia', 'ectopic atrial tachycardia, unspecified', 'unspecified ectopic atrial tachycardia', 'ectopic atrial tachycardia, unifocal', 'unifocal ectopic atrial tachycardia', 'ectopic atrial tachycardia, multifocal', 'multifocal ectopic atrial tachycardia'], 'narrow_qrs_tachycardia': ['narrow qrs tachycardia', 'narrow qrs tachycardia', 'tachycardia narrow qrs', 'narrow complex tachycardia'], 'retrograde_atrial_activation': ['retrograde atrial activation'], 'svt': ['sinus arrhythmia accelerated atrioventricular junctional rhythm', 'supraventricular tachycardia', 'accelerated atrioventricular nodal rhythm', 'accelerated nodal rhythm', 'atrial tachycardia', 'av nodal reentry tachycardia', 'atrioventricular nodal reentry tachycardia', 'avnrt', 'atrioventricular reentrant tachycardia ', 'av reentrant tachycardia ', 'avrt'], 'unspecified': ['junctional tachycardia', 'atrial arrhythmia', 'technically poor tracing ', 'accelerated idioventricular rhythm', 'atrial activity is indistinct', 'rhythm uncertain', 'rhythm unclear', 'uncertain rhythm', 'undetermined rhythm', 'supraventricular rhythm']}))

partners_ecg_read_md_clean_sinus_rhythm = TensorMap('partners_ecg_read_md_clean_sinus_rhythm', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest', channel_map={'sinus_arrhythmia': 0, 'unspecified': 1}, tensor_from_file=make_partners_ecg_label(keys='read_md_clean', dict_of_list = {'sinus_arrhythmia': ['marked sinus arrhythmia', 'marked sinus arrhythmia', 'sinus arrhythmia', 'sinus arrhythmia']}))

partners_ecg_read_pc_clean_sinus_rhythm = TensorMap('partners_ecg_read_pc_clean_sinus_rhythm', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest', channel_map={'sinus_arrhythmia': 0, 'unspecified': 1}, tensor_from_file=make_partners_ecg_label(keys='read_pc_clean', dict_of_list = {'sinus_arrhythmia': ['marked sinus arrhythmia', 'marked sinus arrhythmia', 'sinus arrhythmia', 'sinus arrhythmia']}))

partners_ecg_read_md_clean_atrial_flutter = TensorMap('partners_ecg_read_md_clean_atrial_flutter', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest', channel_map={'fixed_block': 0, 'variable_block': 1, 'unspecified': 2}, tensor_from_file=make_partners_ecg_label(keys='read_md_clean', dict_of_list = {'unspecified': ['fibrillation/flutter', 'fibrillation/flutter', 'aflutter', 'aflutter', 'atrial flutter', 'atrial flutter', 'probable flutter', 'probable flutter', 'atrial flutter unspecified block', 'atrial flutter unspecified block'], 'fixed_block': ['atrial flutter fixed block', 'atrial flutter fixed block'], 'variable_block': ['atrial flutter variable block', 'atrial flutter variable block']}))

partners_ecg_read_pc_clean_atrial_flutter = TensorMap('partners_ecg_read_pc_clean_atrial_flutter', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest', channel_map={'fixed_block': 0, 'variable_block': 1, 'unspecified': 2}, tensor_from_file=make_partners_ecg_label(keys='read_pc_clean', dict_of_list = {'unspecified': ['fibrillation/flutter', 'fibrillation/flutter', 'aflutter', 'aflutter', 'atrial flutter', 'atrial flutter', 'probable flutter', 'probable flutter', 'atrial flutter unspecified block', 'atrial flutter unspecified block'], 'fixed_block': ['atrial flutter fixed block', 'atrial flutter fixed block'], 'variable_block': ['atrial flutter variable block', 'atrial flutter variable block']}))

partners_ecg_read_md_clean_ectopic_atrial_rhythm = TensorMap('partners_ecg_read_md_clean_ectopic_atrial_rhythm', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest', channel_map={'wandering': 0, 'unifocal': 1, 'multifocal': 2, 'unspecified': 3}, tensor_from_file=make_partners_ecg_label(keys='read_md_clean', dict_of_list = {'unifocal': ['unifocal atrial tachycardia', 'unifocal atrial tachycardia', 'unifocal ectopic atrial rhythm', 'unifocal ectopic atrial rhythm', 'unifocal ear', 'unifocal ear', 'unusual p wave axis', 'unusual p wave axis', 'low atrial pacer', 'low atrial pacer'], 'multifocal': ['multifocal atrial tachycardia', 'multifocal atrial tachycardia', 'multifocal ectopic atrial rhythm', 'multifocal ectopic atrial rhythm', 'multifocal ear', 'multifocal ear', 'dual atrial foci ', 'dual atrial foci ', 'multiple atrial foci', 'multiple atrial foci', 'multifocal atrial rhythm', 'multifocal atrial rhythm', 'multifocal atrialrhythm', 'multifocal atrialrhythm'], 'wandering': ['wandering atrial tachycardia', 'wandering atrial tachycardia', 'wandering ectopic atrial rhythm', 'wandering ectopic atrial rhythm', 'wandering ear', 'wandering ear', 'wandering atrial pacemaker', 'wandering atrial pacemaker'], 'unspecified': ['atrial rhythm', 'atrial rhythm', 'ectopic atrial rhythm', 'ectopic atrial rhythm', 'ear', 'ear', 'ectopic atrial bradycardia', 'ectopic atrial bradycardia', 'abnormal p vector', 'abnormal p vector', 'nonsinus atrial mechanism', 'nonsinus atrial mechanism', 'p wave axis suggests atrial rather than sinus mechanism', 'p wave axis suggests atrial rather than sinus mechanism', 'low atrial bradycardia', 'low atrial bradycardia']}))

partners_ecg_read_pc_clean_ectopic_atrial_rhythm = TensorMap('partners_ecg_read_pc_clean_ectopic_atrial_rhythm', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest', channel_map={'wandering': 0, 'unifocal': 1, 'multifocal': 2, 'unspecified': 3}, tensor_from_file=make_partners_ecg_label(keys='read_pc_clean', dict_of_list = {'unifocal': ['unifocal atrial tachycardia', 'unifocal atrial tachycardia', 'unifocal ectopic atrial rhythm', 'unifocal ectopic atrial rhythm', 'unifocal ear', 'unifocal ear', 'unusual p wave axis', 'unusual p wave axis', 'low atrial pacer', 'low atrial pacer'], 'multifocal': ['multifocal atrial tachycardia', 'multifocal atrial tachycardia', 'multifocal ectopic atrial rhythm', 'multifocal ectopic atrial rhythm', 'multifocal ear', 'multifocal ear', 'dual atrial foci ', 'dual atrial foci ', 'multiple atrial foci', 'multiple atrial foci', 'multifocal atrial rhythm', 'multifocal atrial rhythm', 'multifocal atrialrhythm', 'multifocal atrialrhythm'], 'wandering': ['wandering atrial tachycardia', 'wandering atrial tachycardia', 'wandering ectopic atrial rhythm', 'wandering ectopic atrial rhythm', 'wandering ear', 'wandering ear', 'wandering atrial pacemaker', 'wandering atrial pacemaker'], 'unspecified': ['atrial rhythm', 'atrial rhythm', 'ectopic atrial rhythm', 'ectopic atrial rhythm', 'ear', 'ear', 'ectopic atrial bradycardia', 'ectopic atrial bradycardia', 'abnormal p vector', 'abnormal p vector', 'nonsinus atrial mechanism', 'nonsinus atrial mechanism', 'p wave axis suggests atrial rather than sinus mechanism', 'p wave axis suggests atrial rather than sinus mechanism', 'low atrial bradycardia', 'low atrial bradycardia']}))

partners_ecg_read_md_clean_ectopic_atrial_tachycardia = TensorMap('partners_ecg_read_md_clean_ectopic_atrial_tachycardia', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest', channel_map={'unifocal': 0, 'multifocal': 1, 'unspecified': 2}, tensor_from_file=make_partners_ecg_label(keys='read_md_clean', dict_of_list = {'unspecified': ['ectopic atrial tachycardia', 'ectopic atrial tachycardia', 'ectopic atrial tachycardia, unspecified', 'ectopic atrial tachycardia, unspecified', 'unspecified ectopic atrial tachycardia', 'unspecified ectopic atrial tachycardia'], 'unifocal': ['ectopic atrial tachycardia, unifocal', 'ectopic atrial tachycardia, unifocal', 'unifocal ectopic atrial tachycardia', 'unifocal ectopic atrial tachycardia'], 'multifocal': ['ectopic atrial tachycardia, multifocal', 'ectopic atrial tachycardia, multifocal', 'multifocal ectopic atrial tachycardia', 'multifocal ectopic atrial tachycardia']}))

partners_ecg_read_pc_clean_ectopic_atrial_tachycardia = TensorMap('partners_ecg_read_pc_clean_ectopic_atrial_tachycardia', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest', channel_map={'unifocal': 0, 'multifocal': 1, 'unspecified': 2}, tensor_from_file=make_partners_ecg_label(keys='read_pc_clean', dict_of_list = {'unspecified': ['ectopic atrial tachycardia', 'ectopic atrial tachycardia', 'ectopic atrial tachycardia, unspecified', 'ectopic atrial tachycardia, unspecified', 'unspecified ectopic atrial tachycardia', 'unspecified ectopic atrial tachycardia'], 'unifocal': ['ectopic atrial tachycardia, unifocal', 'ectopic atrial tachycardia, unifocal', 'unifocal ectopic atrial tachycardia', 'unifocal ectopic atrial tachycardia'], 'multifocal': ['ectopic atrial tachycardia, multifocal', 'ectopic atrial tachycardia, multifocal', 'multifocal ectopic atrial tachycardia', 'multifocal ectopic atrial tachycardia']}))

partners_ecg_read_md_clean_svt = TensorMap('partners_ecg_read_md_clean_svt', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest', channel_map={'avrt': 0, 'avnrt': 1, 'unspecified': 2}, tensor_from_file=make_partners_ecg_label(keys='read_md_clean', dict_of_list = {'avnrt': ['av nodal reentry tachycardia', 'av nodal reentry tachycardia', 'atrioventricular nodal reentry tachycardia', 'atrioventricular nodal reentry tachycardia', 'avnrt', 'avnrt'], 'avrt': ['atrioventricular reentrant tachycardia ', 'atrioventricular reentrant tachycardia ', 'av reentrant tachycardia ', 'av reentrant tachycardia ', 'avrt', 'avrt']}))

partners_ecg_read_pc_clean_svt = TensorMap('partners_ecg_read_pc_clean_svt', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest', channel_map={'avrt': 0, 'avnrt': 1, 'unspecified': 2}, tensor_from_file=make_partners_ecg_label(keys='read_pc_clean', dict_of_list = {'avnrt': ['av nodal reentry tachycardia', 'av nodal reentry tachycardia', 'atrioventricular nodal reentry tachycardia', 'atrioventricular nodal reentry tachycardia', 'avnrt', 'avnrt'], 'avrt': ['atrioventricular reentrant tachycardia ', 'atrioventricular reentrant tachycardia ', 'av reentrant tachycardia ', 'av reentrant tachycardia ', 'avrt', 'avrt']}))

partners_ecg_afib = TensorMap('partners_ecg_afib', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest', channel_map={'no_atrial_fibrillation': 0, 'atrial_fibrillation': 1}, tensor_from_file=make_partners_ecg_label(keys='read_md_clean', not_found_key='no_atrial_fibrillation', dict_of_list = {'atrial_fibrillation': ['atrial fib', 'afib', 'atrial fibrillation', 'atrial  fibrillation', 'afibrillation']}))
partners_ecg_afib_all = TensorMap('partners_ecg_afib_all', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest', channel_map={'no_atrial_fibrillation': 0, 'atrial_fibrillation': 1}, tensor_from_file=make_partners_ecg_label(keys=['read_md_clean', 'read_pc_clean'], not_found_key='no_atrial_fibrillation', dict_of_list = {'atrial_fibrillation': ['atrial fib', 'afib', 'atrial fibrillation', 'atrial  fibrillation', 'afibrillation']}))
partners_ecg_sinus_rhythm = TensorMap('partners_sinus_rhythm', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest', channel_map={'sinus_arrhythmia': 0, 'unspecified': 1}, tensor_from_file=make_partners_ecg_label(keys=['read_md_clean', 'read_pc_clean'], dict_of_list = {'sinus_arrhythmia': ['marked sinus arrhythmia', 'marked sinus arrhythmia', 'sinus arrhythmia', 'sinus arrhythmia']}))
partners_ecg_supranodal = TensorMap(
    'partners_ecg_supranodal', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest',
    channel_map={'retrograde_atrial_activation': 0, 'narrow_qrs_tachycardia': 1, 'unspecified': 2, 'supraventricular_tachycardia': 3, 'ectopic_atrial_tachycardia': 4, 'ectopic_atrial_rhythm': 5, 'atrial_fibrillation': 6, 'atrial_flutter': 7, 'sinus_rhythm': 8},
    tensor_from_file=make_partners_ecg_label(
        keys=['read_md_clean', 'read_pc_clean'], dict_of_list = {
            'atrial_fibrillation': ['atrial fib', 'afib', 'atrial fibrillation', 'atrial  fibrillation', 'afibrillation'],
            'sinus_rhythm': ['marked sinus arrhythmia', 'sinus arrhythmia', 'normal ecg', 'atrial bigeminal rhythm', 'atrial bigeminy and ventricular bigeminy', 'atrial trigeminy', 'atrialbigeminy', 'normal sinus rhythm', 'normal when compared with ecg of', 'rhythm has reverted to normal', 'rhythm is now clearly sinus', 'sinus bradycardia', 'sinus rhythm', 'sinus rhythm at a rate', 'sinus tachycardia', 'tracing within normal limits', 'tracing is within normal limits', 'tracing within normal limits', 'sinoatrial block', 'sa block', 'sinoatrial block, type ii', 'sa block, type i', 'type i sinoatrial block', 'type i sa block', 'sinoatrial block, type ii', 'sa block, type i', 'type ii sinoatrial block', 'type ii sa block', 'sinus pause', 'sinus arrest', 'sinus slowing', 'with occasional native sinus beats', 'frequent native sinus beats', 'conducted sinus impulses', 'sinus mechanism has replaced', 'rhythm is normal sinus', 'rhythm remains normal sinus'],
            'atrial_flutter': ['fibrillation/flutter', 'aflutter', 'atrial flutter', 'probable flutter', 'atrial flutter fixed block', 'atrial flutter variable block', 'atrial flutter unspecified block', 'tachycardia possibly flutter'],
            'ectopic_atrial_rhythm': ['unifocal atrial tachycardia', 'unifocal ectopic atrial rhythm', 'unifocal ear', 'multifocal atrial tachycardia', 'multifocal ectopic atrial rhythm', 'multifocal ear', 'wandering atrial tachycardia', 'wandering ectopic atrial rhythm', 'wandering ear', 'wandering atrial pacemaker', 'atrial rhythm', 'ectopic atrial rhythm', 'ear', 'dual atrial foci ', 'ectopic atrial bradycardia', 'multiple atrial foci', 'abnormal p vector', 'nonsinus atrial mechanism', 'p wave axis suggests atrial rather than sinus mechanism', 'low atrial bradycardia', 'multifocal atrial rhythm', 'multifocal atrialrhythm', 'unusual p wave axis', 'low atrial pacer'],
            'ectopic_atrial_tachycardia': ['ectopic atrial tachycardia', 'ectopic atrial tachycardia, unspecified', 'unspecified ectopic atrial tachycardia', 'ectopic atrial tachycardia, unifocal', 'unifocal ectopic atrial tachycardia', 'ectopic atrial tachycardia, multifocal', 'multifocal ectopic atrial tachycardia'],
            'narrow_qrs_tachycardia': ['narrow qrs tachycardia', 'narrow qrs tachycardia', 'tachycardia narrow qrs', 'narrow complex tachycardia'],
            'retrograde_atrial_activation': ['retrograde atrial activation'],
            'supraventricular_tachycardia': ['sinus arrhythmia accelerated atrioventricular junctional rhythm', 'supraventricular tachycardia', 'accelerated atrioventricular nodal rhythm', 'accelerated nodal rhythm', 'atrial tachycardia', 'av nodal reentry tachycardia', 'atrioventricular nodal reentry tachycardia', 'avnrt', 'atrioventricular reentrant tachycardia ', 'av reentrant tachycardia ', 'avrt'],
            'unspecified': ['junctional tachycardia', 'atrial arrhythmia', 'technically poor tracing ', 'accelerated idioventricular rhythm', 'atrial activity is indistinct', 'rhythm uncertain', 'rhythm unclear', 'uncertain rhythm', 'undetermined rhythm', 'supraventricular rhythm'],
        },
    ),
)

partners_ecg_supranodal_weighted = TensorMap(
    'partners_ecg_supranodal', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest',
    channel_map={'retrograde_atrial_activation': 0, 'narrow_qrs_tachycardia': 1, 'unspecified': 2, 'supraventricular_tachycardia': 3, 'ectopic_atrial_tachycardia': 4, 'ectopic_atrial_rhythm': 5, 'atrial_fibrillation': 6, 'atrial_flutter': 7, 'sinus_rhythm': 8},
    loss=weighted_crossentropy([20.0, 20.0, 10.0, 20.0, 20.0, 4.0, 1.0, 10.0, 1.0], 'partners_supranodal'),
    tensor_from_file=make_partners_ecg_label(
        keys=['read_md_clean', 'read_pc_clean'],
        dict_of_list = {
            'atrial_fibrillation': ['atrial fib', 'afib', 'atrial fibrillation', 'atrial  fibrillation', 'afibrillation'],
            'sinus_rhythm': ['marked sinus arrhythmia', 'sinus arrhythmia', 'normal ecg', 'atrial bigeminal rhythm', 'atrial bigeminy and ventricular bigeminy', 'atrial trigeminy', 'atrialbigeminy', 'normal sinus rhythm', 'normal when compared with ecg of', 'rhythm has reverted to normal', 'rhythm is now clearly sinus', 'sinus bradycardia', 'sinus rhythm', 'sinus rhythm at a rate', 'sinus tachycardia', 'tracing within normal limits', 'tracing is within normal limits', 'tracing within normal limits', 'sinoatrial block', 'sa block', 'sinoatrial block, type ii', 'sa block, type i', 'type i sinoatrial block', 'type i sa block', 'sinoatrial block, type ii', 'sa block, type i', 'type ii sinoatrial block', 'type ii sa block', 'sinus pause', 'sinus arrest', 'sinus slowing', 'with occasional native sinus beats', 'frequent native sinus beats', 'conducted sinus impulses', 'sinus mechanism has replaced', 'rhythm is normal sinus', 'rhythm remains normal sinus'],
            'atrial_flutter': ['fibrillation/flutter', 'aflutter', 'atrial flutter', 'probable flutter', 'atrial flutter fixed block', 'atrial flutter variable block', 'atrial flutter unspecified block', 'tachycardia possibly flutter'],
            'ectopic_atrial_rhythm': ['unifocal atrial tachycardia', 'unifocal ectopic atrial rhythm', 'unifocal ear', 'multifocal atrial tachycardia', 'multifocal ectopic atrial rhythm', 'multifocal ear', 'wandering atrial tachycardia', 'wandering ectopic atrial rhythm', 'wandering ear', 'wandering atrial pacemaker', 'atrial rhythm', 'ectopic atrial rhythm', 'ear', 'dual atrial foci ', 'ectopic atrial bradycardia', 'multiple atrial foci', 'abnormal p vector', 'nonsinus atrial mechanism', 'p wave axis suggests atrial rather than sinus mechanism', 'low atrial bradycardia', 'multifocal atrial rhythm', 'multifocal atrialrhythm', 'unusual p wave axis', 'low atrial pacer'],
            'ectopic_atrial_tachycardia': ['ectopic atrial tachycardia', 'ectopic atrial tachycardia, unspecified', 'unspecified ectopic atrial tachycardia', 'ectopic atrial tachycardia, unifocal', 'unifocal ectopic atrial tachycardia', 'ectopic atrial tachycardia, multifocal', 'multifocal ectopic atrial tachycardia'],
            'narrow_qrs_tachycardia': ['narrow qrs tachycardia', 'narrow qrs tachycardia', 'tachycardia narrow qrs', 'narrow complex tachycardia'],
            'retrograde_atrial_activation': ['retrograde atrial activation'],
            'supraventricular_tachycardia': ['sinus arrhythmia accelerated atrioventricular junctional rhythm', 'supraventricular tachycardia', 'accelerated atrioventricular nodal rhythm', 'accelerated nodal rhythm', 'atrial tachycardia', 'av nodal reentry tachycardia', 'atrioventricular nodal reentry tachycardia', 'avnrt', 'atrioventricular reentrant tachycardia ', 'av reentrant tachycardia ', 'avrt'],
            'unspecified': ['junctional tachycardia', 'atrial arrhythmia', 'technically poor tracing ', 'accelerated idioventricular rhythm', 'atrial activity is indistinct', 'rhythm uncertain', 'rhythm unclear', 'uncertain rhythm', 'undetermined rhythm', 'supraventricular rhythm'],
        },
    ),
)

partners_ecg_supranodal_6 = TensorMap(
    'partners_ecg_supranodal_6', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest',
    channel_map={'supraventricular_tachycardia': 0, 'ectopic_atrial_rhythm': 1, 'atrial_flutter': 2, 'atrial_fibrillation': 3, 'sinus_rhythm': 4, 'unspecified': 5},
    tensor_from_file=make_partners_ecg_label(
        keys=['read_md_clean', 'read_pc_clean'],
        dict_of_list = {
            'atrial_fibrillation': ['atrial fib', 'afib', 'atrial fibrillation', 'atrial  fibrillation', 'afibrillation'],
            'sinus_rhythm': ['marked sinus arrhythmia', 'sinus arrhythmia', 'normal ecg', 'atrial bigeminal rhythm', 'atrial bigeminy and ventricular bigeminy', 'atrial trigeminy', 'atrialbigeminy', 'normal sinus rhythm', 'normal when compared with ecg of', 'rhythm has reverted to normal', 'rhythm is now clearly sinus', 'sinus bradycardia', 'sinus rhythm', 'sinus rhythm at a rate', 'sinus tachycardia', 'tracing within normal limits', 'tracing is within normal limits', 'tracing within normal limits', 'sinoatrial block', 'sa block', 'sinoatrial block, type ii', 'sa block, type i', 'type i sinoatrial block', 'type i sa block', 'sinoatrial block, type ii', 'sa block, type i', 'type ii sinoatrial block', 'type ii sa block', 'sinus pause', 'sinus arrest', 'sinus slowing', 'with occasional native sinus beats', 'frequent native sinus beats', 'conducted sinus impulses', 'sinus mechanism has replaced', 'rhythm is normal sinus', 'rhythm remains normal sinus'],
            'atrial_flutter': ['fibrillation/flutter', 'aflutter', 'atrial flutter', 'probable flutter', 'atrial flutter fixed block', 'atrial flutter variable block', 'atrial flutter unspecified block', 'tachycardia possibly flutter'],
            'ectopic_atrial_rhythm': ['unifocal atrial tachycardia', 'unifocal ectopic atrial rhythm', 'unifocal ear', 'multifocal atrial tachycardia', 'multifocal ectopic atrial rhythm', 'multifocal ear', 'wandering atrial tachycardia', 'wandering ectopic atrial rhythm', 'wandering ear', 'wandering atrial pacemaker', 'atrial rhythm', 'ectopic atrial rhythm', 'ear', 'dual atrial foci ', 'ectopic atrial bradycardia', 'multiple atrial foci', 'abnormal p vector', 'nonsinus atrial mechanism', 'p wave axis suggests atrial rather than sinus mechanism', 'low atrial bradycardia', 'multifocal atrial rhythm', 'multifocal atrialrhythm', 'unusual p wave axis', 'low atrial pacer'],
            'supraventricular_tachycardia': ['sinus arrhythmia accelerated atrioventricular junctional rhythm', 'supraventricular tachycardia', 'accelerated atrioventricular nodal rhythm', 'accelerated nodal rhythm', 'atrial tachycardia', 'av nodal reentry tachycardia', 'atrioventricular nodal reentry tachycardia', 'avnrt', 'atrioventricular reentrant tachycardia ', 'av reentrant tachycardia ', 'avrt'], 'unspecified': ['junctional tachycardia', 'atrial arrhythmia', 'technically poor tracing ', 'accelerated idioventricular rhythm', 'atrial activity is indistinct', 'rhythm uncertain', 'rhythm unclear', 'uncertain rhythm', 'undetermined rhythm', 'supraventricular rhythm'],
        },
    ),
)

partners_ecg_supranodal_6_v6_valid = TensorMap(
    'partners_ecg_supranodal_6', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest',
    validator=v6_zeros_validator,
    channel_map={'supraventricular_tachycardia': 0, 'ectopic_atrial_rhythm': 1, 'atrial_flutter': 2, 'atrial_fibrillation': 3, 'sinus_rhythm': 4, 'unspecified': 5},
    tensor_from_file=make_partners_ecg_label(
        keys=['read_md_clean', 'read_pc_clean', 'testreason'],
        dict_of_list = {
            'atrial_fibrillation': ['atrial fib', 'afib', 'atrial fibrillation', 'atrial  fibrillation', 'afibrillation'],
            'sinus_rhythm': ['marked sinus arrhythmia', 'sinus arrhythmia', 'normal ecg', 'atrial bigeminal rhythm', 'atrial bigeminy and ventricular bigeminy', 'atrial trigeminy', 'atrialbigeminy', 'normal sinus rhythm', 'normal when compared with ecg of', 'rhythm has reverted to normal', 'rhythm is now clearly sinus', 'sinus bradycardia', 'sinus rhythm', 'sinus rhythm at a rate', 'sinus tachycardia', 'tracing within normal limits', 'tracing is within normal limits', 'tracing within normal limits', 'sinoatrial block', 'sa block', 'sinoatrial block, type ii', 'sa block, type i', 'type i sinoatrial block', 'type i sa block', 'sinoatrial block, type ii', 'sa block, type i', 'type ii sinoatrial block', 'type ii sa block', 'sinus pause', 'sinus arrest', 'sinus slowing', 'with occasional native sinus beats', 'frequent native sinus beats', 'conducted sinus impulses', 'sinus mechanism has replaced', 'rhythm is normal sinus', 'rhythm remains normal sinus'],
            'atrial_flutter': ['fibrillation/flutter', 'aflutter', 'atrial flutter', 'probable flutter', 'atrial flutter fixed block', 'atrial flutter variable block', 'atrial flutter unspecified block', 'tachycardia possibly flutter'],
            'ectopic_atrial_rhythm': ['unifocal atrial tachycardia', 'unifocal ectopic atrial rhythm', 'unifocal ear', 'multifocal atrial tachycardia', 'multifocal ectopic atrial rhythm', 'multifocal ear', 'wandering atrial tachycardia', 'wandering ectopic atrial rhythm', 'wandering ear', 'wandering atrial pacemaker', 'atrial rhythm', 'ectopic atrial rhythm', 'ear', 'dual atrial foci ', 'ectopic atrial bradycardia', 'multiple atrial foci', 'abnormal p vector', 'nonsinus atrial mechanism', 'p wave axis suggests atrial rather than sinus mechanism', 'low atrial bradycardia', 'multifocal atrial rhythm', 'multifocal atrialrhythm', 'unusual p wave axis', 'low atrial pacer'],
            'supraventricular_tachycardia': ['sinus arrhythmia accelerated atrioventricular junctional rhythm', 'supraventricular tachycardia', 'accelerated atrioventricular nodal rhythm', 'accelerated nodal rhythm', 'atrial tachycardia', 'av nodal reentry tachycardia', 'atrioventricular nodal reentry tachycardia', 'avnrt', 'atrioventricular reentrant tachycardia ', 'av reentrant tachycardia ', 'avrt'], 'unspecified': ['junctional tachycardia', 'atrial arrhythmia', 'technically poor tracing ', 'accelerated idioventricular rhythm', 'atrial activity is indistinct', 'rhythm uncertain', 'rhythm unclear', 'uncertain rhythm', 'undetermined rhythm', 'supraventricular rhythm'],
        },
    ),
)

partners_ecg_supranodal_6_weighted = TensorMap(
    'partners_ecg_supranodal_6', interpretation=Interpretation.CATEGORICAL, time_series_limit=0, path_prefix='partners_ecg_rest',
    channel_map={'supraventricular_tachycardia': 0, 'ectopic_atrial_rhythm': 1, 'atrial_flutter': 2, 'atrial_fibrillation': 3, 'sinus_rhythm': 4, 'unspecified': 5},
    loss=weighted_crossentropy([70.0, 10.0, 50.0, 7.0, 1.0, 10.0], 'partners_supranodal_6'),
    tensor_from_file=make_partners_ecg_label(
        keys=['read_md_clean', 'read_pc_clean'],
        dict_of_list = {
            'atrial_fibrillation': ['atrial fib', 'afib', 'atrial fibrillation', 'atrial  fibrillation', 'afibrillation'],
            'sinus_rhythm': ['marked sinus arrhythmia', 'sinus arrhythmia', 'normal ecg', 'atrial bigeminal rhythm', 'atrial bigeminy and ventricular bigeminy', 'atrial trigeminy', 'atrialbigeminy', 'normal sinus rhythm', 'normal when compared with ecg of', 'rhythm has reverted to normal', 'rhythm is now clearly sinus', 'sinus bradycardia', 'sinus rhythm', 'sinus rhythm at a rate', 'sinus tachycardia', 'tracing within normal limits', 'tracing is within normal limits', 'tracing within normal limits', 'sinoatrial block', 'sa block', 'sinoatrial block, type ii', 'sa block, type i', 'type i sinoatrial block', 'type i sa block', 'sinoatrial block, type ii', 'sa block, type i', 'type ii sinoatrial block', 'type ii sa block', 'sinus pause', 'sinus arrest', 'sinus slowing', 'with occasional native sinus beats', 'frequent native sinus beats', 'conducted sinus impulses', 'sinus mechanism has replaced', 'rhythm is normal sinus', 'rhythm remains normal sinus'],
            'atrial_flutter': ['fibrillation/flutter', 'aflutter', 'atrial flutter', 'probable flutter', 'atrial flutter fixed block', 'atrial flutter variable block', 'atrial flutter unspecified block', 'tachycardia possibly flutter'],
            'ectopic_atrial_rhythm': ['unifocal atrial tachycardia', 'unifocal ectopic atrial rhythm', 'unifocal ear', 'multifocal atrial tachycardia', 'multifocal ectopic atrial rhythm', 'multifocal ear', 'wandering atrial tachycardia', 'wandering ectopic atrial rhythm', 'wandering ear', 'wandering atrial pacemaker', 'atrial rhythm', 'ectopic atrial rhythm', 'ear', 'dual atrial foci ', 'ectopic atrial bradycardia', 'multiple atrial foci', 'abnormal p vector', 'nonsinus atrial mechanism', 'p wave axis suggests atrial rather than sinus mechanism', 'low atrial bradycardia', 'multifocal atrial rhythm', 'multifocal atrialrhythm', 'unusual p wave axis', 'low atrial pacer'],
            'supraventricular_tachycardia': ['sinus arrhythmia accelerated atrioventricular junctional rhythm', 'supraventricular tachycardia', 'accelerated atrioventricular nodal rhythm', 'accelerated nodal rhythm', 'atrial tachycardia', 'av nodal reentry tachycardia', 'atrioventricular nodal reentry tachycardia', 'avnrt', 'atrioventricular reentrant tachycardia ', 'av reentrant tachycardia ', 'avrt'], 'unspecified': ['junctional tachycardia', 'atrial arrhythmia', 'technically poor tracing ', 'accelerated idioventricular rhythm', 'atrial activity is indistinct', 'rhythm uncertain', 'rhythm unclear', 'uncertain rhythm', 'undetermined rhythm', 'supraventricular rhythm'],
        },
    ),
)
