import os
import csv
import h5py
import logging
import datetime
import numpy as np
from collections import defaultdict
from typing import Dict, List, Callable, Union, Tuple

from ml4cvd.tensor_maps_by_hand import TMAPS
from ml4cvd.defines import ECG_REST_AMP_LEADS, PARTNERS_DATE_FORMAT
from ml4cvd.TensorMap import TensorMap, str2date, Interpretation, make_range_validator, decompress_data, TimeSeriesOrder

YEAR_DAYS = 365.26
INCIDENCE_CSV = '/media/erisone_snf13/lc_outcomes.csv'
PARTNERS_PREFIX = 'partners_ecg_rest'


def _get_ecg_dates(tm, hd5):
    dates = list(hd5[tm.path_prefix])
    if tm.time_series_order == TimeSeriesOrder.NEWEST:
        dates.sort()
    elif tm.time_series_order == TimeSeriesOrder.OLDEST:
        dates.sort(reverse=True)
    elif tm.time_series_order == TimeSeriesOrder.RANDOM:
        np.random.shuffle(dates)
    else:
        raise ValueError(f'Unknown option "{tm.time_series_order}" passed for which tensors to use in multi tensor HD5')
    start_idx = tm.time_series_limit if tm.time_series_limit is not None else 1
    dates = dates[-start_idx:] # If num_tensors is 0, get all tensors
    dates.sort(reverse=True)
    return dates


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


def make_voltage(population_normalize: float = None):
    def get_voltage_from_file(tm, hd5, dependents={}):
        ecg_dates = _get_ecg_dates(tm, hd5)
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        tensor = np.zeros(shape, dtype=np.float32)
        for i, ecg_date in enumerate(ecg_dates):
            for cm in tm.channel_map:
                try:
                    path = _make_hd5_path(tm, ecg_date, cm)
                    voltage = decompress_data(data_compressed=hd5[path][()], dtype=hd5[path].attrs['dtype'])
                    voltage = _resample_voltage(voltage, shape[1] if dynamic else shape[0])
                    slices = (i, ..., tm.channel_map[cm]) if dynamic else (..., tm.channel_map[cm])
                    tensor[slices] = voltage
                except KeyError:
                    pass
        if population_normalize is None:
            tm.normalization = {'zero_mean_std1': True}
        else:
            tensor /= population_normalize
        return tensor
    return get_voltage_from_file


TMAPS['partners_ecg_voltage'] = TensorMap(
    'partners_ecg_voltage',
    shape=(None, 2500, 12),
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_voltage(population_normalize=2000.0),
    channel_map=ECG_REST_AMP_LEADS,
    time_series_limit=0,
)


TMAPS['partners_ecg_voltage_newest'] = TensorMap(
    'partners_ecg_voltage_newest',
    shape=(2500, 12),
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_voltage(population_normalize=2000.0),
    channel_map=ECG_REST_AMP_LEADS,
)


TMAPS['partners_ecg_2500'] = TensorMap('ecg_rest_2500', shape=(None, 2500, 12), path_prefix=PARTNERS_PREFIX, tensor_from_file=make_voltage(), channel_map=ECG_REST_AMP_LEADS, time_series_limit=0)
TMAPS['partners_ecg_5000'] = TensorMap('ecg_rest_5000', shape=(None, 5000, 12), path_prefix=PARTNERS_PREFIX, tensor_from_file=make_voltage(), channel_map=ECG_REST_AMP_LEADS, time_series_limit=0)
TMAPS['partners_ecg_2500_raw'] = TensorMap('ecg_rest_2500_raw', shape=(None, 2500, 12), path_prefix=PARTNERS_PREFIX, tensor_from_file=make_voltage(population_normalize=2000.0), channel_map=ECG_REST_AMP_LEADS, time_series_limit=0)
TMAPS['partners_ecg_5000_raw'] = TensorMap('ecg_rest_5000_raw', shape=(None, 5000, 12), path_prefix=PARTNERS_PREFIX, tensor_from_file=make_voltage(population_normalize=2000.0), channel_map=ECG_REST_AMP_LEADS, time_series_limit=0)
TMAPS['partners_ecg_2500_newest'] = TensorMap('ecg_rest_2500_newest', shape=(2500, 12), path_prefix=PARTNERS_PREFIX, tensor_from_file=make_voltage(), channel_map=ECG_REST_AMP_LEADS)
TMAPS['partners_ecg_5000_newest'] = TensorMap('ecg_rest_5000_newest', shape=(5000, 12), path_prefix=PARTNERS_PREFIX, tensor_from_file=make_voltage(), channel_map=ECG_REST_AMP_LEADS)
TMAPS['partners_ecg_2500_raw_newest'] = TensorMap('ecg_rest_2500_raw_newest', shape=(2500, 12), path_prefix=PARTNERS_PREFIX, tensor_from_file=make_voltage(population_normalize=2000.0), channel_map=ECG_REST_AMP_LEADS)
TMAPS['partners_ecg_5000_raw_newest'] = TensorMap('ecg_rest_5000_raw_newest', shape=(5000, 12), path_prefix=PARTNERS_PREFIX, tensor_from_file=make_voltage(population_normalize=2000.0), channel_map=ECG_REST_AMP_LEADS)


def make_voltage_attr(volt_attr: str = ""):
    def get_voltage_attr_from_file(tm, hd5, dependents={}):
        ecg_dates = _get_ecg_dates(tm, hd5)
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        tensor = np.zeros(shape, dtype=np.float32)
        for i, ecg_date in enumerate(ecg_dates):
            for cm in tm.channel_map:
                try:
                    path = _make_hd5_path(tm, ecg_date, cm)
                    slices = (i, tm.channel_map[cm]) if dynamic else (tm.channel_map,)
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


TMAPS["voltage_len_newest"] = TensorMap(
    "voltage_len_newest",
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_voltage_attr(volt_attr="len"),
    shape=(12,),
    channel_map=ECG_REST_AMP_LEADS,
)


TMAPS["len_i"] = TensorMap("len_i", shape=(None, 1), path_prefix=PARTNERS_PREFIX, tensor_from_file=make_voltage_attr(volt_attr="len"), channel_map={'I': 0}, time_series_limit=0)
TMAPS["len_v6"] = TensorMap("len_v6", shape=(None, 1), path_prefix=PARTNERS_PREFIX, tensor_from_file=make_voltage_attr(volt_attr="len"), channel_map={'V6': 0}, time_series_limit=0)
TMAPS["len_i_newest"] = TensorMap("len_i_newest", shape=(1,), path_prefix=PARTNERS_PREFIX, tensor_from_file=make_voltage_attr(volt_attr="len"), channel_map={'I': 0})
TMAPS["len_v6_newest"] = TensorMap("len_v6_newest", shape=(1,), path_prefix=PARTNERS_PREFIX, tensor_from_file=make_voltage_attr(volt_attr="len"), channel_map={'V6': 0})


def make_partners_ecg_label(keys: Union[str, List[str]] = "read_md_clean", dict_of_list: Dict = dict(), not_found_key: str = "unspecified"):
    if type(keys) == str:
        keys = [keys]
    def get_partners_ecg_label(tm, hd5, dependents={}):
        ecg_dates = _get_ecg_dates(tm, hd5)
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        label_array = np.zeros(shape, dtype=np.float32)
        for i, ecg_date in enumerate(ecg_dates):
            found = False
            for key in keys:
                path = _make_hd5_path(tm, ecg_date, key)
                if path not in hd5:
                    continue
                try:
                    read = decompress_data(data_compressed=hd5[path][()], dtype=hd5[path].attrs['dtype'])
                    for channel, idx in sorted(tm.channel_map.items(), key=lambda cm: cm[1]):
                        if channel in dict_of_list:
                            for string in dict_of_list[channel]:
                                if string in read:
                                    slices = (i, idx) if dynamic else (idx,)
                                    label_array[slices] = 1
                                    found = True
                                if found: break
                        if found: break
                    if found: break
                except KeyError:
                    pass
            if not found:
                slices = (i, tm.channel_map[not_found_key]) if dynamic else (tm.channel_map[not_found_key],)
                label_array[slices] = 1
        return label_array
    return get_partners_ecg_label


def validator_not_empty(tm: TensorMap, tensor: np.ndarray, hd5: h5py.File):
    if any(tensor == ''):
        raise ValueError(f'TensorMap {tm.name} failed empty string check.')


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
    validator=validator_not_empty,
)


task = "partners_ecg_datetime_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=partners_ecg_datetime,
    shape=(1,),
    validator=validator_not_empty,
)


def make_partners_ecg_tensor(key: str, fill: float = 0):
    def get_partners_ecg_tensor(tm, hd5, dependents={}):
        ecg_dates = _get_ecg_dates(tm, hd5)
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        if tm.interpretation == Interpretation.LANGUAGE:
            tensor = np.full(shape, '', dtype=object)
        elif tm.interpretation == Interpretation.CONTINUOUS:
            tensor = np.zeros(shape, dtype=np.float32) if fill == 0 else np.full(shape, fill, dtype=np.float32)
        elif tm.interpretation == Interpretation.CATEGORICAL:
            tensor = np.zeros(shape, dtype=float) if fill == 0 else np.full(shape, fill, dtype=float)

        for i, ecg_date in enumerate(ecg_dates):
            path = _make_hd5_path(tm, ecg_date, key)
            try:
                tensor[i] = decompress_data(data_compressed=hd5[path][()], dtype='str')
            except KeyError:
                logging.debug(f'Could not obtain tensor {tm.name} from ECG on {ecg_date} in {hd5.filename}')

        if tm.interpretation == Interpretation.LANGUAGE:
            tensor = tensor.astype(str)
        return tensor
    return get_partners_ecg_tensor


task = "partners_ecg_read_md_raw"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="read_md_clean"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_not_empty,
)


task = "partners_ecg_read_md_raw"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="read_md_clean"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_not_empty,
)


task = "partners_ecg_read_md_raw_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="read_md_clean"),
    shape=(1,),
    validator=validator_not_empty,
)


task = "partners_ecg_read_pc_raw"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="read_pc_clean"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_not_empty,
)


task = "partners_ecg_read_pc_raw_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="read_pc_clean"),
    shape=(1,),
    validator=validator_not_empty,
)


# TODO do we still need the cross reference tmaps?
def validator_cross_reference(tm: TensorMap, tensor: np.ndarray):
    if int(tensor) not in tm.cross_reference:
        raise ValueError(f"Skipping TensorMap {tm.name} not found in Apollo.")


def create_cross_reference_dict(fpath="/data/apollo/demographics.csv"):
    try:
        with open(fpath, mode="r") as f:
            reader = csv.reader(f)
            next(reader)
            cross_reference_dict = {int(rows[0]):None for rows in reader}
        return cross_reference_dict
    except FileNotFoundError:
        return {}


task = "partners_ecg_patientid_cross_reference_apollo"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="patientid"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_cross_reference,
)


TMAPS[task].cross_reference = create_cross_reference_dict()


task = "partners_ecg_patientid"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="patientid"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_not_empty,
)


task = "partners_ecg_patientid_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="patientid"),
    shape=(1,),
    validator=validator_not_empty,
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


task = "partners_ecg_patientid_clean_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="patientid_clean"),
    shape=(1,),
    validator=validator_clean_mrn,
)


task = "partners_ecg_firstname"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="patientfirstname"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_not_empty,
)


task = "partners_ecg_firstname_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="patientfirstname"),
    shape=(1,),
    validator=validator_not_empty,
)


task = "partners_ecg_lastname"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="patientlastname"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_not_empty,
)


task = "partners_ecg_lastname_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="patientlastname"),
    shape=(1,),
    validator=validator_not_empty,
)


task = "partners_ecg_gender"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="gender"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_not_empty,
)

task = "partners_ecg_date"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="acquisitiondate"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_not_empty,
)


task = "partners_ecg_date_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="acquisitiondate"),
    shape=(1,),
    validator=validator_not_empty,
)


task = "partners_ecg_time"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="acquisitiontime"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_not_empty,
)


task = "partners_ecg_time_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="acquisitiontime"),
    shape=(1,),
    validator=validator_not_empty,
)


task = "partners_ecg_sitename"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="sitename"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_not_empty,
)


task = "partners_ecg_sitename_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="sitename"),
    shape=(1,),
    validator=validator_not_empty,
)


task = "partners_ecg_location"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="location"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_not_empty,
)


task = "partners_ecg_location_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="location"),
    shape=(1,),
    validator=validator_not_empty,
)


task = "partners_ecg_dob"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="dateofbirth"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_not_empty,
)


task = "partners_ecg_dob_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.LANGUAGE,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="dateofbirth"),
    shape=(1,),
    validator=validator_not_empty,
)


task = "partners_ecg_sampling_frequency"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="ecgsamplebase_pc"),
    shape=(None, 1),
    time_series_limit=0,
    validator=validator_not_empty,
)


task = "partners_ecg_sampling_frequency_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    tensor_from_file=make_partners_ecg_tensor(key="ecgsamplebase_pc"),
    shape=(1,),
    validator=validator_not_empty,
)


task = "partners_ecg_rate"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="ventricularrate_pc"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(10, 200),
)


task = "partners_ecg_rate_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="ventricularrate_pc"),
    shape=(1,),
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


task = "partners_ecg_rate_md_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="ventricularrate_md"),
    shape=(1,),
    validator=make_range_validator(10, 200),
)


TMAPS['partners_ventricular_rate'] = TensorMap(
    'VentricularRate',
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="ventricularrate_pc"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(10, 200),
    normalization={'mean': 59.3, 'std': 10.6},
)


TMAPS['partners_ventricular_rate_newest'] = TensorMap(
    'VentricularRate_newest',
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="ventricularrate_pc"),
    shape=(1,),
    validator=make_range_validator(10, 200),
    normalization={'mean': 59.3, 'std': 10.6},
)


task = "partners_ecg_qrs"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="qrsduration_pc"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(20, 400),
)


task = "partners_ecg_qrs_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="qrsduration_pc"),
    shape=(1,),
    validator=make_range_validator(20, 400),
)


task = "partners_ecg_qrs_md"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="qrsduration_md"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(20, 400),
)


task = "partners_ecg_qrs_md_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="qrsduration_md"),
    shape=(1,),
    validator=make_range_validator(20, 400),
)


task = "partners_ecg_pr"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="printerval_pc"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(50, 500),
)


task = "partners_ecg_pr_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="printerval_pc"),
    shape=(1,),
    validator=make_range_validator(50, 500),
)


task = "partners_ecg_pr_md"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="printerval_md"),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(50, 500),
)


task = "partners_ecg_pr_md_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="printerval_md"),
    shape=(1,),
    validator=make_range_validator(50, 500),
)


task = "partners_ecg_qt"
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


task = "partners_ecg_qt_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="qtinterval_pc"),
    shape=(1,),
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


task = "partners_ecg_qt_md_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="qtinterval_md"),
    shape=(1,),
    validator=make_range_validator(100, 800),
)


task = "partners_ecg_qtc"
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


task = "partners_ecg_qtc_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="qtcorrected_pc"),
    shape=(1,),
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


task = "partners_ecg_qtc_md_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="qtcorrected_md"),
    shape=(1,),
    validator=make_range_validator(100, 800),
)


task = "partners_ecg_paxis"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="paxis_pc", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_paxis_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="paxis_pc", fill=999),
    shape=(1,),
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_paxis_md"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="paxis_md", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_paxis_md_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="paxis_md", fill=999),
    shape=(1,),
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_raxis"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="raxis_pc", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_raxis_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="raxis_pc", fill=999),
    shape=(1,),
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_raxis_md"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="raxis_md", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_raxis_md_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="raxis_md", fill=999),
    shape=(1,),
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_taxis"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="taxis_pc", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_taxis_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="taxis_pc", fill=999),
    shape=(1,),
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_taxis_md"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="taxis_md", fill=999),
    shape=(None, 1),
    time_series_limit=0,
    validator=make_range_validator(-180, 180),
)


task = "partners_ecg_taxis_md_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    metrics=['mse'],
    tensor_from_file=make_partners_ecg_tensor(key="taxis_md", fill=999),
    shape=(1,),
    validator=make_range_validator(-180, 180),
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


task = "partners_ecg_weight_lbs_newest"
TMAPS[task] = TensorMap(
    task,
    interpretation=Interpretation.CONTINUOUS,
    path_prefix=PARTNERS_PREFIX,
    loss='logcosh',
    tensor_from_file=make_partners_ecg_tensor(key="weightlbs"),
    shape=(1,),
    validator=make_range_validator(100, 800),
)


def _partners_str2date(d):
    return datetime.datetime.strptime(d, PARTNERS_DATE_FORMAT).date()


def partners_ecg_age(tm, hd5, dependents={}):
    ecg_dates = _get_ecg_dates(tm, hd5)
    dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
    tensor = np.zeros(shape, dtype=float)
    for i, ecg_date in enumerate(ecg_dates):
        path = lambda key: _make_hd5_path(tm, ecg_date, key)
        try:
            birthday = decompress_data(data_compressed=hd5[path('dateofbirth')][()], dtype='str')
            acquisition = decompress_data(data_compressed=hd5[path('acquisitiondate')][()], dtype='str')
            delta = _partners_str2date(acquisition) - _partners_str2date(birthday)
            years = delta.days / YEAR_DAYS
            tensor[i] = years
        except KeyError:
            logging.debug(f'Could not calculate patient age from birthday, defer to age in ECG on {ecg_date} in {hd5.filename}')
            try:
                tensor[i] = decompress_data(data_compressed=hd5[path('patientage')][()], dtype='str')
            except KeyError:
                logging.debug(f'Could not get patient age from ECG on {ecg_date} in {hd5.filename}')
    return tensor


TMAPS['partners_ecg_age'] = TensorMap('partners_ecg_age', path_prefix=PARTNERS_PREFIX, loss='logcosh', tensor_from_file=partners_ecg_age, shape=(None, 1), time_series_limit=0)
TMAPS['partners_ecg_age_newest'] = TensorMap('partners_ecg_age_newest', path_prefix=PARTNERS_PREFIX, loss='logcosh', tensor_from_file=partners_ecg_age, shape=(1,))


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
TMAPS['partners_ecg_acquisition_year_newest'] = TensorMap('partners_ecg_acquisition_year_newest', path_prefix=PARTNERS_PREFIX, loss='logcosh',  tensor_from_file=partners_ecg_acquisition_year, shape=(1,))


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
TMAPS['partners_ecg_bmi_newest'] = TensorMap('partners_ecg_bmi_newest', path_prefix=PARTNERS_PREFIX, channel_map={'bmi': 0}, tensor_from_file=partners_bmi)


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


TMAPS['partners_ecg_race_newest'] = TensorMap(
    'partners_ecg_race_newest', interpretation=Interpretation.CATEGORICAL, path_prefix=PARTNERS_PREFIX, channel_map={'asian': 0, 'black': 1, 'hispanic': 2, 'white': 3, 'unknown': 4},
    tensor_from_file=partners_channel_string('race', race_synonyms),
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


TMAPS['partners_adult_gender'] = TensorMap(
    'adult_gender', interpretation=Interpretation.CATEGORICAL, path_prefix=PARTNERS_PREFIX, channel_map={'female': 0, 'male': 1},
    tensor_from_file=_partners_adult('gender'), time_series_limit=0,
)


TMAPS['partners_adult_gender_newest'] = TensorMap(
    'adult_gender_newest', interpretation=Interpretation.CATEGORICAL, path_prefix=PARTNERS_PREFIX, channel_map={'female': 0, 'male': 1},
    tensor_from_file=_partners_adult('gender'),
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


TMAPS["lead_i_zeros"] = TensorMap("lead_i_zeros", shape=(None, 1), path_prefix=PARTNERS_PREFIX, tensor_from_file=voltage_zeros, channel_map={'I': 0}, time_series_limit=0)
TMAPS["lead_v6_zeros"] = TensorMap("lead_v6_zeros", shape=(None, 1), path_prefix=PARTNERS_PREFIX, tensor_from_file=voltage_zeros, channel_map={'V6': 0}, time_series_limit=0)
TMAPS["lead_i_zeros_newest"] = TensorMap("lead_i_zeros_newest", shape=(1,), path_prefix=PARTNERS_PREFIX, tensor_from_file=voltage_zeros, channel_map={'I': 0})
TMAPS["lead_v6_zeros_newest"] = TensorMap("lead_v6_zeros_newest", shape=(1,), path_prefix=PARTNERS_PREFIX, tensor_from_file=voltage_zeros, channel_map={'V6': 0})


def v6_zeros_validator(tm: TensorMap, tensor: np.ndarray, hd5: h5py.File):
    voltage = decompress_data(data_compressed=hd5['V6'][()], dtype=hd5['V6'].attrs['dtype'])
    if np.count_nonzero(voltage == 0) > 10:
        raise ValueError(f'TensorMap {tm.name} has too many zeros in V6.')


def _loyalty_str2date(date_string: str):
    return str2date(date_string.split(' ')[0])


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
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        categorical_data = np.zeros(shape, dtype=np.float32)
        for i, ecg_date in enumerate(ecg_dates):
            path = lambda key: _make_hd5_path(tm, ecg_date, key)
            file_split = os.path.basename(hd5.filename).split('-')
            mrn = file_split[0]
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
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        tensor = np.zeros(shape, dtype=np.float32)
        for i, ecg_date in enumerate(ecg_dates):
            file_split = os.path.basename(hd5.filename).split('-')
            patient_key_from_ecg = int(file_split[0])

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

            if incidence_only and censor_date <= assess_date:
                raise ValueError(f'{tm.name} only considers incident diagnoses')

            tensor[(i, 0) if dynamic else (0,)] = has_disease
            tensor[(i, 1) if dynamic else (1,)] = (censor_date - assess_date).days
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
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        survival_then_censor = np.zeros(shape, dtype=np.float32)
        for ed, ecg_date in enumerate(ecg_dates):
            file_split = os.path.basename(hd5.filename).split('-')
            patient_key_from_ecg = int(file_split[0])

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
                survival_then_censor[(ed, i) if dynamic else (i,)] = float(cur_date < censor_date)
                survival_then_censor[(ed, intervals+i) if dynamic else (intervals+i,)] = has_disease * float(censor_date <= cur_date < censor_date + datetime.timedelta(days=days_per_interval))
                if i == 0 and censor_date <= cur_date:  # Handle prevalent diseases
                    survival_then_censor[(ed, intervals) if dynamic else (intervals,)] = has_disease
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

    for diagnosis in diagnosis2column:
        # Build diagnosis classification TensorMaps
        name = f'diagnosis_{diagnosis}'
        if name in needed_tensor_maps:
            tensor_from_file_fxn = build_incidence_tensor_from_file(INCIDENCE_CSV, diagnosis_column=diagnosis2column[diagnosis])
            name2tensormap[name] = TensorMap(name, Interpretation.CATEGORICAL, path_prefix=PARTNERS_PREFIX, channel_map=_diagnosis_channels(diagnosis), tensor_from_file=tensor_from_file_fxn, time_series_limit=0)
            name2tensormap[f'{name}_newest'] = TensorMap(f'{name}_newest', Interpretation.CATEGORICAL, path_prefix=PARTNERS_PREFIX, channel_map=_diagnosis_channels(diagnosis), tensor_from_file=tensor_from_file_fxn)
        name = f'incident_diagnosis_{diagnosis}'
        if name in needed_tensor_maps:
            tensor_from_file_fxn = build_incidence_tensor_from_file(INCIDENCE_CSV, diagnosis_column=diagnosis2column[diagnosis], incidence_only=True)
            name2tensormap[name] = TensorMap(name, Interpretation.CATEGORICAL, path_prefix=PARTNERS_PREFIX, channel_map=_diagnosis_channels(diagnosis, incidence_only=True), tensor_from_file=tensor_from_file_fxn, time_series_limit=0)
            name2tensormap[f'{name}_newest'] = TensorMap(f'{name}_newest', Interpretation.CATEGORICAL, path_prefix=PARTNERS_PREFIX, channel_map=_diagnosis_channels(diagnosis, incidence_only=True), tensor_from_file=tensor_from_file_fxn)

        # Build time to event TensorMaps
        name = f'cox_{diagnosis}'
        if name in needed_tensor_maps:
            tff = loyalty_time_to_event(INCIDENCE_CSV, diagnosis_column=diagnosis2column[diagnosis])
            name2tensormap[name] = TensorMap(name, Interpretation.TIME_TO_EVENT, path_prefix=PARTNERS_PREFIX, tensor_from_file=tff, time_series_limit=0)
            name2tensormap[f'{name}_newest'] = TensorMap(f'{name}_newest', Interpretation.TIME_TO_EVENT, path_prefix=PARTNERS_PREFIX, tensor_from_file=tff)
        name = f'incident_cox_{diagnosis}'
        if name in needed_tensor_maps:
            tff = loyalty_time_to_event(INCIDENCE_CSV, diagnosis_column=diagnosis2column[diagnosis], incidence_only=True)
            name2tensormap[name] = TensorMap(name, Interpretation.TIME_TO_EVENT, path_prefix=PARTNERS_PREFIX, tensor_from_file=tff, time_series_limit=0)
            name2tensormap[f'{name}_newest'] = TensorMap(f'{name}_newest', Interpretation.TIME_TO_EVENT, path_prefix=PARTNERS_PREFIX, tensor_from_file=tff)

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
                name2tensormap[needed_name] = TensorMap(needed_name, Interpretation.SURVIVAL_CURVE, path_prefix=PARTNERS_PREFIX, shape=(None, 50), days_window=days_window, tensor_from_file=tff, time_series_limit=0)
                name2tensormap[f'{needed_name}_newest'] = TensorMap(f'{needed_name}_newest', Interpretation.SURVIVAL_CURVE, path_prefix=PARTNERS_PREFIX, shape=(50,), days_window=days_window, tensor_from_file=tff)
            name = f'incident_survival_{diagnosis}'
            if name in needed_name:
                tff = _survival_from_file(days_window, INCIDENCE_CSV, diagnosis_column=diagnosis2column[diagnosis], incidence_only=True)
                name2tensormap[needed_name] = TensorMap(needed_name, Interpretation.SURVIVAL_CURVE, path_prefix=PARTNERS_PREFIX, shape=(None, 50), days_window=days_window, tensor_from_file=tff, time_series_limit=0)
                name2tensormap[f'{needed_name}_newest'] = TensorMap(f'{needed_name}_newest', Interpretation.SURVIVAL_CURVE, path_prefix=PARTNERS_PREFIX, shape=(50,), days_window=days_window, tensor_from_file=tff)

    return name2tensormap


def _cardiac_surgery_str2date(input_date: str) -> datetime.datetime:
    return datetime.datetime.strptime(input_date, "%d%b%Y")


def build_cardiac_surgery_outcome_tensor_from_file(
    file_name: str,
    patient_column: str = "medrecn",
    outcome_column: str = "mtopd",
    start_column: str = "surgdt",
    delimiter: str = ",",
    day_window: int = 30,
) -> Callable:
    """Build a tensor_from_file function for outcomes given CSV of patients.

    The tensor_from_file function returned here should be used
    with CATEGORICAL TensorMaps to classify patients by disease state.

    :param file_name: CSV or TSV file with header of patient IDs (MRNs) dates of enrollment and dates of diagnosis
    :param patient_column: The header name of the column of patient ids
    :param diagnosis_date_column: The header name of the column of disease diagnosis dates
    :param start_column: The header name of the column of surgery dates
    :param delimiter: The delimiter separating columns of the TSV or CSV
    :return: The tensor_from_file function to provide to TensorMap constructors
    """
    error = None
    try:
        with open(file_name, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter=delimiter)
            header = next(reader)
            patient_index = header.index(patient_column)
            date_index = header.index(start_column)
            outcome_index = header.index(outcome_column)
            date_surg_table = {}
            patient_table = {}
            for row in reader:
                try:
                    patient_key = int(row[patient_index])
                    patient_table[patient_key] = int(row[outcome_index])
                    date_surg_table[patient_key] = _cardiac_surgery_str2date(
                        row[date_index],
                    )
                    if len(patient_table) % 1000 == 0:
                        logging.debug(f"Processed: {len(patient_table)} patient rows.")
                except ValueError as e:
                    logging.warning(f"val err {e}")

    except FileNotFoundError as e:
        error = e

    logging.info(
        f"Done processing {outcome_column}. Got {len(patient_table)} patient rows.",
    )

    def tensor_from_file(tm: TensorMap, hd5: h5py.File, dependents=None):
        if error:
            raise error

        ecg_dates = _get_ecg_dates(tm, hd5)
        dynamic, shape = _is_dynamic_shape(tm, len(ecg_dates))
        categorical_data = np.zeros(shape, dtype=np.float32)
        for i, ecg_date in enumerate(ecg_dates):
            file_split = os.path.basename(hd5.filename).split("-")
            mrn = file_split[0]
            mrn_int = int(mrn)

            if mrn_int not in patient_table:
                raise KeyError(f"MRN not in STS outcomes CSV")

            path = _make_hd5_path(tm, ecg_date, 'acquisitiondate')
            ecg_date = _partners_str2date(
                decompress_data(
                    data_compressed=hd5[path][()],
                    dtype=hd5[path].attrs["dtype"],
                ),
            )

            # Convert ecg_date from datetime.date to datetime.datetime
            ecg_date = datetime.datetime.combine(ecg_date, datetime.time.min)

            # If the date of surgery - date of ECG is > time window, skip it
            if date_surg_table[mrn_int] - ecg_date > datetime.timedelta(days=day_window):
                raise ValueError(f"ECG out of time window")

            categorical_data[(i, patient_table[mrn_int]) if dynamic else (patient_table[mrn_int],)] = 1.0
        return categorical_data
    return tensor_from_file


def build_cardiac_surgery_tensor_maps(
    needed_tensor_maps: List[str],
) -> Dict[str, TensorMap]:
    name2tensormap: Dict[str, TensorMap] = {}
    outcome2column = {
        "death": "mtopd",
        "stroke": "cnstrokp",
        "renal_failure": "crenfail",
        "prolonged_ventilation": "crenfail",
        "dsw_infection": "deepsterninf",
        "reoperation": "reop",
        "any_morbidity": "anymorbidity",
        "long_stay": "llos",
    }

    for outcome in outcome2column:
        name = f"outcome_{outcome}"
        if name in needed_tensor_maps:
            tensor_from_file_fxn = build_cardiac_surgery_outcome_tensor_from_file(
                file_name=CARDIAC_SURGERY_OUTCOMES_CSV,
                outcome_column=outcome2column[outcome],
                day_window=30,
            )
            name2tensormap[name] = TensorMap(
                name,
                Interpretation.CATEGORICAL,
                path_prefix=PARTNERS_PREFIX,
                channel_map=_outcome_channels(outcome),
                tensor_from_file=tensor_from_file_fxn,
                time_series_limit=0,
            )
            name2tensormap[f'{name}_newest'] = TensorMap(
                f'{name}_newest',
                Interpretation.CATEGORICAL,
                path_prefix=PARTNERS_PREFIX,
                channel_map=_outcome_channels(outcome),
                tensor_from_file=tensor_from_file_fxn,
            )
    return name2tensormap
