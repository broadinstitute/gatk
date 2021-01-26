# MRI-specific tensormaps
import logging
from typing import Dict, Tuple, Callable

import h5py
import numpy as np
from tensorflow.keras.utils import to_categorical
import cv2
import blosc

from ml4h.metrics import weighted_crossentropy
from ml4h.normalizer import ZeroMeanStd1, Standardize, NonZeroNormalize, TopKNormalize, ImagenetNormalizeTorch
from ml4h.TensorMap import TensorMap, Interpretation, make_range_validator
from ml4h.tensormap.ukb.demographics import is_genetic_man, is_genetic_woman
from ml4h.defines import MRI_TO_SEGMENT, MRI_SEGMENTED, MRI_SEGMENTED_CHANNEL_MAP, MRI_FRAMES, MRI_LVOT_SEGMENTED_CHANNEL_MAP, \
    MRI_LAX_2CH_SEGMENTED_CHANNEL_MAP
from ml4h.tensormap.general import get_tensor_at_first_date, normalized_first_date, pad_or_crop_array_to_shape
from ml4h.defines import MRI_LAX_3CH_SEGMENTED_CHANNEL_MAP, MRI_LAX_4CH_SEGMENTED_CHANNEL_MAP, MRI_SAX_SEGMENTED_CHANNEL_MAP, MRI_AO_SEGMENTED_CHANNEL_MAP, MRI_LIVER_SEGMENTED_CHANNEL_MAP


def _slice_subset_tensor(
    tensor_key,
    start,
    stop,
    step=1,
    dependent_key=None,
    pad_shape=None,
    dtype_override=None,
    allow_channels=True,
    flip_swap=False,
    swap_axes=-1,
):
    def _slice_subset_tensor_from_file(
        tm: TensorMap,
        hd5: h5py.File,
        dependents=None,
    ):
        if dtype_override is not None:
            big_tensor = get_tensor_at_first_date(
                hd5, tm.path_prefix,
                tensor_key,
            )
        else:
            big_tensor = get_tensor_at_first_date(
                hd5, tm.path_prefix,
                tensor_key,
            )

        if flip_swap:
            big_tensor = np.flip(np.swapaxes(big_tensor, 0, swap_axes))

        if pad_shape is not None:
            big_tensor = pad_or_crop_array_to_shape(pad_shape, big_tensor)

        if allow_channels and tm.shape[-1] < (stop - start) // step:
            tensor = big_tensor[..., np.arange(start, stop, step), :]
        else:
            tensor = big_tensor[..., np.arange(start, stop, step)]

        if dependent_key is not None:
            label_tensor = np.array(
                hd5[dependent_key][..., start:stop],
                dtype=np.float32,
            )
            dependents[tm.dependent_map] = to_categorical(
                label_tensor, tm.dependent_map.shape[-1],
            )
        return tensor

    return _slice_subset_tensor_from_file


def _random_slice_tensor(tensor_key, dependent_key=None):
    def _random_slice_tensor_from_file(
        tm: TensorMap,
        hd5: h5py.File,
        dependents=None,
    ):
        big_tensor = get_tensor_at_first_date(hd5, tm.path_prefix, tensor_key)
        cur_slice = np.random.choice(range(big_tensor.shape[-1]))
        tensor = np.zeros(tm.shape, dtype=np.float32)
        tensor[..., 0] = big_tensor[..., cur_slice]
        if dependent_key is not None:
            dependents[tm.dependent_map] = np.zeros(
                tm.dependent_map.shape,
                dtype=np.float32,
            )
            label_tensor = np.array(
                hd5[dependent_key][..., cur_slice],
                dtype=np.float32,
            )
            dependents[tm.dependent_map][:, :, :] = to_categorical(
                label_tensor, tm.dependent_map.shape[-1],
            )
        return tensor

    return _random_slice_tensor_from_file


def _mask_from_file(tm: TensorMap, hd5: h5py.File, dependents=None):
    original = get_tensor_at_first_date(hd5, tm.path_prefix, tm.name)
    reshaped = pad_or_crop_array_to_shape(tm.shape, original)
    tensor = to_categorical(reshaped[..., 0], tm.shape[-1])
    return tensor


def _mask_subset_tensor(tensor_key, start, stop, step=1, pad_shape=None):
    slice_subset_tensor_from_file = _slice_subset_tensor(
        tensor_key,
        start,
        stop,
        step=step,
        pad_shape=pad_shape,
        dtype_override='float_array',
    )

    def mask_subset_from_file(tm: TensorMap, hd5: h5py.File, dependents=None):
        original = slice_subset_tensor_from_file(tm, hd5, dependents)
        tensor = to_categorical(original[..., 0], tm.shape[-1])
        return tensor

    return mask_subset_from_file


def _combined_subset_tensor(
    tensor_keys,
    start,
    stop,
    step=1,
    pad_shape=None,
    flip_swap=False,
):
    slice_subsets = [
        _slice_subset_tensor(
            k,
            start,
            stop,
            step=step,
            pad_shape=pad_shape,
            allow_channels=False,
            flip_swap=flip_swap,
        ) for k in tensor_keys
    ]

    def mask_subset_from_file(tm: TensorMap, hd5: h5py.File, dependents=None):
        tensor = np.zeros(tm.shape, dtype=np.float32)
        for i, slice_subset_tensor_from_file in enumerate(slice_subsets):
            tensor[..., i] = slice_subset_tensor_from_file(tm, hd5, dependents)
        return tensor

    return mask_subset_from_file


def mri_tensor_2d(hd5, name):
    """
    Returns MRI image annotation tensors as 2-D numpy arrays. Useful for annotations that may vary from slice to slice
    """
    if isinstance(hd5[name], h5py.Group):
        nslices = len(hd5[name]) // MRI_FRAMES
        for ann in hd5[name]:
            ann_shape = hd5[name][ann].shape
            break
        shape = (ann_shape[0], nslices)
        arr = np.zeros(shape)
        t = 0
        s = 0
        for k in sorted(hd5[name], key=int):
            t += 1
            if t == MRI_FRAMES:
                arr[:, s] = hd5[name][k]
                s += 1
                t = 0
    elif isinstance(hd5[name], h5py.Dataset):
        nslices = 1
        shape = (hd5[name].shape[0], nslices)
        arr = np.zeros(shape)
        arr[:, 0] = hd5[name]
    else:
        raise ValueError(f'{name} is neither a HD5 Group nor a HD5 dataset')
    return arr


def _make_mri_series_orientation_and_position_from_file(
    population_normalize=None,
):
    def mri_series_orientation_and_position(tm, hd5):
        if len(tm.shape) < 2:
            tensor = np.array(hd5[tm.name], dtype=np.float32)
        else:
            arr = mri_tensor_2d(hd5, tm.name)
            tensor = np.array(arr, dtype=np.float32)
        if population_normalize is not None:
            tensor /= population_normalize
        return tensor

    return mri_series_orientation_and_position


def _slice_tensor(tensor_key, slice_index):
    def _slice_tensor_from_file(tm, hd5, dependents={}):
        if tm.shape[-1] == 1:
            t = pad_or_crop_array_to_shape(
                tm.shape[:-1],
                np.array(hd5[tensor_key][..., slice_index], dtype=np.float32),
            )
            tensor = np.expand_dims(t, axis=-1)
        else:
            tensor = pad_or_crop_array_to_shape(
                tm.shape,
                np.array(hd5[tensor_key][..., slice_index], dtype=np.float32),
            )
        return tensor

    return _slice_tensor_from_file


def _segmented_dicom_slices(dicom_key_prefix, path_prefix='ukb_cardiac_mri', step=1, total_slices=50):
    def _segmented_dicom_tensor_from_file(tm, hd5, dependents={}):
        tensor = np.zeros(tm.shape, dtype=np.float32)
        if tm.axes() == 3 or path_prefix == 'ukb_liver_mri':
            categorical_index_slice = get_tensor_at_first_date(hd5, path_prefix, f'{dicom_key_prefix}1')
            categorical_one_hot = to_categorical(categorical_index_slice, len(tm.channel_map))
            tensor[..., :] = pad_or_crop_array_to_shape(tensor[..., :].shape, categorical_one_hot)
        elif tm.axes() == 4:
            tensor_index = 0
            for i in range(0, total_slices, step):
                categorical_index_slice = get_tensor_at_first_date(hd5, path_prefix, f'{dicom_key_prefix}{i+1}')
                categorical_one_hot = to_categorical(categorical_index_slice, len(tm.channel_map))
                tensor[..., tensor_index, :] = pad_or_crop_array_to_shape(tensor[..., tensor_index, :].shape, categorical_one_hot)
                tensor_index += 1
                if tensor_index >= tensor.shape[-2]:
                    break
        else:
            raise ValueError(f'No method to get segmented slices for TensorMap: {tm}')
        return tensor
    return _segmented_dicom_tensor_from_file


def _mri_slice_blackout_tensor_from_file(tm, hd5, dependents={}):
    cur_slice = np.random.choice(list(hd5[MRI_TO_SEGMENT].keys()))
    tensor = np.zeros(tm.shape, dtype=np.float32)
    dependents[tm.dependent_map] = np.zeros(
        tm.dependent_map.shape,
        dtype=np.float32,
    )
    tensor[:, :, 0] = np.array(
        hd5[MRI_TO_SEGMENT][cur_slice],
        dtype=np.float32,
    )
    label_tensor = np.array(hd5[MRI_SEGMENTED][cur_slice], dtype=np.float32)
    dependents[tm.dependent_map][:, :, :] = to_categorical(
        label_tensor, tm.dependent_map.shape[-1],
    )
    tensor[:, :, 0] *= np.not_equal(label_tensor, 0, dtype=np.float32)
    return tm.zero_mean_std1(tensor)


t2_flair_sag_p2_1mm_fs_ellip_pf78_1 = TensorMap(
    't2_flair_sag_p2_1mm_fs_ellip_pf78_1',
    shape=(256, 256, 192),
    path_prefix='ukb_brain_mri',
    tensor_from_file=normalized_first_date,
    normalization=ZeroMeanStd1(),
)
t2_flair_sag_p2_1mm_fs_ellip_pf78_2 = TensorMap(
    't2_flair_sag_p2_1mm_fs_ellip_pf78_2',
    shape=(256, 256, 192),
    path_prefix='ukb_brain_mri',
    tensor_from_file=normalized_first_date,
    normalization=ZeroMeanStd1(),
)
t2_flair_slice_1 = TensorMap(
    't2_flair_slice_1',
    shape=(256, 256, 1),
    path_prefix='ukb_brain_mri',
    tensor_from_file=_random_slice_tensor(
        't2_flair_sag_p2_1mm_fs_ellip_pf78_1',
    ),
    normalization=ZeroMeanStd1(),
)
t2_flair_slice_2 = TensorMap(
    't2_flair_slice_2',
    shape=(256, 256, 1),
    path_prefix='ukb_brain_mri',
    tensor_from_file=_random_slice_tensor(
        't2_flair_sag_p2_1mm_fs_ellip_pf78_2',
    ),
    normalization=ZeroMeanStd1(),
)
t1_p2_1mm_fov256_sag_ti_880_1 = TensorMap(
    't1_p2_1mm_fov256_sag_ti_880_1',
    shape=(256, 256, 208),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)
t1_p2_1mm_fov256_sag_ti_880_2 = TensorMap(
    't1_p2_1mm_fov256_sag_ti_880_2',
    shape=(256, 256, 208),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)
t1_dicom_30_slices = TensorMap(
    't1_dicom_30_slices',
    shape=(192, 256, 30),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        't1_p2_1mm_fov256_sag_ti_880_1',
        130,
        190,
        2,
        pad_shape=(192, 256, 256),
        flip_swap=True,
    ),
)
t2_dicom_30_slices = TensorMap(
    't2_dicom_30_slices',
    shape=(192, 256, 30),
    path_prefix='ukb_brain_mri/',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        't2_flair_sag_p2_1mm_fs_ellip_pf78_1',
        130,
        190,
        2,
        pad_shape=(192, 256, 256),
        flip_swap=True,
    ),
)

t1_slice_1 = TensorMap(
    't1_slice_1',
    shape=(256, 256, 1),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_random_slice_tensor('t1_p2_1mm_fov256_sag_ti_880_1'),
)
t1_slice_2 = TensorMap(
    't1_slice_2',
    shape=(256, 256, 1),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_random_slice_tensor('t1_p2_1mm_fov256_sag_ti_880_2'),
)
t1_20_slices_1 = TensorMap(
    't1_20_slices_1',
    shape=(256, 256, 20),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        't1_p2_1mm_fov256_sag_ti_880_1', 94,
        114,
    ),
)
t1_20_slices_2 = TensorMap(
    't1_20_slices_2',
    shape=(256, 256, 20),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        't1_p2_1mm_fov256_sag_ti_880_2', 94,
        114,
    ),
)
t2_20_slices_1 = TensorMap(
    't2_20_slices_1',
    shape=(256, 256, 20),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        't2_flair_sag_p2_1mm_fs_ellip_pf78_1', 86, 106,
    ),
)
t2_20_slices_2 = TensorMap(
    't2_20_slices_2',
    shape=(256, 256, 20),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        't2_flair_sag_p2_1mm_fs_ellip_pf78_2', 86, 106,
    ),
)
t1_40_slices_1 = TensorMap(
    't1_40_slices_1',
    shape=(256, 256, 40),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        't1_p2_1mm_fov256_sag_ti_880_1', 64,
        144, 2,
    ),
)
t2_40_slices_1 = TensorMap(
    't2_40_slices_1',
    shape=(256, 256, 40),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        't2_flair_sag_p2_1mm_fs_ellip_pf78_1', 56, 136, 2,
    ),
)
sos_te1 = TensorMap(
    'SOS_TE1',
    shape=(256, 288, 48),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)
sos_te2 = TensorMap(
    'SOS_TE2',
    shape=(256, 288, 48),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)
swi = TensorMap(
    'SWI',
    shape=(256, 288, 48),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)
swi_total_mag = TensorMap(
    'SWI_TOTAL_MAG',
    shape=(256, 288, 48),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)
swi_total_mag_te2_orig = TensorMap(
    'SWI_TOTAL_MAG_TE2_orig',
    shape=(256, 288, 48),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)
swi_total_mag_orig = TensorMap(
    'SWI_TOTAL_MAG_orig',
    shape=(256, 288, 48),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)
t2star = TensorMap(
    'T2star',
    shape=(256, 288, 48),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)
brain_mask_normed = TensorMap(
    'brain_mask_normed',
    shape=(256, 288, 48),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)

filtered_phase = TensorMap(
    'filtered_phase',
    shape=(256, 288, 48),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)
swi_to_t1_40_slices = TensorMap(
    'swi_to_t1_40_slices',
    shape=(173, 231, 40),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor('SWI_TOTAL_MAG_to_T1', 60, 140, 2),
)
t2star_to_t1_40_slices = TensorMap(
    't2star_to_t1_40_slices',
    shape=(173, 231, 40),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor('T2star_to_T1', 60, 140, 2),
)

t1 = TensorMap(
    'T1',
    shape=(192, 256, 256, 1),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)
t1_brain = TensorMap(
    'T1_brain',
    shape=(192, 256, 256, 1),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)
t1_brain_30_slices = TensorMap(
    't1_brain_30_slices',
    shape=(192, 256, 30),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        'T1_brain',
        66,
        126,
        2,
        pad_shape=(192, 256, 256),
    ),
)
t1_30_slices = TensorMap(
    't1_30_slices',
    shape=(192, 256, 30),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        'T1', 90, 150, 2, pad_shape=(192, 256, 256),
    ),
)
t1_30_slices_4d = TensorMap(
    't1_30_slices_4d',
    shape=(192, 256, 30, 1),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        'T1',
        90,
        150,
        2,
        pad_shape=(192, 256, 256, 1),
    ),
)
t1_30_slices_fs = TensorMap(
    't1_30_slices',
    shape=(192, 256, 30),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        'T1',
        90,
        150,
        2,
        pad_shape=(192, 256, 256),
        flip_swap=True,
        swap_axes=1,
    ),
)
t1_30_slices_4d_fs = TensorMap(
    't1_30_slices_4d',
    shape=(192, 256, 30, 1),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        'T1',
        90,
        150,
        2,
        pad_shape=(192, 256, 256, 1),
        flip_swap=True,
        swap_axes=1,
    ),
)

t1_brain_to_mni = TensorMap(
    'T1_brain_to_MNI',
    shape=(192, 256, 256, 1),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)
t1_fast_t1_brain_bias = TensorMap(
    'T1_fast_T1_brain_bias',
    shape=(192, 256, 256, 1),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)

t2_flair = TensorMap(
    'T2_FLAIR',
    shape=(192, 256, 256, 1),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)
t2_flair_brain = TensorMap(
    'T2_FLAIR_brain',
    shape=(192, 256, 256, 1),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)
t2_flair_brain_30_slices = TensorMap(
    't2_flair_brain_30_slices',
    shape=(192, 256, 30),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        'T2_FLAIR_brain',
        66,
        126,
        2,
        pad_shape=(192, 256, 256),
    ),
)
t2_flair_30_slices = TensorMap(
    't2_flair_30_slices',
    shape=(192, 256, 30),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        'T2_FLAIR',
        90,
        150,
        2,
        pad_shape=(192, 256, 256),
    ),
)
t2_flair_30_slices_4d = TensorMap(
    't2_flair_30_slices_4d',
    shape=(192, 256, 30, 1),
    path_prefix='ukb_brain_mri',
    tensor_from_file=_slice_subset_tensor(
        'T2_FLAIR',
        90,
        150,
        2,
        pad_shape=(192, 256, 256, 1),
    ),
    normalization=ZeroMeanStd1(),
)
t2_flair_30_slices_fs = TensorMap(
    't2_flair_30_slices',
    shape=(192, 256, 30),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        'T2_FLAIR',
        90,
        150,
        2,
        pad_shape=(192, 256, 256),
        flip_swap=True,
        swap_axes=1,
    ),
)
t2_flair_30_slices_4d_fs = TensorMap(
    't2_flair_30_slices_4d',
    shape=(192, 256, 30, 1),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_subset_tensor(
        'T2_FLAIR',
        90,
        150,
        2,
        pad_shape=(192, 256, 256, 1),
        flip_swap=True,
        swap_axes=1,
    ),
)
t2_flair_unbiased_brain = TensorMap(
    'T2_FLAIR_unbiased_brain',
    shape=(192, 256, 256, 1),
    path_prefix='ukb_brain_mri',
    normalization=ZeroMeanStd1(),
    tensor_from_file=normalized_first_date,
)

swi_brain_mask = TensorMap(
    'SWI_brain_mask',
    Interpretation.CATEGORICAL,
    shape=(256, 288, 48, 2),
    path_prefix='ukb_brain_mri',
    tensor_from_file=_mask_from_file,
    channel_map={
        'not_brain': 0,
        'brain': 1,
    },
)
t1_brain_mask = TensorMap(
    'T1_brain_mask',
    Interpretation.CATEGORICAL,
    shape=(192, 256, 256, 2),
    path_prefix='ukb_brain_mri',
    tensor_from_file=_mask_from_file,
    channel_map={
        'not_brain': 0,
        'brain': 1,
    },
)
t1_seg = TensorMap(
    'T1_fast_T1_brain_seg',
    Interpretation.CATEGORICAL,
    shape=(192, 256, 256, 4),
    path_prefix='ukb_brain_mri',
    tensor_from_file=_mask_from_file,
    channel_map={
        'not_brain_tissue': 0,
        'csf': 1,
        'grey': 2,
        'white': 3,
    },
)
t1_seg_30_slices = TensorMap(
    'T1_fast_T1_brain_seg_30_slices',
    Interpretation.CATEGORICAL,
    shape=(192, 256, 30, 4),
    path_prefix='ukb_brain_mri',
    tensor_from_file=_mask_subset_tensor(
        'T1_fast_T1_brain_seg',
        90,
        150,
        2,
        pad_shape=(192, 256, 256, 1),
    ),
    channel_map={
        'not_brain_tissue': 0,
        'csf': 1,
        'grey': 2,
        'white': 3,
    },
)
t1_brain_mask_30_slices = TensorMap(
    'T1_brain_mask_30_slices',
    Interpretation.CATEGORICAL,
    shape=(192, 256, 30, 2),
    path_prefix='ukb_brain_mri',
    tensor_from_file=_mask_subset_tensor(
        'T1_brain_mask',
        90,
        150,
        2,
        pad_shape=(192, 256, 256, 1),
    ),
    channel_map={
        'not_brain': 0,
        'brain': 1,
    },
)
lesions = TensorMap(
    'lesions_final_mask',
    Interpretation.CATEGORICAL,
    shape=(192, 256, 256, 2),
    path_prefix='ukb_brain_mri',
    tensor_from_file=_mask_from_file,
    channel_map={
        'not_lesion': 0,
        'lesion': 1,
    },
    loss=weighted_crossentropy([0.01, 10.0], 'lesion'),
)

t1_and_t2_flair_30_slices = TensorMap(
    't1_and_t2_flair_30_slices',
    Interpretation.CONTINUOUS,
    shape=(192, 256, 30, 2),
    path_prefix='ukb_brain_mri',
    tensor_from_file=_combined_subset_tensor(
        ['T1', 'T2_FLAIR'],
        90,
        150,
        2,
        pad_shape=(192, 256, 256),
    ),
    normalization=ZeroMeanStd1(),
)
_dicom_keys = [
    't1_p2_1mm_fov256_sag_ti_880_1', 't2_flair_sag_p2_1mm_fs_ellip_pf78_1',
]
t1_t2_dicom_30_slices = TensorMap(
    't1_t2_dicom_30_slices',
    Interpretation.CONTINUOUS,
    shape=(192, 256, 30, 2),
    path_prefix='ukb_brain_mri',
    tensor_from_file=_combined_subset_tensor(
        _dicom_keys,
        130,
        190,
        2,
        pad_shape=(192, 256, 256),
        flip_swap=True,
    ),
    normalization=ZeroMeanStd1(),
)
t1_t2_dicom_50_slices = TensorMap(
    't1_t2_dicom_50_slices',
    Interpretation.CONTINUOUS,
    shape=(192, 256, 50, 2),
    path_prefix='ukb_brain_mri',
    tensor_from_file=_combined_subset_tensor(
        _dicom_keys,
        100,
        200,
        2,
        pad_shape=(192, 256, 256),
        flip_swap=True,
    ),
    normalization=ZeroMeanStd1(),
)

mri_slice_blackout_segmented_weighted = TensorMap(
    'mri_slice_segmented',
    Interpretation.CATEGORICAL,
    shape=(256, 256, 3),
    channel_map=MRI_SEGMENTED_CHANNEL_MAP,
    loss=weighted_crossentropy(
        [0.1, 25.0, 25.0],
        'mri_slice_blackout_segmented',
    ),
)
mri_slice_blackout = TensorMap(
    'mri_slice_blackout',
    Interpretation.CONTINUOUS,
    shape=(256, 256, 1),
    tensor_from_file=_mri_slice_blackout_tensor_from_file,
    dependent_map=mri_slice_blackout_segmented_weighted,
)

mri_patient_orientation_cine_segmented_lax_2ch = TensorMap(
    'mri_patient_orientation_cine_segmented_lax_2ch',
    Interpretation.CONTINUOUS,
    shape=(6,),
    path_prefix='mri_orientation',
    tensor_from_file=_make_mri_series_orientation_and_position_from_file(),
)
mri_patient_orientation_cine_segmented_lax_3ch = TensorMap(
    'mri_patient_orientation_cine_segmented_lax_3ch',
    Interpretation.CONTINUOUS,
    shape=(6,),
    path_prefix='mri_orientation',
    tensor_from_file=_make_mri_series_orientation_and_position_from_file(),
)
mri_patient_orientation_cine_segmented_lax_4ch = TensorMap(
    'mri_patient_orientation_cine_segmented_lax_4ch',
    Interpretation.CONTINUOUS,
    shape=(6,),
    path_prefix='mri_orientation',
    tensor_from_file=_make_mri_series_orientation_and_position_from_file(),
)
mri_patient_orientation_cine_segmented_sax_b1 = TensorMap(
    'mri_patient_orientation_cine_segmented_sax_b1',
    Interpretation.CONTINUOUS,
    shape=(6,),
    path_prefix='mri_orientation',
    tensor_from_file=_make_mri_series_orientation_and_position_from_file(),
)
mri_patient_orientation_cine_segmented_sax_inlinevf = TensorMap(
    'mri_patient_orientation_cine_segmented_sax_inlinevf',
    Interpretation.CONTINUOUS,
    shape=(6, 750),
    path_prefix='mri_orientation',
    tensor_from_file=_make_mri_series_orientation_and_position_from_file(),
)
mri_patient_position_cine_segmented_lax_2ch = TensorMap(
    'mri_patient_position_cine_segmented_lax_2ch',
    Interpretation.CONTINUOUS,
    shape=(3,),
    path_prefix='mri_position',
    tensor_from_file=_make_mri_series_orientation_and_position_from_file(),
)
mri_patient_position_cine_segmented_lax_3ch = TensorMap(
    'mri_patient_position_cine_segmented_lax_3ch',
    Interpretation.CONTINUOUS,
    shape=(3,),
    path_prefix='mri_position',
    tensor_from_file=_make_mri_series_orientation_and_position_from_file(),
)
mri_patient_position_cine_segmented_lax_4ch = TensorMap(
    'mri_patient_position_cine_segmented_lax_4ch',
    Interpretation.CONTINUOUS,
    shape=(3,),
    path_prefix='mri_position',
    tensor_from_file=_make_mri_series_orientation_and_position_from_file(),
)
mri_patient_position_cine_segmented_sax_b1 = TensorMap(
    'mri_patient_position_cine_segmented_sax_b1',
    Interpretation.CONTINUOUS,
    shape=(3,),
    path_prefix='mri_position',
    tensor_from_file=_make_mri_series_orientation_and_position_from_file(),
)
mri_patient_position_cine_segmented_sax_inlinevf = TensorMap(
    'mri_patient_position_cine_segmented_sax_inlinevf',
    Interpretation.CONTINUOUS,
    shape=(3, 750),
    path_prefix='mri_position',
    tensor_from_file=_make_mri_series_orientation_and_position_from_file(),
)

lax_4ch_diastole_slice0_3d = TensorMap(
    'lax_4ch_diastole_slice0_3d',
    Interpretation.CONTINUOUS,
    shape=(200, 160, 1),
    loss='logcosh',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_tensor(
        'ukb_cardiac_mri/cine_segmented_lax_4ch/instance_0', 0,
    ),
)
lax_4ch_diastole_slice0_224_3d = TensorMap(
    'lax_4ch_diastole_slice0_224_3d', Interpretation.CONTINUOUS, shape=(160, 224, 1), loss='logcosh',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_tensor('ukb_cardiac_mri/cine_segmented_lax_4ch/instance_0', 0),
)
lax_4ch_diastole_slice0_256_3d = TensorMap(
    'lax_4ch_diastole_slice0_256_3d', Interpretation.CONTINUOUS, shape=(192, 256, 1),
    normalization=ZeroMeanStd1(), tensor_from_file=_slice_tensor('ukb_cardiac_mri/cine_segmented_lax_4ch/instance_0', 0),
)
lax_2ch_diastole_slice0_3d = TensorMap(
    'lax_2ch_diastole_slice0_3d',
    Interpretation.CONTINUOUS,
    shape=(200, 160, 1),
    loss='logcosh',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_tensor(
        'ukb_cardiac_mri/cine_segmented_lax_2ch/instance_0', 0,
    ),
)
lax_3ch_diastole_slice0_3d = TensorMap(
    'lax_3ch_diastole_slice0_3d',
    Interpretation.CONTINUOUS,
    shape=(200, 160, 1),
    loss='logcosh',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_tensor(
        'ukb_cardiac_mri/cine_segmented_lax_3ch/instance_0', 0,
    ),
)
cine_segmented_ao_dist_slice0_3d = TensorMap(
    'cine_segmented_ao_dist_slice0_3d',
    Interpretation.CONTINUOUS,
    shape=(256, 256, 1),
    loss='logcosh',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_tensor(
        'ukb_cardiac_mri/cine_segmented_ao_dist/instance_0', 0,
    ),
)
lax_4ch_diastole_slice0 = TensorMap(
    'lax_4ch_diastole_slice0',
    Interpretation.CONTINUOUS,
    shape=(256, 256),
    loss='logcosh',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_tensor(
        'ukb_cardiac_mri/cine_segmented_lax_4ch/instance_0', 0,
    ),
)
cine_segmented_ao_dist_slice0 = TensorMap(
    'cine_segmented_ao_dist_slice0',
    Interpretation.CONTINUOUS,
    shape=(256, 256),
    loss='logcosh',
    normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_tensor(
        'ukb_cardiac_mri/cine_segmented_ao_dist/instance_0', 0,
    ),
)
aorta_diastole_slice0_3d = TensorMap(
    'aorta_diastole_slice0_3d', Interpretation.CONTINUOUS, shape=(192, 256, 1), loss='logcosh',
    normalization=ZeroMeanStd1(), tensor_from_file=_slice_tensor('ukb_cardiac_mri/cine_segmented_ao_dist/instance_0', 0),
)
cine_segmented_lvot_slice0_3d = TensorMap(
    'cine_segmented_lvot_slice0_3d', Interpretation.CONTINUOUS, shape=(208, 160, 1), loss='logcosh',
    normalization=ZeroMeanStd1(), tensor_from_file=_slice_tensor('ukb_cardiac_mri/cine_segmented_lvot/instance_0', 0),
)


def _pad_crop_tensor(tm, hd5, dependents={}):
    return pad_or_crop_array_to_shape(
        tm.shape,
        np.array(
            tm.hd5_first_dataset_in_group(hd5, tm.hd5_key_guess()),
            dtype=np.float32,
        ),
    )


cine_lax_3ch_192 = TensorMap(
    'cine_segmented_lax_3ch',
    Interpretation.CONTINUOUS,
    shape=(192, 192, 50),
    path_prefix='ukb_cardiac_mri',
    tensor_from_file=_pad_crop_tensor,
    normalization=ZeroMeanStd1(),
)
cine_lax_3ch_160_1 = TensorMap(
    'cine_segmented_lax_3ch',
    Interpretation.CONTINUOUS,
    shape=(160, 160, 50, 1),
    path_prefix='ukb_cardiac_mri',
    tensor_from_file=_pad_crop_tensor,
    normalization=ZeroMeanStd1(),
)
cine_lax_3ch_192_160_1 = TensorMap(
    'cine_segmented_lax_3ch',
    Interpretation.CONTINUOUS,
    shape=(192, 160, 50, 1),
    path_prefix='ukb_cardiac_mri',
    tensor_from_file=_pad_crop_tensor,
    normalization=ZeroMeanStd1(),
)
cine_ao_dist_4d = TensorMap(
    'cine_ao_dist_4d',
    Interpretation.CONTINUOUS,
    shape=(160, 192, 100, 1),
    path_prefix='ukb_cardiac_mri',
    tensor_from_file=_pad_crop_tensor,
    normalization=ZeroMeanStd1(),
)
cine_lax_4ch_192 = TensorMap(
    'cine_segmented_lax_3ch',
    Interpretation.CONTINUOUS,
    shape=(192, 192, 50),
    path_prefix='ukb_cardiac_mri',
    tensor_from_file=_pad_crop_tensor,
    normalization=ZeroMeanStd1(),
)
cine_lax_4ch_192_1 = TensorMap(
    'cine_segmented_lax_3ch',
    Interpretation.CONTINUOUS,
    shape=(192, 192, 50, 1),
    path_prefix='ukb_cardiac_mri',
    tensor_from_file=_pad_crop_tensor,
    normalization=ZeroMeanStd1(),
)
cine_sax_b6_192 = TensorMap(
    'cine_segmented_sax_b6',
    Interpretation.CONTINUOUS,
    shape=(192, 192, 50),
    path_prefix='ukb_cardiac_mri',
    tensor_from_file=_pad_crop_tensor,
    normalization=ZeroMeanStd1(),
)
cine_sax_b6_192_1 = TensorMap(
    'cine_segmented_sax_b6',
    Interpretation.CONTINUOUS,
    shape=(192, 192, 50, 1),
    path_prefix='ukb_cardiac_mri',
    tensor_from_file=_pad_crop_tensor,
    normalization=ZeroMeanStd1(),
)
flow_250_tp_aov_bh_epat = TensorMap(
    'flow_250_tp_aov_bh_epat', Interpretation.CONTINUOUS, shape=(192, 192, 30), path_prefix='ukb_cardiac_mri',
    tensor_from_file=_slice_subset_tensor('flow_250_tp_aov_bh_epat@c', 0, 30, pad_shape=(192, 192, 30)),
    normalization=ZeroMeanStd1(),
)
flow_250_tp_aov_bh_epat_mag = TensorMap(
    'flow_250_tp_aov_bh_epat_mag', Interpretation.CONTINUOUS, shape=(192, 192, 30), path_prefix='ukb_cardiac_mri',
    tensor_from_file=_slice_subset_tensor('flow_250_tp_aov_bh_epat@c_mag', 0, 30, pad_shape=(192, 192, 30)),
    normalization=ZeroMeanStd1(),
)
flow_250_tp_aov_bh_epat_p = TensorMap(
    'flow_250_tp_aov_bh_epat_p', Interpretation.CONTINUOUS, shape=(192, 192, 30), path_prefix='ukb_cardiac_mri',
    tensor_from_file=_slice_subset_tensor('flow_250_tp_aov_bh_epat@c_p', 0, 30, pad_shape=(192, 192, 30)),
    normalization=ZeroMeanStd1(),
)
flow_250_tp_aov_bh_epat_4d = TensorMap(
    'flow_250_tp_aov_bh_epat', Interpretation.CONTINUOUS, shape=(192, 192, 30, 1), path_prefix='ukb_cardiac_mri',
    tensor_from_file=_slice_subset_tensor('flow_250_tp_aov_bh_epat@c', 0, 30, pad_shape=(192, 192, 30, 1)),
    normalization=ZeroMeanStd1(),
)
cine_lax_2ch_192_16_3 = TensorMap(
    'cine_lax_2ch_192_16_3', Interpretation.CONTINUOUS, shape=(192, 160, 16), path_prefix='ukb_cardiac_mri',
    tensor_from_file=_slice_subset_tensor('cine_segmented_lax_2ch', 0, 48, 3, pad_shape=(192, 160, 48)),
    normalization=ZeroMeanStd1(),
)
cine_lax_3ch_192_16_3 = TensorMap(
    'cine_lax_3ch_192_16_3', Interpretation.CONTINUOUS, shape=(192, 160, 16), path_prefix='ukb_cardiac_mri',
    tensor_from_file=_slice_subset_tensor('cine_segmented_lax_3ch', 0, 48, 3, pad_shape=(192, 160, 48)),
    normalization=ZeroMeanStd1(),
)
cine_lax_4ch_192_16_3 = TensorMap(
    'cine_lax_4ch_192_16_3', Interpretation.CONTINUOUS, shape=(192, 160, 16), path_prefix='ukb_cardiac_mri',
    tensor_from_file=_slice_subset_tensor('cine_segmented_lax_4ch', 0, 48, 3, pad_shape=(192, 160, 48)),
    normalization=ZeroMeanStd1(),
)
cine_lax_2ch_192_16_3_4d = TensorMap(
    'cine_lax_2ch_192_16_3_4d', Interpretation.CONTINUOUS, shape=(192, 160, 16, 1), path_prefix='ukb_cardiac_mri',
    tensor_from_file=_slice_subset_tensor('cine_segmented_lax_2ch', 0, 48, 3, pad_shape=(192, 160, 48, 1)),
    normalization=ZeroMeanStd1(),
)
cine_lax_3ch_192_16_3_4d = TensorMap(
    'cine_lax_3ch_192_16_3_4d', Interpretation.CONTINUOUS, shape=(192, 160, 16, 1), path_prefix='ukb_cardiac_mri',
    tensor_from_file=_slice_subset_tensor('cine_segmented_lax_3ch', 0, 48, 3, pad_shape=(192, 160, 48, 1)),
    normalization=ZeroMeanStd1(),
)
cine_lax_4ch_192_16_3_4d = TensorMap(
    'cine_lax_4ch_192_16_3_4d', Interpretation.CONTINUOUS, shape=(192, 160, 16, 1), path_prefix='ukb_cardiac_mri',
    tensor_from_file=_slice_subset_tensor('cine_segmented_lax_4ch', 0, 48, 3, pad_shape=(192, 160, 48, 1)),
    normalization=ZeroMeanStd1(),
)
cine_lax_4ch_224_16_3_4d = TensorMap(
    'cine_lax_4ch_224_16_3_4d', Interpretation.CONTINUOUS, shape=(160, 224, 16, 1), path_prefix='ukb_cardiac_mri',
    tensor_from_file=_slice_subset_tensor('cine_segmented_lax_4ch', 0, 48, 3, pad_shape=(160, 224, 48, 1)),
    normalization=ZeroMeanStd1(),
)
cine_lvot_208_16_3_4d = TensorMap(
    'cine_lvot_208_16_3_4d', Interpretation.CONTINUOUS, shape=(208, 192, 16, 1), path_prefix='ukb_cardiac_mri',
    tensor_from_file=_slice_subset_tensor('cine_segmented_lvot', 0, 48, 3, pad_shape=(208, 192, 48, 1)),
    normalization=ZeroMeanStd1(),
)
cine_lvot_192_16_3 = TensorMap(
    'cine_lvot_192_16_3', Interpretation.CONTINUOUS, shape=(192, 160, 16), path_prefix='ukb_cardiac_mri',
    tensor_from_file=_slice_subset_tensor('cine_segmented_lvot', 0, 48, 3, pad_shape=(192, 160, 48)),
    normalization=ZeroMeanStd1(),
)
cine_lvot_192_16_3_4d = TensorMap(
    'cine_lvot_192_16_3_4d', Interpretation.CONTINUOUS, shape=(192, 160, 16, 1), path_prefix='ukb_cardiac_mri',
    tensor_from_file=_slice_subset_tensor('cine_segmented_lvot', 0, 48, 3, pad_shape=(192, 160, 48, 1)),
    normalization=ZeroMeanStd1(),
)


lax_2ch_segmented_192_16_3 = TensorMap(
    'lax_2ch_segmented_192_16_3', Interpretation.CATEGORICAL, shape=(192, 160, 16, 13),
    tensor_from_file=_segmented_dicom_slices('cine_segmented_lax_2ch_annotated_', step=3),
    channel_map=MRI_LAX_2CH_SEGMENTED_CHANNEL_MAP,
)
lax_2ch_segmented_192 = TensorMap(
    'lax_2ch_segmented_192',
    Interpretation.CATEGORICAL,
    shape=(192, 192, 50, 6),
    tensor_from_file=_segmented_dicom_slices('cine_segmented_lax_2ch_annotated_'),
    channel_map=MRI_LAX_2CH_SEGMENTED_CHANNEL_MAP,
)
lax_3ch_segmented = TensorMap(
    'lax_3ch_segmented',
    Interpretation.CATEGORICAL,
    shape=(256, 256, 50, 6),
    tensor_from_file=_segmented_dicom_slices('cine_segmented_lax_3ch_annotated_'),
    channel_map=MRI_LAX_3CH_SEGMENTED_CHANNEL_MAP,
)
lax_3ch_segmented_192_160 = TensorMap(
    'lax_3ch_segmented_192_160',
    Interpretation.CATEGORICAL,
    shape=(192, 160, 50, 6),
    tensor_from_file=_segmented_dicom_slices(
        'cine_segmented_lax_3ch_annotated_',
    ),
    channel_map=MRI_LAX_3CH_SEGMENTED_CHANNEL_MAP,
)
lax_3ch_segmented_192_16_3 = TensorMap(
    'lax_3ch_segmented_192_16_3', Interpretation.CATEGORICAL, shape=(192, 160, 16, 6),
    tensor_from_file=_segmented_dicom_slices('cine_segmented_lax_3ch_annotated_', step=3),
    channel_map=MRI_LAX_3CH_SEGMENTED_CHANNEL_MAP,
)
lax_4ch_segmented = TensorMap(
    'lax_4ch_segmented',
    Interpretation.CATEGORICAL,
    shape=(256, 256, 50, len(MRI_LAX_4CH_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slices('cine_segmented_lax_4ch_annotated_'),
    channel_map=MRI_LAX_4CH_SEGMENTED_CHANNEL_MAP,
)
lax_4ch_segmented_224_16_3 = TensorMap(
    'lax_4ch_segmented_224_16_3', Interpretation.CATEGORICAL, shape=(160, 224, 16, len(MRI_LAX_4CH_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slices('cine_segmented_lax_4ch_annotated_', step=3),
    channel_map=MRI_LAX_4CH_SEGMENTED_CHANNEL_MAP,
)
lax_4ch_segmented_224_16_3_w = TensorMap(
    'lax_4ch_segmented_224_16_3', Interpretation.CATEGORICAL, shape=(160, 224, 16, len(MRI_LAX_4CH_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slices('cine_segmented_lax_4ch_annotated_', step=3),
    channel_map=MRI_LAX_4CH_SEGMENTED_CHANNEL_MAP,
    loss=weighted_crossentropy([0.01, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 1.0, 1.0, 1.0, 1.0, 5.0, 0.5, 10.0, 10.0]),
)
sax_segmented_b6 = TensorMap(
    'sax_segmented_b6',
    Interpretation.CATEGORICAL,
    shape=(256, 256, 50, 11),
    tensor_from_file=_segmented_dicom_slices(
        'cine_segmented_sax_b6_annotated_',
    ),
    channel_map=MRI_SAX_SEGMENTED_CHANNEL_MAP,
)
sax_segmented_b6_192 = TensorMap(
    'sax_segmented_b6',
    Interpretation.CATEGORICAL,
    shape=(192, 192, 50, 11),
    tensor_from_file=_segmented_dicom_slices(
        'cine_segmented_sax_b6_annotated_',
    ),
    channel_map=MRI_SAX_SEGMENTED_CHANNEL_MAP,
)

segmented_aorta_diastole = TensorMap(
    'segmented_aorta_diastole', Interpretation.CATEGORICAL, shape=(192, 256, len(MRI_AO_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slices('cine_segmented_ao_dist_annotated_'), channel_map=MRI_AO_SEGMENTED_CHANNEL_MAP,
)
cine_segmented_ao_dist = TensorMap(
    'cine_segmented_ao_dist', Interpretation.CATEGORICAL, shape=(160, 192, 100, len(MRI_AO_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slices('cine_segmented_ao_dist_annotated_'), channel_map=MRI_AO_SEGMENTED_CHANNEL_MAP,
)
cine_segmented_lvot = TensorMap(
    'cine_segmented_lvot', Interpretation.CATEGORICAL, shape=(208, 160, 50, len(MRI_LVOT_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slices('cine_segmented_lvot_annotated_'), channel_map=MRI_LVOT_SEGMENTED_CHANNEL_MAP,
)
cine_segmented_lvot_208_192_16_3 = TensorMap(
    'cine_segmented_lvot_208_192_16_3', Interpretation.CATEGORICAL, shape=(208, 192, 16, len(MRI_LVOT_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slices('cine_segmented_lvot_annotated_', step=3), channel_map=MRI_LVOT_SEGMENTED_CHANNEL_MAP,
)
flow_segmented = TensorMap(
    'flow_segmented', Interpretation.CATEGORICAL, shape=(192, 192, 30, len(MRI_AO_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slices('flow_250_tp_aov_bh_epat_annotated_'), channel_map=MRI_AO_SEGMENTED_CHANNEL_MAP,
)

liver_shmolli_segmented = TensorMap(
    'liver_shmolli_segmented',
    Interpretation.CATEGORICAL,
    shape=(288, 384, len(MRI_LIVER_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slices(
        'liver_shmolli_segmented_annotated_', path_prefix='ukb_liver_mri',
    ),
    channel_map=MRI_LIVER_SEGMENTED_CHANNEL_MAP,
)


def sax_tensor(b_series_prefix):
    def sax_tensor_from_file(tm, hd5, dependents={}):
        missing = 0
        tensor = np.zeros(tm.shape, dtype=np.float32)
        if tm.dependent_map is not None:
            dependents[tm.dependent_map] = np.zeros(tm.dependent_map.shape, dtype=np.float32)
        if tm.axes() == 3:
            for b in range(tm.shape[-1]):
                try:
                    tm_shape = (tm.shape[0], tm.shape[1])
                    tensor[:, :, b] = pad_or_crop_array_to_shape(tm_shape, np.array(hd5[f'{tm.path_prefix}/{b_series_prefix}/{(50*b)+1}/instance_0'], dtype=np.float32))
                except KeyError:
                    missing += 1
                    tensor[:, :, b] = 0
        else:
            for b in range(tm.shape[-2]):
                try:
                    tm_shape = (tm.shape[0], tm.shape[1])
                    tensor[:, :, b, 0] = pad_or_crop_array_to_shape(tm_shape, np.array(hd5[f'{tm.path_prefix}/{b_series_prefix}/{(50*b)+1}/instance_0'], dtype=np.float32))
                    if tm.dependent_map is not None:
                        index_tensor = pad_or_crop_array_to_shape(tm_shape, np.array(hd5[f'{tm.path_prefix}/{b_series_prefix}_mask_b{b}/instance_0'], dtype=np.float32))
                        dependents[tm.dependent_map][:, :, b, :] = to_categorical(index_tensor, tm.dependent_map.shape[-1])
                except KeyError:
                    missing += 1
                    tensor[:, :, b, 0] = 0
                    if tm.dependent_map is not None:
                        dependents[tm.dependent_map][:, :, b, MRI_SEGMENTED_CHANNEL_MAP['background']] = 1
            if missing == tm.shape[-2]:
                raise ValueError(f'Could not find any slices in {tm.name} was hoping for {tm.shape[-2]} looked at: {tm.path_prefix}/{b_series_prefix}')
        return tensor
    return sax_tensor_from_file


sax_all_diastole_segmented = TensorMap(
    'sax_all_diastole_segmented', Interpretation.CATEGORICAL, shape=(256, 256, 13, 3),
    channel_map=MRI_SEGMENTED_CHANNEL_MAP,
)
sax_all_diastole_segmented_weighted = TensorMap(
    'sax_all_diastole_segmented', Interpretation.CATEGORICAL, shape=(256, 256, 13, 3),
    channel_map=MRI_SEGMENTED_CHANNEL_MAP,
    loss=weighted_crossentropy(
        [1.0, 40.0, 40.0], 'sax_all_diastole_segmented',
    ),
)
sax_all_diastole_192_segmented_weighted = TensorMap(
    'sax_all_diastole_segmented', Interpretation.CATEGORICAL, shape=(192, 192, 13, 3),
    channel_map=MRI_SEGMENTED_CHANNEL_MAP,
    loss=weighted_crossentropy(
        [1.0, 40.0, 40.0], 'sax_all_diastole_segmented',
    ),
)

sax_all_diastole = TensorMap(
    'sax_all_diastole', shape=(256, 256, 13, 1), tensor_from_file=sax_tensor('diastole'),
    path_prefix='ukb_cardiac_mri',
)
sax_all_diastole_weighted = TensorMap(
    'sax_all_diastole', shape=(256, 256, 13, 1), tensor_from_file=sax_tensor('diastole'),
    dependent_map=sax_all_diastole_segmented_weighted, path_prefix='ukb_cardiac_mri',
)

sax_all_diastole_192 = TensorMap(
    'sax_all_diastole', shape=(192, 192, 13, 1), tensor_from_file=sax_tensor('diastole'),
    dependent_map=sax_all_diastole_segmented, path_prefix='ukb_cardiac_mri',
)
sax_all_diastole_192_weighted = TensorMap(
    'sax_all_diastole', shape=(192, 192, 13, 1), tensor_from_file=sax_tensor('diastole'),
    dependent_map=sax_all_diastole_segmented_weighted, path_prefix='ukb_cardiac_mri',
)

sax_all_systole_segmented = TensorMap(
    'sax_all_systole_segmented', Interpretation.CATEGORICAL, shape=(256, 256, 13, 3),
    channel_map=MRI_SEGMENTED_CHANNEL_MAP,
)
sax_all_systole_segmented_weighted = TensorMap(
    'sax_all_systole_segmented_weighted', Interpretation.CATEGORICAL, shape=(256, 256, 13, 3),
    channel_map=MRI_SEGMENTED_CHANNEL_MAP,
    loss=weighted_crossentropy([1.0, 40.0, 40.0], 'sax_all_systole_segmented'),
)


sax_all_systole = TensorMap(
    'sax_all_systole', shape=(256, 256, 13, 1), tensor_from_file=sax_tensor('systole'),
    dependent_map=sax_all_systole_segmented,
)
sax_all_systole_weighted = TensorMap(
    'sax_all_systole_weighted', shape=(256, 256, 13, 1), tensor_from_file=sax_tensor('systole'),
    dependent_map=sax_all_systole_segmented_weighted,
)


def all_sax_tensor(total_b_slices=13):
    def sax_tensor_from_file(tm, hd5, dependents={}):
        missing = 0
        tensor = np.zeros(tm.shape, dtype=np.float32)
        dependents[tm.dependent_map] = np.zeros(tm.dependent_map.shape, dtype=np.float32)
        for b in range(total_b_slices):
            try:
                tm_shape = (tm.shape[0], tm.shape[1])
                tensor[:, :, b, 0] = pad_or_crop_array_to_shape(tm_shape, np.array(hd5[f'{tm.hd5_key_guess()}diastole_frame_b{b}/instance_0'], dtype=np.float32))
                index_tensor = pad_or_crop_array_to_shape(tm_shape, np.array(hd5[f'{tm.hd5_key_guess()}diastole_mask_b{b}/instance_0'], dtype=np.float32))
                dependents[tm.dependent_map][:, :, b, :] = to_categorical(index_tensor, tm.dependent_map.shape[-1])
                tensor[:, :, b + total_b_slices, 0] = pad_or_crop_array_to_shape(tm_shape, np.array(hd5[f'{tm.hd5_key_guess()}systole_frame_b{b}/instance_0'], dtype=np.float32))
                index_tensor = pad_or_crop_array_to_shape(tm_shape, np.array(hd5[f'{tm.hd5_key_guess()}systole_mask_b{b}/instance_0'], dtype=np.float32))
                dependents[tm.dependent_map][:, :, b + total_b_slices, :] = to_categorical(index_tensor, tm.dependent_map.shape[-1])
            except KeyError as e:
                missing += 1
                tensor[:, :, b, 0] = 0
                dependents[tm.dependent_map][:, :, b, MRI_SEGMENTED_CHANNEL_MAP['background']] = 1
        if missing == tm.shape[-2]:
            raise ValueError(f'Could not find any slices in {tm.name} was hoping for {tm.shape[-2]}')
        return tensor
    return sax_tensor_from_file


sax_all_segmented = TensorMap(
    'sax_all_segmented', Interpretation.CATEGORICAL, shape=(
    256, 256, 26, 3,
    ), channel_map=MRI_SEGMENTED_CHANNEL_MAP,
)
sax_all_segmented_weighted = TensorMap(
    'sax_all_segmented_weighted', Interpretation.CATEGORICAL, shape=(256, 256, 26, 3),
    channel_map=MRI_SEGMENTED_CHANNEL_MAP, loss=weighted_crossentropy(
        [1.0, 40.0, 40.0], 'sax_all_segmented',
    ),
)

sax_all = TensorMap(
    'sax_all', shape=(
    256, 256, 26, 1,
    ), tensor_from_file=all_sax_tensor(), dependent_map=sax_all_segmented,
)
sax_all_weighted = TensorMap(
    'sax_all_weighted', shape=(
    256, 256, 26, 1,
    ), tensor_from_file=all_sax_tensor(), dependent_map=sax_all_segmented_weighted,
)


def _slice_tensor_with_segmentation(tensor_key, segmentation_key, path_prefix='ukb_cardiac_mri', max_slices=100):
    def _slice_tensor_from_file(tm, hd5, dependents={}):
        for i in range(max_slices):
            if f'/{path_prefix}/{segmentation_key}{i + 1}' in hd5:
                if tm.shape[-1] == 1:
                    t = pad_or_crop_array_to_shape(tm.shape[:-1], np.array(hd5[f'{path_prefix}/{tensor_key}'][..., i], dtype=np.float32))
                    tensor = np.expand_dims(t, axis=-1)
                else:
                    tensor = pad_or_crop_array_to_shape(tm.shape, np.array(hd5[f'{path_prefix}/{tensor_key}'][..., i], dtype=np.float32))
                return tensor
        raise ValueError(f'No segmented slice found for {tm.name} prefix {segmentation_key}')
    return _slice_tensor_from_file


aorta_slice_jamesp = TensorMap(
    'aorta_slice_jamesp', shape=(200, 240, 1), normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_tensor_with_segmentation('cine_segmented_ao_dist/instance_0', 'cine_segmented_ao_dist_jamesp_annotated_'),
)
aorta_slice_nekoui = TensorMap(
    'aorta_slice_nekoui', shape=(200, 240, 1), normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_tensor_with_segmentation('cine_segmented_ao_dist/instance_0', 'cine_segmented_ao_dist_nekoui_annotated_'),
)
lvot_slice_jamesp = TensorMap(
    'lvot_slice_jamesp', shape=(200, 240, 1), normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_tensor_with_segmentation('cine_segmented_lvot/instance_0', 'cine_segmented_lvot_jamesp_annotated_'),
)
lvot_slice_nekoui = TensorMap(
    'lvot_slice_nekoui', shape=(200, 240, 1), normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_tensor_with_segmentation('cine_segmented_lvot/instance_0', 'cine_segmented_lvot_nekoui_annotated_'),
)
lax_2ch_slice_jamesp = TensorMap(
    'lax_2ch_slice_jamesp', shape=(192, 160, 1), normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_tensor_with_segmentation('cine_segmented_lax_2ch/instance_0', 'cine_segmented_lax_2ch_jamesp_annotated_'),
)
lax_3ch_slice_jamesp = TensorMap(
    'lax_3ch_slice_jamesp', shape=(192, 160, 1), normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_tensor_with_segmentation('cine_segmented_lax_3ch/instance_0', 'cine_segmented_lax_3ch_jamesp_annotated_'),
)
lax_4ch_slice_jamesp = TensorMap(
    'lax_4ch_slice_jamesp', shape=(160, 224, 1), normalization=ZeroMeanStd1(),
    tensor_from_file=_slice_tensor_with_segmentation('cine_segmented_lax_4ch/instance_0', 'cine_segmented_lax_4ch_jamesp_annotated_'),
)


def _segmented_dicom_slice(dicom_key_prefix, path_prefix='ukb_cardiac_mri', max_slices=100):
    def _segmented_dicom_tensor_from_file(tm, hd5, dependents={}):
        tensor = np.zeros(tm.shape, dtype=np.float32)
        for i in range(max_slices):
            slice_key = f'{dicom_key_prefix}{i + 1}'
            if f'{path_prefix}/{slice_key}' in hd5:
                categorical_index_slice = get_tensor_at_first_date(hd5, path_prefix, slice_key)
                categorical_one_hot = to_categorical(categorical_index_slice, len(tm.channel_map))
                tensor[..., :] = pad_or_crop_array_to_shape(tensor[..., :].shape, categorical_one_hot)
                return tensor
        raise ValueError(f'No segmented slice found for {tm.name} prefix {dicom_key_prefix}')
    return _segmented_dicom_tensor_from_file


cine_segmented_ao_dist_jamesp = TensorMap(
    'cine_segmented_ao_dist', Interpretation.CATEGORICAL, shape=(200, 240, len(MRI_AO_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slice('cine_segmented_ao_dist_jamesp_annotated_'), channel_map=MRI_AO_SEGMENTED_CHANNEL_MAP,
)
cine_segmented_ao_dist_nekoui = TensorMap(
    'cine_segmented_ao_dist', Interpretation.CATEGORICAL, shape=(200, 240, len(MRI_AO_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slice('cine_segmented_ao_dist_nekoui_annotated_'), channel_map=MRI_AO_SEGMENTED_CHANNEL_MAP,
)
cine_segmented_lvot_jamesp = TensorMap(
    'cine_segmented_lvot', Interpretation.CATEGORICAL, shape=(200, 240, len(MRI_LVOT_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slice('cine_segmented_lvot_jamesp_annotated_'), channel_map=MRI_LVOT_SEGMENTED_CHANNEL_MAP,
)
cine_segmented_lvot_nekoui = TensorMap(
    'cine_segmented_lvot', Interpretation.CATEGORICAL, shape=(200, 240, len(MRI_LVOT_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slice('cine_segmented_lvot_nekoui_annotated_'), channel_map=MRI_LVOT_SEGMENTED_CHANNEL_MAP,
)
cine_segmented_lax_2ch_jamesp = TensorMap(
    'cine_segmented_lax_2ch_slice', Interpretation.CATEGORICAL, shape=(192, 160, len(MRI_LAX_2CH_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slice('cine_segmented_lax_2ch_jamesp_annotated_'), channel_map=MRI_LAX_2CH_SEGMENTED_CHANNEL_MAP,
)
cine_segmented_lax_3ch_jamesp = TensorMap(
    'cine_segmented_lax_3ch_slice', Interpretation.CATEGORICAL, shape=(192, 160, len(MRI_LAX_3CH_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slice('cine_segmented_lax_3ch_jamesp_annotated_'), channel_map=MRI_LAX_3CH_SEGMENTED_CHANNEL_MAP,
)
cine_segmented_lax_4ch_jamesp = TensorMap(
    'cine_segmented_lax_4ch_slice', Interpretation.CATEGORICAL, shape=(160, 224, len(MRI_LAX_4CH_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slice('cine_segmented_lax_4ch_jamesp_annotated_'), channel_map=MRI_LAX_4CH_SEGMENTED_CHANNEL_MAP,
)
cine_segmented_lax_4ch_diastole = TensorMap(
    'cine_segmented_lax_4ch_diastole', Interpretation.CATEGORICAL, shape=(160, 224, len(MRI_LAX_4CH_SEGMENTED_CHANNEL_MAP)),
    tensor_from_file=_segmented_dicom_slice('cine_segmented_lax_4ch_annotated_'), channel_map=MRI_LAX_4CH_SEGMENTED_CHANNEL_MAP,
)


def _segmented_index_slices(key_prefix: str, shape: Tuple[int], path_prefix: str ='ukb_cardiac_mri') -> Callable:
    """Get semantic segmentation with label index as pixel values for an MRI slice"""
    def _segmented_dicom_tensor_from_file(tm, hd5, dependents={}):
        tensor = np.zeros(shape, dtype=np.float32)
        for i in range(shape[-1]):
            categorical_index_slice = get_tensor_at_first_date(
                hd5, path_prefix, key_prefix + str(i + 1),
            )
            tensor[..., i] = pad_or_crop_array_to_shape(
                shape[:-1], categorical_index_slice,
            )
        return tensor
    return _segmented_dicom_tensor_from_file


def _bounding_box_from_categorical(segmented_shape: Tuple[int], segmented_key: str, class_index: int) -> Callable:
    """Given an hd5 key of a semantic segmentation return a bounding box that covers the extent of a given class
    :param segmented_shape: The shape of each segmentation
    :param segmented_key: The key for the HD5 file where the segmentation is stored
    :param class_index: The index in the segmentation asssociated with the class we will find the bounding box for
    :return: A np.ndarray encoding a bounding box with min coordinates followed by max coordinates
            For example, a 2D bounding box will be returned as a 1D tensor of 4 numbers: [min_x, min_y, max_x, max_y]
            for a 3d bounding box we would get 6 numbers: [min_x, min_y, min_z max_x, max_y, max_z]
    """
    def bbox_from_file(tm, hd5, dependents={}):
        tensor = np.zeros(tm.shape, dtype=np.float32)
        index_tensor = pad_or_crop_array_to_shape(
            segmented_shape, np.array(hd5[segmented_key], dtype=np.float32),
        )
        bitmask = np.where(index_tensor == class_index)
        # Divide by 2 because we need min and max for each axis
        total_axes = tm.shape[-1] // 2
        for i in range(total_axes):
            tensor[i] = np.min(bitmask[i])
            tensor[i+total_axes] = np.max(bitmask[i])
        return tensor
    return bbox_from_file


def _bounding_box_from_callable(class_index: int, tensor_from_file_fxn: Callable) -> Callable:
    """Given a tensor_from_file function that returns a semantic segmentation find the bounding box that covers the extent of a given class
    :param class_index: The index in the segmentation asssociated with the class we will find the bounding box for
    :param tensor_from_file_fxn: A tensor from file function that returns a class index tensor.
            This tensor should NOT be one hot, ie the segmentation before `to_categorical` has been applied.
    :return: A np.ndarray encoding a bounding box with min coordinates followed by max coordinates
            For example, a 2D bounding box will be returned as a 1D tensor of 4 numbers: [min_x, min_y, max_x, max_y]
            for a 3d bounding box we would get 6 numbers: [min_x, min_y, min_z max_x, max_y, max_z]
    """
    def bbox_from_file(tm, hd5, dependents={}):
        tensor = np.zeros(tm.shape, dtype=np.float32)
        index_tensor = tensor_from_file_fxn(None, hd5)
        bitmask = np.where(index_tensor == class_index)
        # Divide by 2 because we need min and max for each axis
        total_axes = tm.shape[-1] // 2
        for i in range(total_axes):
            tensor[i] = np.min(bitmask[i])
            tensor[i+total_axes] = np.max(bitmask[i])
        return tensor
    return bbox_from_file


def _bounding_box_channel_map(total_axes: int) -> Dict[str, int]:
    channel_map = {}
    for i in range(total_axes):
        channel_map[f'min_axis_{i}'] = i
        channel_map[f'max_axis_{i}'] = i + total_axes
    return channel_map


def _make_index_tensor_from_file(index_map_name):
    def indexed_lvmass_tensor_from_file(tm, hd5, dependents={}):
        tensor = np.zeros(tm.shape, dtype=np.float32)
        for k in tm.channel_map:
            tensor = np.array(hd5[tm.path_prefix][k], dtype=np.float32)
        index = np.array(hd5[tm.path_prefix][index_map_name], dtype=np.float32)
        return tensor / index
    return indexed_lvmass_tensor_from_file


def _select_tensor_from_file(selection_predicate: Callable):
    def selected_tensor_from_file(tm, hd5, dependents={}):
        if not selection_predicate(hd5):
            raise ValueError(
                f'Tensor did not meet selection criteria:{selection_predicate.__name__} with Tensor Map:{tm.name}',
            )
        tensor = np.zeros(tm.shape, dtype=np.float32)
        for k in tm.channel_map:
            tensor = np.array(hd5[tm.path_prefix][k], dtype=np.float32)
        return tensor
    return selected_tensor_from_file


def _make_lvh_from_lvm_tensor_from_file(lvm_key, group_key='continuous', male_lvh_threshold=72, female_lvh_threshold=55):
    def lvh_from_lvm_tensor_from_file(tm, hd5, dependents={}):
        tensor = np.zeros(tm.shape, dtype=np.float32)
        lvm_indexed = float(hd5[group_key][lvm_key][0])
        index = 0
        if is_genetic_man(hd5) and lvm_indexed > male_lvh_threshold:
            index = 1
        elif is_genetic_woman(hd5) and lvm_indexed > female_lvh_threshold:
            index = 1
        tensor[index] = 1
        return tensor
    return lvh_from_lvm_tensor_from_file


def _make_fallback_tensor_from_file(tensor_keys):
    def fallback_tensor_from_file(tm, hd5, dependents={}):
        for k in tensor_keys:
            if k in hd5:
                return pad_or_crop_array_to_shape(tm.shape, np.array(hd5[k], dtype=np.float32))
        raise ValueError(f'No fallback tensor found from keys: {tensor_keys}')
    return fallback_tensor_from_file


lv_mass_dubois_index = TensorMap(
    'lv_mass_dubois_index', Interpretation.CONTINUOUS, activation='linear', loss='logcosh', loss_weight=1.0,
    tensor_from_file=_make_index_tensor_from_file('bsa_dubois'),
    channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8},
)
lv_mass_mosteller_index = TensorMap(
    'lv_mass_mosteller_index', Interpretation.CONTINUOUS, activation='linear', loss='logcosh', loss_weight=1.0,
    tensor_from_file=_make_index_tensor_from_file('bsa_mosteller'),
    channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8},
)
lv_mass_dubois_index_sentinel = TensorMap(
    'lv_mass_dubois_index', Interpretation.CONTINUOUS, activation='linear', sentinel=0, loss_weight=1.0,
    tensor_from_file=_make_index_tensor_from_file('bsa_dubois'),
    channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8},
)
lv_mass_mosteller_index_sentinel = TensorMap(
    'lv_mass_mosteller_index', Interpretation.CONTINUOUS, activation='linear', sentinel=0, loss_weight=1.0,
    tensor_from_file=_make_index_tensor_from_file('bsa_mosteller'),
    channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8},
)

lvm_dubois_index = TensorMap(
    'lvm_dubois_index', Interpretation.CONTINUOUS, activation='linear', loss='logcosh', loss_weight=1.0,
    tensor_from_file=_make_index_tensor_from_file('bsa_dubois'),
    channel_map={'LVM': 0}, normalization={'mean': 89.7, 'std': 24.8},
)
lvm_mosteller_index = TensorMap(
    'lvm_mosteller_index', Interpretation.CONTINUOUS, activation='linear', loss='logcosh', loss_weight=1.0,
    tensor_from_file=_make_index_tensor_from_file('bsa_mosteller'),
    channel_map={'LVM': 0}, normalization={'mean': 89.7, 'std': 24.8},
)
lvm_dubois_index_w4 = TensorMap(
    'lvm_dubois_index', Interpretation.CONTINUOUS, activation='linear', loss='logcosh', loss_weight=4.0,
    tensor_from_file=_make_index_tensor_from_file('bsa_dubois'),
    channel_map={'LVM': 0}, normalization={'mean': 89.7, 'std': 24.8},
)
lvm_mosteller_index_w4 = TensorMap(
    'lvm_mosteller_index', Interpretation.CONTINUOUS, activation='linear', loss='logcosh', loss_weight=4.0,
    tensor_from_file=_make_index_tensor_from_file('bsa_mosteller'),
    channel_map={'LVM': 0}, normalization={'mean': 89.7, 'std': 24.8},
)
lvm_dubois_index_sentinel = TensorMap(
    'lvm_dubois_index', Interpretation.CONTINUOUS, activation='linear', sentinel=0, loss_weight=1.0,
    tensor_from_file=_make_index_tensor_from_file('bsa_dubois'),
    channel_map={'LVM': 0}, normalization={'mean': 89.7, 'std': 24.8},
)
lvm_mosteller_index_sentinel = TensorMap(
    'lvm_mosteller_index', Interpretation.CONTINUOUS, activation='linear', sentinel=0, loss_weight=1.0,
    tensor_from_file=_make_index_tensor_from_file('bsa_mosteller'),
    channel_map={'LVM': 0}, normalization={'mean': 89.7, 'std': 24.8},
)


myocardial_mass_noheritable_men_only = TensorMap(
    'inferred_myocardial_mass_noheritable', Interpretation.CONTINUOUS, activation='linear', loss='logcosh',
    tensor_from_file=_select_tensor_from_file(is_genetic_man),
    channel_map={'inferred_myocardial_mass_noheritable': 0}, normalization={'mean': 100.0, 'std': 18.0},
)
myocardial_mass_noheritable_women_only = TensorMap(
    'inferred_myocardial_mass_noheritable', Interpretation.CONTINUOUS, activation='linear', loss='logcosh',
    tensor_from_file=_select_tensor_from_file(is_genetic_woman),
    channel_map={'inferred_myocardial_mass_noheritable': 0}, normalization={'mean': 78.0, 'std': 16.0},
)


lvh_from_indexed_lvm = TensorMap(
    'lvh_from_indexed_lvm', Interpretation.CATEGORICAL, channel_map={'no_lvh': 0, 'left_ventricular_hypertrophy': 1},
    tensor_from_file=_make_lvh_from_lvm_tensor_from_file(
        'adjusted_myocardium_mass_indexed',
    ),
)
lvh_from_indexed_lvm_weighted = TensorMap(
    'lvh_from_indexed_lvm', Interpretation.CATEGORICAL, channel_map={'no_lvh': 0, 'left_ventricular_hypertrophy': 1},
    tensor_from_file=_make_lvh_from_lvm_tensor_from_file(
        'adjusted_myocardium_mass_indexed',
    ),
    loss=weighted_crossentropy([1.0, 25.0], 'lvh_from_indexed_lvm'),
)
adjusted_myocardium_mass = TensorMap(
    'adjusted_myocardium_mass', Interpretation.CONTINUOUS, validator=make_range_validator(0, 400), path_prefix='continuous',
    loss='logcosh', channel_map={'adjusted_myocardium_mass': 0}, normalization={'mean': 89.70, 'std': 24.80},
)
adjusted_myocardium_mass_indexed = TensorMap(
    'adjusted_myocardium_mass_indexed', Interpretation.CONTINUOUS, validator=make_range_validator(0, 400),
    loss='logcosh', channel_map={'adjusted_myocardium_mass_indexed': 0}, path_prefix='continuous',
    normalization={'mean': 89.70, 'std': 24.80},
)
lvh_from_indexed_lvm_parented = TensorMap(
    'lvh_from_indexed_lvm', Interpretation.CATEGORICAL, channel_map={'no_lvh': 0, 'left_ventricular_hypertrophy': 1},
    tensor_from_file=_make_lvh_from_lvm_tensor_from_file(
        'adjusted_myocardium_mass_indexed',
    ),
    loss=weighted_crossentropy([1.0, 25.0], 'lvh_from_indexed_lvm_parented'),
    parents=[
        adjusted_myocardium_mass_indexed,
        adjusted_myocardium_mass,
    ],
)

shmolli_192i_both = TensorMap(
    'shmolli_192i', Interpretation.CONTINUOUS, shape=(288, 384, 7),
    tensor_from_file=_make_fallback_tensor_from_file(
        ['shmolli_192i', 'shmolli_192i_liver'],
    ),
)
shmolli_192i_both_4d = TensorMap(
    'shmolli_192i', Interpretation.CONTINUOUS, shape=(288, 384, 7, 1),
    tensor_from_file=_make_fallback_tensor_from_file(
        ['shmolli_192i', 'shmolli_192i_liver'],
    ),
)
shmolli_192i_both_instance1 = TensorMap(
    'shmolli_192i_instance1', Interpretation.CONTINUOUS, shape=(288, 384, 1),
    tensor_from_file=_make_fallback_tensor_from_file(
        ['shmolli_192i', 'shmolli_192i_liver'],
    ),
)
shmolli_192i_liver_only = TensorMap(
    'shmolli_192i', Interpretation.CONTINUOUS, shape=(288, 384, 7),
    tensor_from_file=_make_fallback_tensor_from_file(['shmolli_192i_liver']),
)

lax_3ch_lv_cavity_bbox_slice0 = TensorMap(
    'lax_3ch_lv_cavity_bbox_slice0', Interpretation.MESH, shape=(4,),
    tensor_from_file=_bounding_box_from_categorical(
        (160, 160), 'ukb_cardiac_mri/cine_segmented_lax_3ch_annotated_1/instance_0', MRI_LAX_3CH_SEGMENTED_CHANNEL_MAP['LV_Cavity'],
    ),
    channel_map=_bounding_box_channel_map(2),
)
lax_3ch_left_atrium_bbox_slice0 = TensorMap(
    'lax_3ch_left_atrium_bbox_slice0', Interpretation.MESH, shape=(4,),
    tensor_from_file=_bounding_box_from_categorical(
        (160, 160), 'ukb_cardiac_mri/cine_segmented_lax_3ch_annotated_1/instance_0', MRI_LAX_3CH_SEGMENTED_CHANNEL_MAP['left_atrium'],
    ),
    channel_map=_bounding_box_channel_map(2),
)

aorta_descending_tff = _bounding_box_from_categorical(
    (192, 224), 'ukb_cardiac_mri/cine_segmented_ao_dist_annotated_1/instance_0', MRI_AO_SEGMENTED_CHANNEL_MAP['descending_aorta'],
)
cine_segmented_ao_descending_aorta_bbox_slice0 = TensorMap(
    'cine_segmented_ao_descending_aorta_bbox_slice0', Interpretation.MESH, shape=(4,),
    tensor_from_file=aorta_descending_tff, channel_map=_bounding_box_channel_map(
        2,
    ),
)
aorta_ascending_tff = _bounding_box_from_categorical(
    (192, 224), 'ukb_cardiac_mri/cine_segmented_ao_dist_annotated_1/instance_0', MRI_AO_SEGMENTED_CHANNEL_MAP['ascending_aorta'],
)
cine_segmented_ao_ascending_aorta_bbox_slice0 = TensorMap(
    'cine_segmented_ao_ascending_aorta_bbox_slice0', Interpretation.MESH, shape=(4,),
    tensor_from_file=aorta_ascending_tff, channel_map=_bounding_box_channel_map(
        2,
    ),
)

lax_3ch_lv_cavity_bbox = TensorMap(
    'lax_3ch_lv_cavity_bbox', Interpretation.MESH, shape=(6,), channel_map=_bounding_box_channel_map(3),
    tensor_from_file=_bounding_box_from_callable(
        5, _segmented_index_slices(
        'cine_segmented_lax_3ch_annotated_', (192, 160, 50),
        ),
    ),
)

bbfc = _bounding_box_from_callable(
    MRI_AO_SEGMENTED_CHANNEL_MAP['descending_aorta'], _segmented_index_slices(
    'cine_segmented_ao_dist_annotated_', (192, 224, 100),
    ),
)
cine_segmented_ao_descending_aorta_bbox = TensorMap(
    'cine_segmented_ao_descending_aorta_bbox', Interpretation.MESH, shape=(6,), tensor_from_file=bbfc,
    channel_map=_bounding_box_channel_map(3),
)

abbfc = _bounding_box_from_callable(
    MRI_AO_SEGMENTED_CHANNEL_MAP['ascending_aorta'], _segmented_index_slices(
    'cine_segmented_ao_dist_annotated_', (192, 224, 100),
    ),
)

cine_segmented_ao_ascending_aorta_bbox = TensorMap(
    'cine_segmented_ao_ascending_aorta_bbox', Interpretation.MESH, shape=(6,), tensor_from_file=abbfc,
    channel_map=_bounding_box_channel_map(3),
)


###
lv_mass = TensorMap(
    'lv_mass', Interpretation.CONTINUOUS, activation='linear', loss='logcosh', validator=make_range_validator(0, 500),
    channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8},
)

# lv_mass_no0 = TensorMap(
#     'lv_mass', Interpretation.CONTINUOUS, activation='linear', loss=ignore_zeros_logcosh,
#     channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8},
# )

lv_mass_sentinel = TensorMap(
    'lv_mass', Interpretation.CONTINUOUS, activation='linear', sentinel=0,
    channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8},
)
LVM_sentinel = TensorMap(
    'LVM',  Interpretation.CONTINUOUS, normalization={'mean': 89.70372484725051, 'std': 24.803669503436304}, sentinel=0,
    validator=make_range_validator(-1, 300), channel_map={'LVM': 0},
)
lv_mass_prediction = TensorMap(
    'lv_mass_sentinel_prediction', Interpretation.CONTINUOUS, activation='linear', loss='logcosh', loss_weight=10.0,
    validator=make_range_validator(0, 300), channel_map={'lv_mass_sentinel_prediction': 0},
    normalization={'mean': 89.7, 'std': 24.8},
)
lv_mass_dubois_index_prediction = TensorMap(
    'lv_mass_dubois_index_sentinel_prediction', Interpretation.CONTINUOUS, activation='linear', loss='logcosh',
    validator=make_range_validator(0, 300), loss_weight=10.0,
    channel_map={'lv_mass_dubois_index_sentinel_prediction': 0}, normalization={'mean': 89.7, 'std': 24.8},
)
lv_mass_mosteller_index_prediction = TensorMap(
    'lv_mass_mosteller_index_sentinel_prediction', Interpretation.CONTINUOUS, activation='linear', loss='logcosh',
    validator=make_range_validator(0, 300), loss_weight=10.0,
    channel_map={'lv_mass_mosteller_index_sentinel_prediction': 0},
    normalization={'mean': 89.7, 'std': 24.8},
)

LVM_prediction = TensorMap(
    'LVM_sentinel_prediction',  Interpretation.CONTINUOUS, normalization={'mean': 89.70372484725051, 'std': 24.803669503436304},
    validator=make_range_validator(0, 300), channel_map={'LVM_sentinel_prediction': 0},
)

lvm_dubois_index_prediction = TensorMap(
    'lvm_dubois_index_sentinel_prediction', Interpretation.CONTINUOUS, activation='linear', loss='logcosh',
    validator=make_range_validator(0, 300), channel_map={'lvm_dubois_index_sentinel_prediction': 0},
    normalization={'mean': 42.0, 'std': 8.0},
)
lvm_mosteller_index_prediction = TensorMap(
    'lvm_mosteller_index_sentinel_prediction', Interpretation.CONTINUOUS, activation='linear', loss='logcosh',
    validator=make_range_validator(0, 300), channel_map={'lvm_mosteller_index_sentinel_prediction': 0},
    normalization={'mean': 42.0, 'std': 8.0},
)

LVM_prediction_sentinel = TensorMap(
    'LVM_sentinel_prediction',  Interpretation.CONTINUOUS, sentinel=0, channel_map={'LVM_sentinel_prediction': 0},
    normalization={'mean': 89.70372484725051, 'std': 24.803669503436304},
)
lvm_dubois_index_prediction_sentinel = TensorMap(
    'lvm_dubois_index_sentinel_prediction', Interpretation.CONTINUOUS, activation='linear', loss='logcosh',
    sentinel=0, channel_map={'lvm_dubois_index_sentinel_prediction': 0},
    normalization={'mean': 89.7, 'std': 24.8},
)
lvm_mosteller_index_prediction_sentinel = TensorMap(
    'lvm_mosteller_index_sentinel_prediction', Interpretation.CONTINUOUS, activation='linear', loss='logcosh',
    sentinel=0, channel_map={'lvm_mosteller_index_sentinel_prediction': 0},
    normalization={'mean': 89.7, 'std': 24.8},
)



end_systole_volume = TensorMap(
    'end_systole_volume', Interpretation.CONTINUOUS, activation='linear', validator=make_range_validator(0, 300),
    loss='logcosh', channel_map={'end_systole_volume': 0},
    normalization={'mean': 47.0, 'std': 10.0},
)
end_diastole_volume = TensorMap(
    'end_diastole_volume', Interpretation.CONTINUOUS, activation='linear', validator=make_range_validator(0, 400),
    loss='logcosh', channel_map={'end_diastole_volume': 0},
    normalization={'mean': 142.0, 'std': 21.0},
)
ejection_fraction = TensorMap(
    'ejection_fraction', Interpretation.CONTINUOUS, activation='linear', validator=make_range_validator(0.2, 0.9),
    normalization={'mean': 0.50, 'std': 0.046},
    loss='logcosh', loss_weight=1.0, channel_map={'ejection_fraction': 0},
)


# Apply correction from Sanghvi et al.Journal of Cardiovascular Magnetic Resonance 2016
corrected_extracted_lvedv = TensorMap(
    'corrected_extracted_lvedv', Interpretation.CONTINUOUS, activation='linear', validator=make_range_validator(0, 400),
    loss='logcosh', channel_map={'corrected_extracted_lvedv': 0},
    normalization={'mean': 142.0, 'std': 21.0},
)
corrected_extracted_lvef = TensorMap(
    'corrected_extracted_lvef', Interpretation.CONTINUOUS, activation='linear', validator=make_range_validator(0.2, 0.9),
    normalization={'mean': 0.50, 'std': 0.046},
    loss='logcosh', channel_map={'corrected_extracted_lvef': 0},
)
corrected_extracted_lvesv = TensorMap(
    'corrected_extracted_lvesv', Interpretation.CONTINUOUS, activation='linear', validator=make_range_validator(0, 300),
    loss='logcosh', channel_map={'corrected_extracted_lvesv': 0},
    normalization={'mean': 47.0, 'std': 10.0},
)

corrected_extracted_lvesv_sentinel = TensorMap(
    'corrected_extracted_lvesv', Interpretation.CONTINUOUS, activation='linear', sentinel=0.0,
    channel_map={'corrected_extracted_lvesv': 0}, normalization={'mean': 47.0, 'std': 10.0},
)
corrected_extracted_lvedv_sentinel = TensorMap(
    'corrected_extracted_lvedv', Interpretation.CONTINUOUS, activation='linear', sentinel=0.0,
    channel_map={'corrected_extracted_lvedv': 0}, normalization={'mean': 142.0, 'std': 21.0},
)
corrected_extracted_lvef_sentinel = TensorMap(
    'corrected_extracted_lvef', Interpretation.CONTINUOUS, activation='linear', sentinel=0.0,
    normalization={'mean': 0.50, 'std': 0.046}, channel_map={'corrected_extracted_lvef': 0},
)
corrected_extracted_lvef_sentinel = TensorMap(
    'corrected_extracted_lvef', Interpretation.CONTINUOUS, activation='linear', sentinel=0.0,
    normalization={'mean': 0.50, 'std': 0.046}, channel_map={'corrected_extracted_lvef': 0},
)

LA_2Ch_vol_max = TensorMap(
    'LA_2Ch_vol_max',  Interpretation.CONTINUOUS, normalization={'mean': 63.45582391534391, 'std': 22.548034481265972},
    validator=make_range_validator(0, 400), loss='logcosh', channel_map={'LA_2Ch_vol_max': 0}, path_prefix='continuous',
)
LA_2Ch_vol_min = TensorMap(
    'LA_2Ch_vol_min',  Interpretation.CONTINUOUS, normalization={'mean': 28.308681904761904, 'std': 15.842444310837582},
    validator=make_range_validator(0, 200), loss='logcosh', channel_map={'LA_2Ch_vol_min': 0}, path_prefix='continuous',
)
LA_4Ch_vol_max = TensorMap(
    'LA_4Ch_vol_max',  Interpretation.CONTINUOUS, normalization={'mean': 74.53903305263158, 'std': 25.448756860639776},
    validator=make_range_validator(0, 400), loss='logcosh', channel_map={'LA_4Ch_vol_max': 0}, path_prefix='continuous',
)
LA_4Ch_vol_min = TensorMap(
    'LA_4Ch_vol_min',  Interpretation.CONTINUOUS, normalization={'mean': 31.014961894736846, 'std': 17.146722819760804},
    validator=make_range_validator(0, 200), loss='logcosh', channel_map={'LA_4Ch_vol_min': 0}, path_prefix='continuous',
)
LA_Biplan_vol_max = TensorMap(
    'LA_Biplan_vol_max',  Interpretation.CONTINUOUS, normalization={'mean': 67.86355108225109, 'std': 21.793845470012105},
    validator=make_range_validator(0, 400), loss='logcosh', channel_map={'LA_Biplan_vol_max': 0}, path_prefix='continuous',
)
LA_Biplan_vol_min = TensorMap(
    'LA_Biplan_vol_min',  Interpretation.CONTINUOUS, normalization={'mean': 28.79685670995671, 'std': 15.43219634139272},
    validator=make_range_validator(0, 300), loss='logcosh', channel_map={'LA_Biplan_vol_min': 0}, path_prefix='continuous',
)
LVEDV = TensorMap(
    'LVEDV',  Interpretation.CONTINUOUS, normalization={'mean': 144.1479505192425, 'std': 34.39409859908663}, loss='logcosh',
    validator=make_range_validator(0, 500), channel_map={'LVEDV': 0}, path_prefix='continuous',
)
LVEF = TensorMap(
    'LVEF',  Interpretation.CONTINUOUS, normalization={'mean': 47.0, 'std': 10.0}, loss='logcosh', path_prefix='continuous',
    validator=make_range_validator(0, 500), channel_map={'LVEF': 0},
)
LVESV = TensorMap(
    'LVESV',  Interpretation.CONTINUOUS, normalization={'mean': 59.58324862553452, 'std': 21.186976544044025}, loss='logcosh',
    validator=make_range_validator(0, 400), channel_map={'LVESV': 0}, path_prefix='continuous',
)
LVM = TensorMap(
    'LVM',  Interpretation.CONTINUOUS, normalization={'mean': 89.70372484725051, 'std': 24.803669503436304}, loss='logcosh',
    validator=make_range_validator(0, 400), channel_map={'LVM': 0}, path_prefix='continuous',
)
LVSV = TensorMap(
    'LVSV',  Interpretation.CONTINUOUS, normalization={'mean': 84.85198120147119, 'std': 19.2700091046526}, loss='logcosh',
    validator=make_range_validator(0, 400), channel_map={'LVSV': 0}, path_prefix='continuous',
)
RA_4Ch_vol_max = TensorMap(
    'RA_4Ch_vol_max',  Interpretation.CONTINUOUS, normalization={'mean': 79.22289586811351, 'std': 26.504015552539048},
    validator=make_range_validator(0, 500), loss='logcosh', channel_map={'RA_4Ch_vol_max': 0}, path_prefix='continuous',
)
RA_4Ch_vol_min = TensorMap(
    'RA_4Ch_vol_min',  Interpretation.CONTINUOUS, normalization={'mean': 46.25831176961603, 'std': 20.002160080524803},
    validator=make_range_validator(0, 400), loss='logcosh', channel_map={'RA_4Ch_vol_min': 0}, path_prefix='continuous',
)
RVEDV = TensorMap(
    'RVEDV',  Interpretation.CONTINUOUS, normalization={'mean': 152.41239853151131, 'std': 37.15198900632509}, loss='logcosh',
    validator=make_range_validator(0, 500), channel_map={'RVEDV': 0}, path_prefix='continuous',
)
RVEF = TensorMap(
    'RVEF',  Interpretation.CONTINUOUS, normalization={'mean': 56.404863078182565, 'std': 6.526231365539632}, loss='logcosh',
    validator=make_range_validator(10, 200), channel_map={'RVEF': 0}, path_prefix='continuous',
)
RVESV = TensorMap(
    'RVESV',  Interpretation.CONTINUOUS, normalization={'mean': 67.61379869467673, 'std': 22.853189258914284}, loss='logcosh',
    validator=make_range_validator(0, 300), channel_map={'RVESV': 0}, path_prefix='continuous',
)
RVSV = TensorMap(
    'RVSV',  Interpretation.CONTINUOUS, normalization={'mean': 85.0908258288989, 'std': 19.30893645374548}, loss='logcosh',
    validator=make_range_validator(0, 200), channel_map={'RVSV': 0}, path_prefix='continuous',
)
LAQC = TensorMap(
    'LAQC',  Interpretation.CONTINUOUS, normalization={'mean': 1.2657977883096367, 'std': 0.5561369836438385}, loss='logcosh',
    validator=make_range_validator(0, 200), channel_map={'LAQC': 0}, path_prefix='continuous',
)
LVQC = TensorMap(
    'LVQC',  Interpretation.CONTINUOUS, normalization={'mean': 1.1737756714060033, 'std': 0.4620420984104567}, loss='logcosh',
    validator=make_range_validator(0, 200), channel_map={'LVQC': 0}, path_prefix='continuous',
)
RAQC = TensorMap(
    'RAQC',  Interpretation.CONTINUOUS, normalization={'mean': 1.1860189573459716, 'std': 0.4791815490882246}, loss='logcosh',
    validator=make_range_validator(0, 200), channel_map={'RAQC': 0}, path_prefix='continuous',
)
RVQC = TensorMap(
    'RVQC',  Interpretation.CONTINUOUS, normalization={'mean': 1.179699842022117, 'std': 0.4648958893626213}, loss='logcosh',
    validator=make_range_validator(0, 200), channel_map={'RVQC': 0}, path_prefix='continuous',
)

myocardial_mass = TensorMap(
    'myocardium_mass',  Interpretation.CONTINUOUS, validator=make_range_validator(0, 400), loss='logcosh', path_prefix='continuous',
    channel_map={'myocardium_mass': 0}, normalization={'mean': 89.70, 'std': 24.80},
)
myocardial_mass_noheritable = TensorMap(
    'inferred_myocardial_mass_noheritable',  Interpretation.CONTINUOUS, path_prefix='continuous',
    loss='logcosh', validator=make_range_validator(0, 400), normalization={'mean': 89.70, 'std': 24.80},
    channel_map={'inferred_myocardial_mass_noheritable': 0},
)
myocardial_mass_noheritable_sentinel = TensorMap(
    'inferred_myocardial_mass_noheritable',  Interpretation.CONTINUOUS, sentinel=0, loss='logcosh',
    normalization={'mean': 89.70, 'std': 24.80}, path_prefix='continuous',
    channel_map={'inferred_myocardial_mass_noheritable': 0},
)

myocardial_mass = TensorMap(
    'myocardium_mass',  Interpretation.CONTINUOUS, validator=make_range_validator(0, 400), loss='logcosh', path_prefix='continuous',
    channel_map={'myocardium_mass': 0}, normalization={'mean': 89.70, 'std': 24.80},
)


adjusted_myocardium_mass_sentinel = TensorMap(
    'adjusted_myocardium_mass', Interpretation.CONTINUOUS, validator=make_range_validator(0, 400), path_prefix='continuous',
    loss='logcosh', channel_map={'adjusted_myocardium_mass': 0}, normalization={'mean': 89.70, 'std': 24.80},
    sentinel=0.0,
)

adjusted_myocardium_mass_mse = TensorMap(
    'adjusted_myocardium_mass', Interpretation.CONTINUOUS, validator=make_range_validator(0, 400), path_prefix='continuous',
    loss='mse', channel_map={'adjusted_myocardium_mass': 0}, normalization={'mean': 89.70, 'std': 24.80},
)

# adjusted_myocardium_mass_y_true_mse = TensorMap(
#     'adjusted_myocardium_mass', Interpretation.CONTINUOUS, validator=make_range_validator(0, 400), path_prefix='continuous',
#     loss=y_true_times_mse, channel_map={'adjusted_myocardium_mass': 0}, normalization={'mean': 89.70, 'std': 24.80},
# )

# adjusted_myocardium_mass_y_true_sqr_mse = TensorMap(
#     'adjusted_myocardium_mass', Interpretation.CONTINUOUS, validator=make_range_validator(0, 400), path_prefix='continuous',
#     loss=y_true_squared_times_mse, channel_map={'adjusted_myocardium_mass': 0}, normalization={'mean': 89.70, 'std': 24.80},
# )

# adjusted_myocardium_mass_y_true_cube_mse = TensorMap(
#     'adjusted_myocardium_mass', Interpretation.CONTINUOUS, validator=make_range_validator(0, 400), path_prefix='continuous',
#     loss=y_true_cubed_times_mse, channel_map={'adjusted_myocardium_mass': 0}, normalization={'mean': 89.70, 'std': 24.80},
# )


# adjusted_myocardium_mass_y_true_sqr_logcosh = TensorMap(
#     'adjusted_myocardium_mass', Interpretation.CONTINUOUS, validator=make_range_validator(0, 400),
#     loss=y_true_squared_times_logcosh, channel_map={'adjusted_myocardium_mass': 0}, path_prefix='continuous',
#     normalization={'mean': 89.70, 'std': 24.80},
# )


proton_fat = TensorMap(
    '22402_Proton-density-fat-fraction-PDFF_2_0', Interpretation.CONTINUOUS, channel_map={'22402_Proton-density-fat-fraction-PDFF_2_0': 0},
    activation='linear', loss='logcosh',  annotation_units=1, path_prefix='continuous',
    validator=make_range_validator(0, 100), normalization={'mean': 3.91012, 'std': 4.64437},
)
liver_fat = TensorMap(
    '22402_Liver-fat-percentage_2_0', Interpretation.CONTINUOUS, channel_map={'22402_Liver-fat-percentage_2_0': 0},
    activation='linear', loss='logcosh',  annotation_units=1, path_prefix='continuous',
    validator=make_range_validator(0, 100), normalization={'mean': 3.91012, 'std': 4.64437},
)
liver_fat_sentinel = TensorMap(
    '22402_Liver-fat-percentage_2_0', Interpretation.CONTINUOUS, channel_map={'22402_Liver-fat-percentage_2_0': 0},
    normalization={'mean': 3.91012, 'std': 4.64437}, activation='linear', sentinel=0.0, path_prefix='continuous',
)
liver_fat_echo_predicted = TensorMap(
    'liver_fat_sentinel_prediction', Interpretation.CONTINUOUS, channel_map={'liver_fat_sentinel_prediction': 0},
    validator=make_range_validator(0, 100), normalization={'mean': 3.91012, 'std': 4.64437}, path_prefix='continuous', activation='linear', loss='logcosh',
)
liver_fat_echo_predicted_sentinel = TensorMap(
    'liver_fat_sentinel_prediction', Interpretation.CONTINUOUS, channel_map={'liver_fat_sentinel_prediction': 0},
    normalization={'mean': 3.91012, 'std': 4.64437}, activation='linear', path_prefix='continuous', sentinel=0.0,
)

gre_mullti_echo_10_te_liver = TensorMap('gre_mullti_echo_10_te_liver', shape=(160, 160, 10), loss='logcosh', normalization={'zero_mean_std1': 1.0})
gre_mullti_echo_10_te_liver_12bit = TensorMap('gre_mullti_echo_10_te_liver_12bit', shape=(160, 160, 10), loss='logcosh', normalization={'zero_mean_std1': 1.0})
lms_ideal_optimised_low_flip_6dyn = TensorMap('lms_ideal_optimised_low_flip_6dyn', shape=(232, 256, 36), loss='logcosh', normalization={'zero_mean_std1': 1.0})
lms_ideal_optimised_low_flip_6dyn_12bit = TensorMap('lms_ideal_optimised_low_flip_6dyn_12bit', shape=(232, 256, 36), loss='logcosh', normalization={'zero_mean_std1': 1.0})
lms_ideal_optimised_low_flip_6dyn_4slice = TensorMap('lms_ideal_optimised_low_flip_6dyn_4slice', shape=(232, 256, 4), loss='logcosh', normalization={'zero_mean_std1': 1.0})

shmolli_192i = TensorMap('shmolli_192i', shape=(288, 384, 7), normalization={'zero_mean_std1': 1.0})
shmolli_192i_liver = TensorMap('shmolli_192i_liver', shape=(288, 384, 7), normalization={'zero_mean_std1': 1.0})
shmolli_192i_12bit = TensorMap('shmolli_192i_12bit', shape=(288, 384, 7), normalization={'zero_mean_std1': 1.0})
shmolli_192i_fitparams = TensorMap('shmolli_192i_fitparams', shape=(288, 384, 7), normalization={'zero_mean_std1': 1.0})
shmolli_192i_t1map = TensorMap('shmolli_192i_t1map', shape=(288, 384, 2), normalization={'zero_mean_std1': 1.0})

sax_pixel_width = TensorMap(
    'mri_pixel_width_cine_segmented_sax_inlinevf', Interpretation.CONTINUOUS, annotation_units=2, channel_map={'mri_pixel_width_cine_segmented_sax_inlinevf': 0},
    validator=make_range_validator(0, 4), normalization={'mean': 1.83, 'std': 0.1},
)
sax_pixel_height = TensorMap(
    'mri_pixel_height_segmented_sax_inlinevf', Interpretation.CONTINUOUS, annotation_units=2, channel_map={'mri_pixel_height_cine_segmented_sax_inlinevf': 0},
    validator=make_range_validator(0, 4), normalization={'mean': 1.83, 'std': 0.1},
)

ejection_fractionp = TensorMap(
    'ejection_fraction', Interpretation.CONTINUOUS, activation='linear',
    normalization={'mean': 0.50, 'std': 0.046},
    loss='logcosh', loss_weight=1.0, channel_map={'ejection_fraction': 0},
    parents=[end_systole_volume, end_diastole_volume],
)

cine_segmented_sax_b1 = TensorMap('cine_segmented_sax_b1', shape=(256, 256, 50), loss='mse')
cine_segmented_sax_b2 = TensorMap('cine_segmented_sax_b2', shape=(256, 256, 50), loss='mse')
cine_segmented_sax_b4 = TensorMap('cine_segmented_sax_b4', shape=(256, 256, 50), loss='mse')
cine_segmented_sax_b6 = TensorMap('cine_segmented_sax_b6', shape=(256, 256, 50), loss='mse')

cine_segmented_lax_2ch = TensorMap('cine_segmented_lax_2ch', shape=(256, 256, 50), normalization=ZeroMeanStd1())
cine_segmented_lax_3ch = TensorMap('cine_segmented_lax_3ch', shape=(256, 256, 50), normalization=ZeroMeanStd1())
cine_segmented_lax_4ch = TensorMap('cine_segmented_lax_4ch', shape=(256, 256, 50), normalization=ZeroMeanStd1())

cine_segmented_lax_2ch_4d = TensorMap('cine_segmented_lax_2ch_4d', shape=(256, 256, 50, 1), normalization=ZeroMeanStd1())
cine_segmented_lax_3ch_4d = TensorMap('cine_segmented_lax_3ch_4d', shape=(256, 256, 50, 1), normalization=ZeroMeanStd1())
cine_segmented_lax_4ch_4d = TensorMap('cine_segmented_lax_4ch_4d', shape=(256, 256, 50, 1), normalization=ZeroMeanStd1())

lax_view_detect = TensorMap(
    'lax_view_detect', Interpretation.CATEGORICAL,
    channel_map={
        'cine_segmented_lax_2ch': 0, 'cine_segmented_lax_3ch': 1,
        'cine_segmented_lax_4ch': 2,
    },
)

sax_view_detect = TensorMap(
    'sax_view_detect', Interpretation.CATEGORICAL,
    channel_map={
        'cine_segmented_sax_b1': 0, 'cine_segmented_sax_b2': 1,
        'cine_segmented_sax_b3': 2, 'cine_segmented_sax_b4': 3,
        'cine_segmented_sax_b5': 4, 'cine_segmented_sax_b6': 5,
        'cine_segmented_sax_b7': 6, 'cine_segmented_sax_b8': 7,
        'cine_segmented_sax_b9': 8, 'cine_segmented_sax_b10': 9,
        'cine_segmented_sax_b11': 10,
    },
)

slax_view_detect = TensorMap(
    'slax_view_detect', Interpretation.CATEGORICAL,
    channel_map={
        'cine_segmented_lax_2ch': 11, 'cine_segmented_lax_3ch': 12,
        'cine_segmented_lax_4ch': 13, 'cine_segmented_sax_b1': 0,
        'cine_segmented_sax_b2': 1, 'cine_segmented_sax_b3': 2,
        'cine_segmented_sax_b4': 3, 'cine_segmented_sax_b5': 4,
        'cine_segmented_sax_b6': 5, 'cine_segmented_sax_b7': 6,
        'cine_segmented_sax_b8': 7, 'cine_segmented_sax_b9': 8,
        'cine_segmented_sax_b10': 9, 'cine_segmented_sax_b11': 10,
    },
)



def mri_adiposity_rotate_image(image, angle):
    image_center = tuple(np.array(image.shape[1::-1]) / 2)
    rot_mat = cv2.getRotationMatrix2D(image_center, angle, 1.0)
    result  = cv2.warpAffine(image, rot_mat, image.shape[1::-1], flags=cv2.INTER_LINEAR)
    return result


def mri_adiposity_translate_image(image, shape, steps):
    M = np.float32([[1, 0, steps], [0, 1, steps]])
    image = cv2.warpAffine(image, M, shape)
    return image


def mri_adiposity_uncompress_data(compressed_data: h5py.Dataset) -> np.ndarray:
    return np.frombuffer(
        blosc.decompress(compressed_data[()]), dtype=np.uint16
    ).reshape(compressed_data.attrs["shape"]).astype(np.float32)


def mdrk_adiposity_mri_2dprojection_view_mixup(data: np.ndarray, alpha: float = 1.0, midpoint: int = 144) -> np.ndarray:
    """Mixup on the CPU: this special version computes two separate mixup
    operations on a single input tensor where the sagittal and coronal
    2D projections are placed side-by-side. Failure to compute mixup this
    way for side-by-side images does not give the expected results.

    Args:
        data (np.ndarray): Input tensor
        alpha (float, optional): Mixing factor: passed as a parameter to the Beta distribution. Defaults to 1.0.
        midpoint (int, optional): Number of pixels for the image to the left. Defaults to 144.

    Returns:
        np.ndarray: Side-by-side mixup of input tensor
    """
    batch_size = len(data)
    weights1 = np.random.beta(alpha, alpha, batch_size)
    weights2 = np.random.beta(alpha, alpha, batch_size)
    index1   = np.random.permutation(batch_size)
    index2   = np.random.permutation(batch_size)
    x2 = np.zeros(data.shape)
    # Midpoint is scaled and should be 224/1.5546875 = 144
    for i in range(len(weights1)):
        x2[i][:, 0:midpoint, :] = data[i][:, 0:midpoint, :] * weights1[i] + data[index1[i]][:, 0:midpoint, :] * (1 - weights1[i])
        x2[i][:, midpoint:, :]  = data[i][:, midpoint:, :]  * weights2[i] + data[index2[i]][:, midpoint:, :]  * (1 - weights2[i])
    return x2


def mdrk_projection_single(field: str, instance: int = 2, augment: bool = False):
    def _mdrk_projection_single(tm, hd5, dependents={}):
        try:
            compressed_data = hd5["instance"][str(instance)][field]
        except Exception as e:
            raise Exception(e)
        #
        try:
            tensor = mri_adiposity_uncompress_data(compressed_data)
        except Exception as e:
            raise Exception(e)
        #
        tensor = cv2.resize(tensor, (224, 368))
        if augment:
            if np.random.random() > 0.5:
                M = np.float32(
                    [
                        [1, 0, np.random.randint(-15, 15)],
                        [0, 1, np.random.randint(-15, 15)],
                    ]
                )
                tensor = cv2.warpAffine(tensor, M, (224, 368))
            if np.random.random() > 0.5:
                tensor = cv2.flip(tensor, 1)
        tensor = np.expand_dims(tensor, axis=-1)
        return tensor
    return _mdrk_projection_single


def mdrk_projection_single_both_views_all_stationwide_normalization(
    instance: int = 2,
    augment: bool = False,
    stationwise_normalization=True,
    normalize_histogram=True,
    clahe_amount=5,
    clahe_clip=2.0,
):
    """This function wrapper constructs a new image with the coronal and sagittal 2D
    projections side-by-side and each capture/reconstruciton modality stacked in the
    channels. Returns a (237, 256, 4) tensor. This function does *NOT* respect the
    desired shape provided in the TensorMap instance.

    Requirements:

    This subroutine requires that the target HDF5 file has the following datasets:
    * /instance/{instance}/w_sagittal and /instance/{instance}/w_coronal
    * /instance/{instance}/f_sagittal and /instance/{instance}/f_coronal
    * /instance/{instance}/in_sagittal and /instance/{instance}/in_coronal
    * /instance/{instance}/opp_sagittal and /instance/{instance}/opp_coronal

    that are compressed with blosc as 16-bit unsigned integers. Each dataset must also
    have the attribute `shape`.

    Example ML4H usage:

    ```python
    >>> actual_train_tm = ml4h.TensorMap(
        'mdrk_projection_single_both_views_all_stationwide_normalization',
        tensor_from_file=mdrk_projection_single_both_views_all_stationwide_normalization(instance = 2, augment=True),
        shape=(237, 256, 4),
        normalization=ml4h.ZeroMeanStd1(),
    )
    ```

    Args:
        instance (int, optional): UK Biobank instance numbering. Defaults to 2.
        augment (bool, optional): Augment data: includes a translation, rotation, and axis flip. Defaults to False.
        stationwise_normalization (bool, optional): Normalize each station separately before appending as a channel. Defaults to True.
        normalize_histogram (bool, optional): Normalize intensity histogram using CLAHE. Defaults to True.
        clahe_amount (int, optional): Size of CLAHE kernel. Defaults to 5.
        clahe_clip (float, optional): Clip limit for the CLAHE kernel. Defaults to 2.0.
    """
    def _mdrk_projection_single_both_views_all_stationwide_normalization(
        tm, hd5, dependents={}
    ):
        # 174 + 224 = 398 -> (368, 398)
        # map to (237, 256)
        clahe = cv2.createCLAHE(
            clipLimit=clahe_clip, 
            tileGridSize=(clahe_amount, clahe_amount)
        )

        do_augment = False
        do_flip = False
        rand_angle = 0.0
        rand_move = 0.0
        if augment:
            do_augment = True
            if np.random.random() > 0.5:
                do_flip = True
            if np.random.random() > 0.5:
                rand_angle = np.random.randint(-5, 5)
            if np.random.random() > 0.5:
                rand_move = np.random.randint(-16, 16)
        
        prefixes = ["w", "f", "in", "opp"]
        tensor = np.zeros((368, 174 + 224, len(prefixes)))
        
        for p, i in zip(prefixes, range(len(prefixes))):
            # Coronal view
            compressed_data = hd5["instance"][str(instance)][f"{p}_coronal"]
            tensor_coronal = mri_adiposity_uncompress_data(compressed_data)
            
            if stationwise_normalization:
                tensor_coronal = TopKNormalize(50).normalize(tensor_coronal)
            if normalize_histogram:
                tensor_coronal = (
                    clahe.apply((tensor_coronal * 255.0).astype(np.uint8)).astype(
                        np.float32
                    )
                    / 255.0
                )
            
            tensor_coronal = cv2.resize(tensor_coronal, (224, 368))
            if do_augment:
                if do_flip:
                    tensor_coronal = cv2.flip(tensor_coronal, 1)
                tensor_coronal = mri_adiposity_translate_image(tensor_coronal, (224, 368), rand_move)
                tensor_coronal = mri_adiposity_rotate_image(tensor_coronal, rand_angle)
            tensor[..., 0:224, i] = tensor_coronal
            
            compressed_data = hd5["instance"][str(instance)][f"{p}_sagittal"]
            tensor_sagittal = mri_adiposity_uncompress_data(compressed_data)
            
            if stationwise_normalization:
                tensor_sagittal = TopKNormalize(50).normalize(tensor_sagittal)
            if normalize_histogram:
                tensor_sagittal = (
                    clahe.apply((tensor_sagittal * 255.0).astype(np.uint8)).astype(
                        np.float32
                    )
                    / 255.0
                )
            #
            tensor_sagittal = cv2.resize(tensor_sagittal, (174, 368))
            if do_augment:
                tensor_sagittal = mri_adiposity_translate_image(
                    tensor_sagittal, (174, 368), rand_move
                )
                tensor_sagittal = mri_adiposity_rotate_image(tensor_sagittal, rand_angle)
            tensor[..., 224:, i] = tensor_sagittal
        tensor = cv2.resize(tensor, (256, 237))
        return tensor
    return _mdrk_projection_single_both_views_all_stationwide_normalization


def mdrk_projection_both_views_pretrained(
    instance: int = 2,
    augment: bool = False,
    stationwise_normalization=True,
    normalize_histogram=True,
    clahe_amount=5,
    clahe_clip=2.0,
):
    """This function wrapper constructs a new image with the coronal and sagittal 2D
    projections side-by-side for the water/fat reconstructions stacked in the
    channels followed by an empty channel to fit the expectations of pretrained image
    models. Returns a (237, 256, 3) tensor. This function does *NOT* respect the
    desired shape provided in the TensorMap instance.

    Requirements:

    This subroutine requires that the target HDF5 file has the following datasets:
    * /instance/{instance}/w_sagittal and /instance/{instance}/w_coronal
    * /instance/{instance}/f_sagittal and /instance/{instance}/f_coronal

    that are compressed with blosc as 16-bit unsigned integers. Each dataset must also
    have the attribute `shape`.

    Example ML4H usage:

    ```python
    >>> actual_train_tm = ml4h.TensorMap(
        'mdrk_projection_both_views_pretrained',
        tensor_from_file=mdrk_projection_both_views_pretrained(instance = 2, augment=True),
        shape=(237, 256, 2),
        normalization=ml4h.ZeroMeanStd1(),
    )
    ```

    Args:
        instance (int, optional): UK Biobank instance numbering. Defaults to 2.
        augment (bool, optional): Augment data: includes a translation, rotation, and axis flip. Defaults to False.
        stationwise_normalization (bool, optional): Normalize each station separately before appending as a channel. Defaults to True.
        normalize_histogram (bool, optional): Normalize intensity histogram using CLAHE. Defaults to True.
        clahe_amount (int, optional): Size of CLAHE kernel. Defaults to 5.
        clahe_clip (float, optional): Clip limit for the CLAHE kernel. Defaults to 2.0.
    """
    def _mdrk_projection_both_views_pretrained(tm, hd5, dependents={}):
        do_augment = False
        do_flip = False
        rand_angle = 0.0
        rand_move = 0.0
        cclip = clahe_clip
        camount = clahe_amount
        if augment:
            do_augment = True
            if np.random.random() > 0.5:
                do_flip = True
            rand_angle = np.random.randint(-5, 5)
            rand_move = np.random.randint(-16, 16)
            cclip = np.random.randint(0, 5)
            camount = np.random.randint(1, 10)
        
        clahe = cv2.createCLAHE(
            clipLimit=cclip, tileGridSize=(camount, camount)
        )
        prefixes = ["w", "f"]
        tensor = np.zeros((368, 174 + 224, 3), dtype=np.float32)
        
        for p, i in zip(prefixes, range(len(prefixes))):
            # Coronal view
            compressed_data = hd5["instance"][str(instance)][f"{p}_coronal"]
            tensor_coronal = mri_adiposity_uncompress_data(compressed_data)
            
            if stationwise_normalization:
                tensor_coronal = TopKNormalize(50).normalize(tensor_coronal)
            if normalize_histogram:
                tensor_coronal = (
                    clahe.apply((tensor_coronal * 255.0).astype(np.uint8)).astype(
                        np.float32
                    )
                    / 255.0
                )
            
            tensor_coronal = cv2.resize(tensor_coronal, (224, 368))
            if do_augment:
                if do_flip:
                    tensor_coronal = cv2.flip(tensor_coronal, 1)
                tensor_coronal = mri_adiposity_translate_image(tensor_coronal, (224, 368), rand_move)
                tensor_coronal = mri_adiposity_rotate_image(tensor_coronal, rand_angle)
            tensor[..., 0:224, i] = tensor_coronal
            
            compressed_data = hd5["instance"][str(instance)][f"{p}_sagittal"]
            # Sagittal view
            tensor_sagittal = mri_adiposity_uncompress_data(compressed_data)
            
            if stationwise_normalization:
                tensor_sagittal = TopKNormalize(50).normalize(tensor_sagittal)
            if normalize_histogram:
                tensor_sagittal = (
                    clahe.apply((tensor_sagittal * 255.0).astype(np.uint8)).astype(
                        np.float32
                    )
                    / 255.0
                )
            
            tensor_sagittal = cv2.resize(tensor_sagittal, (174, 368))
            if do_augment:
                tensor_sagittal = mri_adiposity_translate_image(
                    tensor_sagittal, (174, 368), rand_move
                )
                tensor_sagittal = mri_adiposity_rotate_image(tensor_sagittal, rand_angle)
            tensor[..., 224:, i] = tensor_sagittal
        tensor = cv2.resize(tensor, (256, 237))
        return tensor
    return _mdrk_projection_both_views_pretrained


mdrk_adiposity_mri_2dprojection_actual_train_tm = TensorMap(
    "mdrk_projection_single_both_views_all_stationwide_normalization",
    tensor_from_file=mdrk_projection_both_views_pretrained(
        instance=2, augment=True
    ),
    # shape=(368,174+224, 3), reshaped to (237, 256, 3)
    shape=(237, 256, 3),
    normalization=ZeroMeanStd1()
)

mdrk_adiposity_mri_2dprojection_actual_test_tm = TensorMap(
    "mdrk_projection_single_both_views_all_stationwide_normalization",
    tensor_from_file=mdrk_projection_both_views_pretrained(instance=2, augment=False),
    # shape=(368,174+224, 3), reshaped to (237, 256, 3)
    shape=(237, 256, 3),
    normalization=ZeroMeanStd1()
)

# Fake TMAP for compatibility with ML4H constructor.
mdrk_adiposity_mri_2dprojection_scalar_output_fake = TensorMap(
    "mdrk_adiposity_scalar_output_fake",
    shape=(1,),
    normalization=None,
    tensor_from_file=None
)
