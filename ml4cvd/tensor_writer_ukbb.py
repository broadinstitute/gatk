# tensor_writer_ukbb.py
#
# UK Biobank-specific tensor writing, SQL querying, data munging goes here
#

# Imports
import os
import re
import csv
import glob
import time
import shutil
import logging
import datetime
import operator
import tempfile
import traceback
from functools import partial
from itertools import product
from collections import Counter, defaultdict
from typing import Dict, List, Tuple, Optional, Any

import h5py
import imageio
import pydicom
import sqlite3
import zipfile
import matplotlib
matplotlib.use('Agg')
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw  # Polygon to mask
import xml.etree.ElementTree as et
from timeit import default_timer as timer
from scipy.ndimage.morphology import binary_closing, binary_erosion  # Morphological operator

from ml4cvd.plots import plot_value_counter, plot_histograms
from ml4cvd.defines import ECG_BIKE_LEADS, ECG_BIKE_MEDIAN_SIZE, ECG_BIKE_STRIP_SIZE, ECG_BIKE_FULL_SIZE, MRI_FRAMES
from ml4cvd.defines import MRI_TO_SEGMENT, MRI_SEGMENTED_CHANNEL_MAP, MRI_ANNOTATION_CHANNEL_MAP, MRI_ANNOTATION_NAME
from ml4cvd.defines import StorageType, IMAGE_EXT, TENSOR_EXT, DICOM_EXT, JOIN_CHAR, CONCAT_CHAR, HD5_GROUP_CHAR, DATE_FORMAT
from ml4cvd.defines import MRI_PIXEL_WIDTH, MRI_PIXEL_HEIGHT, MRI_SLICE_THICKNESS, MRI_PATIENT_ORIENTATION, MRI_PATIENT_POSITION



MRI_MIN_RADIUS = 2
MRI_MAX_MYOCARDIUM = 20
MRI_BIG_RADIUS_FACTOR = 0.9
MRI_SMALL_RADIUS_FACTOR = 0.19
MRI_MITRAL_VALVE_THICKNESS = 6
MRI_SWI_SLICES_TO_AXIS_SHIFT = 48
MRI_CARDIAC_SERIES = [
    'cine_segmented_lax_2ch', 'cine_segmented_lax_3ch', 'cine_segmented_lax_4ch', 'cine_segmented_sax_b1', 'cine_segmented_sax_b2',
    'cine_segmented_sax_b3', 'cine_segmented_sax_b4', 'cine_segmented_sax_b5', 'cine_segmented_sax_b6', 'cine_segmented_sax_b7',
    'cine_segmented_sax_b8', 'cine_segmented_sax_b9', 'cine_segmented_sax_b10', 'cine_segmented_sax_b11', 'cine_segmented_sax_b12',
    'cine_segmented_sax_b13', 'cine_segmented_sax_inlinevf', 'cine_segmented_lax_inlinevf', 'cine_segmented_ao_dist',
    'cine_segmented_lvot', 'flow_250_tp_aov_bh_epat@c_p', 'flow_250_tp_aov_bh_epat@c', 'flow_250_tp_aov_bh_epat@c_mag',
]
MRI_CARDIAC_SERIES_SEGMENTED = [series+'_segmented' for series in MRI_CARDIAC_SERIES]
MRI_BRAIN_SERIES = ['t1_p2_1mm_fov256_sag_ti_880', 't2_flair_sag_p2_1mm_fs_ellip_pf78']
MRI_NIFTI_FIELD_ID_TO_ROOT = {'20251': 'SWI', '20252': 'T1', '20253': 'T2_FLAIR'}
MRI_LIVER_SERIES = ['gre_mullti_echo_10_te_liver', 'lms_ideal_optimised_low_flip_6dyn', 'shmolli_192i', 'shmolli_192i_liver', 'shmolli_192i_fitparams', 'shmolli_192i_t1map']
MRI_LIVER_SERIES_12BIT = ['gre_mullti_echo_10_te_liver_12bit', 'lms_ideal_optimised_low_flip_6dyn_12bit', 'shmolli_192i_12bit', 'shmolli_192i_liver_12bit']
MRI_LIVER_IDEAL_PROTOCOL = ['lms_ideal_optimised_low_flip_6dyn', 'lms_ideal_optimised_low_flip_6dyn_12bit']

DICOM_MRI_FIELDS = ['20209', '20208', '20210', '20212', '20213', '20204', '20203', '20254', '20216', '20220', '20218', '20227', '20225', '20217']

ECG_BIKE_FIELD = '6025'
ECG_REST_FIELD = '20205'
ECG_TABLE_TAGS = ['RAmplitude', 'SAmplitude']
ECG_TAGS_TO_WRITE = [
    'VentricularRate', 'PQInterval', 'PDuration', 'QRSDuration', 'QTInterval', 'QTCInterval', 'RRInterval', 'PPInterval',
    'SokolovLVHIndex', 'PAxis', 'RAxis', 'TAxis', 'QTDispersion', 'QTDispersionBazett', 'QRSNum', 'POnset', 'POffset', 'QOnset',
    'QOffset', 'TOffset',
]
ECG_BIKE_SAMPLE_RATE = 500
ECG_BIKE_NUM_LEADS = 3
SECONDS_PER_MINUTE = 60


def write_tensors(
    a_id: str,
    xml_folder: str,
    zip_folder: str,
    output_folder: str,
    tensors: str,
    mri_unzip: str,
    mri_field_ids: List[int],
    xml_field_ids: List[int],
    zoom_x: int,
    zoom_y: int,
    zoom_width: int,
    zoom_height: int,
    write_pngs: bool,
    min_sample_id: int,
    max_sample_id: int,
    min_values_to_print: int,
) -> None:
    """Write tensors as HD5 files containing any kind of data from UK BioBank

    One HD5 file is generated per sample.  Each file may contain many tensor encodings of data including:
     survey responses, MRI, and ECG.

    :param a_id: User chosen string to identify this run
    :param xml_folder: Path to folder containing ECG XML files
    :param zip_folder: Path to folder containing zipped DICOM files
    :param output_folder: Folder to write outputs to (mostly for debugging)
    :param tensors: Folder to populate with HD5 tensors
    :param mri_unzip: Folder where zipped DICOM will be decompressed
    :param mri_field_ids: List of MRI field IDs from UKBB
    :param xml_field_ids: List of ECG field IDs from UKBB
    :param x: Maximum x dimension of MRIs
    :param y: Maximum y dimension of MRIs
    :param z: Maximum z dimension of MRIs
    :param zoom_x: x coordinate of the zoom
    :param zoom_y: y coordinate of the zoom
    :param zoom_width: width of the zoom
    :param zoom_height: height of the zoom
    :param write_pngs: write MRIs as PNG images for debugging
    :param min_sample_id: Minimum sample id to generate, for parallelization
    :param max_sample_id: Maximum sample id to generate, for parallelization
    :param min_values_to_print: Minimum number of samples that have responded to question for it to be included in the
            categorical or continuous dictionaries printed after tensor generation

    :return: None
    """
    stats = Counter()
    continuous_stats = defaultdict(list)
    sample_ids = range(min_sample_id, max_sample_id)
    for sample_id in sorted(sample_ids):

        start_time = timer()  # Keep track of elapsed execution time
        tp = os.path.join(tensors, str(sample_id) + TENSOR_EXT)
        if not os.path.exists(os.path.dirname(tp)):
            os.makedirs(os.path.dirname(tp))
        if _prune_sample(sample_id, min_sample_id, max_sample_id, mri_field_ids, xml_field_ids, zip_folder, xml_folder):
            continue
        try:
            with h5py.File(tp, 'w') as hd5:
                _write_tensors_from_zipped_dicoms(zoom_x, zoom_y, zoom_width, zoom_height, write_pngs, tensors, mri_unzip, mri_field_ids, zip_folder, hd5, sample_id, stats)
                _write_tensors_from_zipped_niftis(zip_folder, mri_field_ids, hd5, sample_id, stats)
                _write_tensors_from_xml(xml_field_ids, xml_folder, hd5, sample_id, write_pngs, stats, continuous_stats)
                stats['Tensors written'] += 1
        except AttributeError:
            logging.exception('Encountered AttributeError trying to write a UKBB tensor at path:{}'.format(tp))
            logging.info('Deleting attempted tensor at path:{}'.format(tp))
            os.remove(tp)
        except ValueError:
            logging.exception('Encountered ValueError trying to write a UKBB tensor at path:{}'.format(tp))
            logging.info('Deleting attempted tensor at path:{}'.format(tp))
            os.remove(tp)
        except RuntimeError:
            logging.exception('Encountered RuntimeError trying to write a UKBB tensor at path:{}'.format(tp))
            logging.info('Deleting attempted tensor at path:{}'.format(tp))
            os.remove(tp)
        except IndexError:
            logging.exception('Encountered IndexError trying to write a UKBB tensor at path:{}'.format(tp))
            logging.info('Deleting attempted tensor at path:{}'.format(tp))
            os.remove(tp)
        except OSError:
            logging.exception('Encountered OSError trying to write a UKBB tensor at path:{}'.format(tp))
            logging.info('Deleting attempted tensor at path:{}'.format(tp))
            os.remove(tp)

        end_time = timer()
        elapsed_time = end_time - start_time
        logging.info("Populated {} in {} seconds.".format(tp, elapsed_time))

    _dicts_and_plots_from_tensorization(a_id, output_folder, min_values_to_print, write_pngs, continuous_stats, stats)


def write_tensors_from_dicom_pngs(
    tensors, png_path, manifest_tsv, series, min_sample_id, max_sample_id, x=256, y=256,
    sample_header='sample_id', dicom_header='dicom_file',
    instance_header='instance_number', png_postfix='.png.mask.png',
    path_prefix='ukb_cardiac_mri',
):
    stats = Counter()
    reader = csv.reader(open(manifest_tsv), delimiter='\t')
    header = next(reader)
    logging.info(f"DICOM Manifest Header is:{header}")
    instance_index = header.index(instance_header)
    sample_index = header.index(sample_header)
    dicom_index = header.index(dicom_header)
    for row in reader:
        sample_id = row[sample_index]
        if not min_sample_id <= int(sample_id) < max_sample_id:
            continue
        stats[sample_header + '_' + sample_id] += 1
        dicom_file = row[dicom_index]
        try:
            png = imageio.imread(os.path.join(png_path, dicom_file + png_postfix))
            full_tensor = np.zeros((x, y), dtype=np.float32)
            full_tensor[:png.shape[0], :png.shape[1]] = png
            tensor_file = os.path.join(tensors, str(sample_id) + TENSOR_EXT)
            if not os.path.exists(os.path.dirname(tensor_file)):
                os.makedirs(os.path.dirname(tensor_file))
            with h5py.File(tensor_file, 'a') as hd5:
                tensor_name = series + '_annotated_' + row[instance_index]
                tp = tensor_path(path_prefix, tensor_name)
                if tp in hd5:
                    tensor = first_dataset_at_path(hd5, tp)
                    tensor[:] = full_tensor
                    stats['updated'] += 1
                else:
                    create_tensor_in_hd5(hd5, path_prefix, tensor_name, full_tensor, stats)
                    stats['created'] += 1

        except FileNotFoundError:
            logging.warning(f'Could not find file: {os.path.join(png_path, dicom_file + png_postfix)}')
            stats['File not found error'] += 1
    for k in stats:
        if sample_header in k and stats[k] == 50:
            continue
        logging.info(f'{k} has {stats[k]}')


def write_tensors_from_ecg_pngs(tensors, png_path, min_sample_id, max_sample_id, name='ecg_segmented', png_postfix='.png.mask.png', path_prefix='ecg_rest'):
    stats = Counter()
    for png_file in os.listdir(png_path):
        sample_id = png_file.split('.')[0]
        offset_seconds = float(png_file.split('_')[-1].replace(png_postfix, ''))
        if not min_sample_id <= int(sample_id) < max_sample_id:
            continue

        png = imageio.imread(os.path.join(png_path, png_file))
        stats[f'png shape {png.shape}'] += 1
        ecg_1d_segmentation = png[0, :, 0]
        logging.info(f' ecg1d segmentation unique: {np.unique(ecg_1d_segmentation)} \nstart {ecg_1d_segmentation[:60]} ')
        stats[f'ecg_1d_segmentation shape {ecg_1d_segmentation.shape}'] += 1
        tensor_file = os.path.join(tensors, str(sample_id) + TENSOR_EXT)
        if not os.path.exists(os.path.dirname(tensor_file)):
            os.makedirs(os.path.dirname(tensor_file))
        with h5py.File(tensor_file, 'a') as hd5:
            tp = tensor_path(path_prefix, name)
            if tp in hd5:
                tensor_dataset = first_dataset_at_path(hd5, tp)
                tensor_dataset[:] = ecg_1d_segmentation
                tensor_dataset.attrs['offset_seconds'] = offset_seconds
                stats['updated'] += 1
            else:
                create_tensor_in_hd5(hd5, path_prefix, name, ecg_1d_segmentation, stats, attributes={'offset_seconds': offset_seconds})
                stats['created'] += 1

    for k in stats:
        logging.info(f'{k} has {stats[k]}')


def _load_meta_data_for_tensor_writing(volume_csv: str, lv_mass_csv: str, min_sample_id: int, max_sample_id: int) -> Tuple[Dict[int, Dict[str, float]], List[int]]:
    """ Gather metadata necessary to write tensors from UK biobank

    Loads the field IDs of survey data, dates of assessment, diagnosis of diseases,
    ejection fractions, diastolic volumes, systolic volumes, and sample IDs to make tensors from.

    :param volume_csv: CSV containing systole and diastole volumes and ejection fraction for samples with MRI
    :param lv_mass_csv: TSV containing left ventricular mass and other cardiac MRI readouts on ~5000 people returned from app 2964
    :param min_sample_id: Minimum sample id to generate, for parallelization
    :param max_sample_id: Maximum sample id to generate, for parallelization
    :return: Tuple of metadata containers
        nested_dictionary: Dictionary mapping sample IDs (as ints) to dictionaries mapping strings to values from the CSV
        sample_ids: List of sample IDs (as ints) to generate tensors for
    """
    nested_dictionary = defaultdict(dict)
    with open(volume_csv, 'r') as volumes:
        lol = list(csv.reader(volumes, delimiter='\t'))
        logging.info(f"CSV of MRI volumes header:{list(enumerate(lol[0]))}")
        fields = lol[0][1:]  # Assumes sample id is the first field
        for row in lol[1:]:
            sample_id = int(row[0])
            if min_sample_id <= sample_id <= max_sample_id:
                nested_dictionary[sample_id] = {fields[i].strip().lower(): row[i+1] for i in range(len(fields))}

    with open(lv_mass_csv, 'r') as lvm:
        lol = list(csv.reader(lvm, delimiter='\t'))
        logging.info('CSV of returned MRI mass, etc, header:{}'.format(list(enumerate(lol[0]))))
        for row in lol[1:]:
            # column 0 is the original app's sample ID. Column 1 is app 7089's sample ID.
            sample_id = int(row[1])
            if min_sample_id <= sample_id <= max_sample_id and len(row[13]) > 0 and row[13] != 'NA':
                # Zero-based column #13 is the LV mass
                nested_dictionary[sample_id]['lv_mass'] = float(row[13])

    return nested_dictionary


def _sample_has_dicom_mris(zip_folder, sample_id) -> bool:
    sample_str = str(sample_id)
    return any([os.path.exists(os.path.join(zip_folder, f'{sample_str}_{mri_f}_2_0.zip')) for mri_f in DICOM_MRI_FIELDS])


def _sample_has_nifti_mris(zip_folder, sample_id) -> bool:
    sample_str = str(sample_id)
    return any([os.path.exists(os.path.join(zip_folder, f'{sample_str}_{mri_f}_2_0.zip')) for mri_f in MRI_NIFTI_FIELD_ID_TO_ROOT])


def _sample_has_ecgs(xml_folder, xml_field_ids, sample_id) -> bool:
    sample_str = str(sample_id)
    for xml_id in xml_field_ids:
        if os.path.exists(xml_folder + sample_str + '_' + xml_id + '_0_0.xml'):
            return True
        if os.path.exists(xml_folder + sample_str + '_' + xml_id + '_1_0.xml'):
            return True
        if os.path.exists(xml_folder + sample_str + '_' + xml_id + '_2_0.xml'):
            return True
    return False


def _dicts_and_plots_from_tensorization(
    a_id: str,
    output_folder: str,
    min_values_to_print: int,
    write_pngs: bool,
    continuous_stats: Dict[str, List[float]],
    stats: Dict[str, int],
) -> None:
    """Print out dictionaries of data encountered during tensorization. Optionally make plots of this data.

    :param a_id: User chosen string to identify this run
    :param output_folder: Folder to write outputs to (mostly for debugging)
    :param min_values_to_print: Minimum number of samples that have responded to question for it to be included in the
            categorical or continuous dictionaries printed after tensor generation
    :param write_pngs: write MRIs as PNG images for debugging
    :param continuous_stats: Dictionary mapping field meanings to the list of continuous values found for them
    :param stats: Dictionary mapping strings to ints keeps track of categorical responses and other info.
    :return: None
    """
    categories = {}
    continuous = {}
    value_counter = Counter()
    for k in sorted(list(stats.keys())):
        logging.info("{} has {}".format(k, stats[k]))

        if 'categorical' not in k and 'continuous' not in k:
            continue
        parts = k.split(JOIN_CHAR)
        if len(parts) > 1:
            column_key = k.replace(HD5_GROUP_CHAR, JOIN_CHAR)  # flatten the hd5 group
            value_counter[column_key] += stats[k]
            if min_values_to_print < value_counter[column_key]:
                if 'categorical' in k and column_key not in categories:
                    categories[column_key.replace('categorical_', '')] = len(categories)
                if 'continuous' in k and column_key not in continuous:
                    continuous[column_key.replace('continuous_', '')] = len(continuous)

    if write_pngs:
        plot_value_counter(list(categories.keys()), value_counter, a_id + '_v_count', os.path.join(output_folder, a_id))
        plot_histograms(continuous_stats, a_id, os.path.join(output_folder, a_id))

    logging.info("Continuous tensor map: {}".format(continuous))
    logging.info("Continuous Columns: {}".format(len(continuous)))
    logging.info("Category tensor map: {}".format(categories))
    logging.info("Categories Columns: {}".format(len(categories)))


def _to_float_or_false(s):
    try:
        return float(s)
    except ValueError:
        return False


def _to_float_or_nan(s):
    try:
        return float(s)
    except ValueError:
        return np.nan


def _write_tensors_from_zipped_dicoms(
    zoom_x: int,
    zoom_y: int,
    zoom_width: int,
    zoom_height: int,
    write_pngs: bool,
    tensors: str,
    dicoms: str,
    mri_field_ids: List[str],
    zip_folder: str,
    hd5: h5py.File,
    sample_id: int,
    stats: Dict[str, int],
) -> None:
    sample_str = str(sample_id)
    for mri_field in set(mri_field_ids).intersection(DICOM_MRI_FIELDS):
        mris = glob.glob(zip_folder + sample_str + '_' + mri_field + '*.zip')
        for zipped in mris:
            logging.info("Got zipped dicoms for sample: {} with MRI field: {}".format(sample_id, mri_field))
            dicom_folder = os.path.join(dicoms, sample_str, mri_field)
            if not os.path.exists(dicom_folder):
                os.makedirs(dicom_folder)
            with zipfile.ZipFile(zipped, "r") as zip_ref:
                zip_ref.extractall(dicom_folder)
                _write_tensors_from_dicoms(
                    zoom_x, zoom_y, zoom_width, zoom_height, write_pngs, tensors, dicom_folder,
                    hd5, sample_str, stats,
                )
                stats['MRI fields written'] += 1
            shutil.rmtree(dicom_folder)


def _write_tensors_from_zipped_niftis(zip_folder: str, mri_field_ids: List[str], hd5: h5py.File, sample_id: str, stats: Dict[str, int]) -> None:
    for mri_field in set(mri_field_ids).intersection(set(MRI_NIFTI_FIELD_ID_TO_ROOT.keys())):
        mris = glob.glob(os.path.join(zip_folder, f'{sample_id}_{mri_field}*.zip'))
        for zipped in mris:
            logging.info(f"Got zipped niftis for sample: {sample_id} with MRI field: {mri_field}")
            with tempfile.TemporaryDirectory() as temp_folder, zipfile.ZipFile(zipped, "r") as zip_ref:
                zip_ref.extractall(temp_folder)
                _write_tensors_from_niftis(temp_folder, hd5, mri_field, stats)
                stats['MRI fields written'] += 1


def _write_tensors_from_dicoms(
    zoom_x: int, zoom_y: int, zoom_width: int, zoom_height: int, write_pngs: bool, tensors: str,
    dicom_folder: str, hd5: h5py.File, sample_str: str, stats: Dict[str, int],
) -> None:
    """Convert a folder of DICOMs from a sample into tensors for each series

    Segmented dicoms require special processing and are written to tensor per-slice

    Arguments
        :param x: Width of the tensors (actual MRI width will be padded with 0s or cropped to this number)
        :param y: Height of the tensors (actual MRI width will be padded with 0s or cropped to this number)
        :param z: Minimum number of slices to include in the each tensor if more slices are found they will be kept
        :param zoom_x: x coordinate of the zoom
        :param zoom_y: y coordinate of the zoom
        :param zoom_width: width of the zoom
        :param zoom_height: height of the zoom
        :param write_pngs: write MRIs as PNG images for debugging
        :param tensors: Folder where hd5 tensor files are being written
        :param dicom_folder: Folder with all dicoms associated with one sample.
        :param hd5: Tensor file in which to create datasets for each series and each segmented slice
        :param sample_str: The current sample ID as a string
        :param stats: Counter to keep track of summary statistics

    """
    views = defaultdict(list)
    min_ideal_series = 9e9
    for dicom in os.listdir(dicom_folder):
        if os.path.splitext(dicom)[-1] != DICOM_EXT:
            continue
        d = pydicom.read_file(os.path.join(dicom_folder, dicom))
        series = d.SeriesDescription.lower().replace(' ', '_')
        if series + '_12bit' in MRI_LIVER_SERIES_12BIT and d.LargestImagePixelValue > 2048:
            views[series + '_12bit'].append(d)
            stats[series + '_12bit'] += 1
        elif series in MRI_LIVER_SERIES + MRI_CARDIAC_SERIES + MRI_BRAIN_SERIES:
            views[series].append(d)
            stats[series] += 1
        if series in MRI_LIVER_IDEAL_PROTOCOL:
            min_ideal_series = min(min_ideal_series, int(d.SeriesNumber))

    for v in views:
        mri_shape = (views[v][0].Rows, views[v][0].Columns, len(views[v]))
        mri_date = _datetime_from_dicom(views[v][0])
        stats[v + ' mri shape:' + str(mri_shape)] += 1
        if v in MRI_BRAIN_SERIES:
            mri_group = 'ukb_brain_mri'
        elif v in MRI_LIVER_SERIES + MRI_LIVER_SERIES_12BIT:
            mri_group = 'ukb_liver_mri'
        elif v in MRI_CARDIAC_SERIES + MRI_CARDIAC_SERIES_SEGMENTED:
            mri_group = 'ukb_cardiac_mri'
        else:
            mri_group = 'ukb_mri'

        if v == MRI_TO_SEGMENT:
            _tensorize_short_and_long_axis_segmented_cardiac_mri(views[v], v, zoom_x, zoom_y, zoom_width, zoom_height, write_pngs, tensors, hd5, mri_date, mri_group, stats)
        elif v in MRI_BRAIN_SERIES:
            _tensorize_brain_mri(views[v], v, mri_date, mri_group, hd5)
        else:
            mri_data = np.zeros((views[v][0].Rows, views[v][0].Columns, len(views[v])), dtype=np.float32)
            for slicer in views[v]:
                _save_pixel_dimensions_if_missing(slicer, v, hd5)
                _save_slice_thickness_if_missing(slicer, v, hd5)
                _save_series_orientation_and_position_if_missing(slicer, v, hd5)
                slice_index = slicer.InstanceNumber - 1
                if v in MRI_LIVER_IDEAL_PROTOCOL:
                    slice_index = _slice_index_from_ideal_protocol(slicer, min_ideal_series)
                mri_data[..., slice_index] = slicer.pixel_array.astype(np.float32)
            create_tensor_in_hd5(hd5, mri_group, v, mri_data, stats, mri_date)


def _tensorize_short_and_long_axis_segmented_cardiac_mri(
    slices: List[pydicom.Dataset], series: str, zoom_x: int, zoom_y: int,
    zoom_width: int, zoom_height: int, write_pngs: bool, tensors: str,
    hd5: h5py.File, mri_date: datetime.datetime, mri_group: str,
    stats: Dict[str, int],
) -> None:
    systoles = {}
    diastoles = {}
    systoles_pix = {}
    systoles_masks = {}
    diastoles_masks = {}

    for slicer in slices:
        full_mask = np.zeros((slicer.Rows, slicer.Columns), dtype=np.float32)
        full_slice = np.zeros((slicer.Rows, slicer.Columns), dtype=np.float32)

        if _has_overlay(slicer):
            if _is_mitral_valve_segmentation(slicer):
                series = series.replace('sax', 'lax')
            else:
                series = series.replace('lax', 'sax')
            series_segmented = f'{series}_segmented'
            series_zoom = f'{series}_zoom'
            series_zoom_segmented = f'{series}_zoom_segmented'

            try:
                overlay, mask, ventricle_pixels, _ = _get_overlay_from_dicom(slicer)
            except KeyError:
                logging.exception(f'Got key error trying to make anatomical mask, skipping.')
                continue

            _save_pixel_dimensions_if_missing(slicer, series, hd5)
            _save_slice_thickness_if_missing(slicer, series, hd5)
            _save_series_orientation_and_position_if_missing(slicer, series, hd5, str(slicer.InstanceNumber))
            _save_pixel_dimensions_if_missing(slicer, series_segmented, hd5)
            _save_slice_thickness_if_missing(slicer, series_segmented, hd5)
            _save_series_orientation_and_position_if_missing(slicer, series_segmented, hd5, str(slicer.InstanceNumber))

            cur_angle = (slicer.InstanceNumber - 1) // MRI_FRAMES  # dicom InstanceNumber is 1-based
            full_slice[:] = slicer.pixel_array.astype(np.float32)
            create_tensor_in_hd5(hd5, mri_group, f'{series}{HD5_GROUP_CHAR}{slicer.InstanceNumber}', full_slice, stats, mri_date)
            create_tensor_in_hd5(hd5, mri_group, f'{series_zoom_segmented}{HD5_GROUP_CHAR}{slicer.InstanceNumber}', mask, stats, mri_date)

            zoom_slice = full_slice[zoom_x: zoom_x + zoom_width, zoom_y: zoom_y + zoom_height]
            zoom_mask = mask[zoom_x: zoom_x + zoom_width, zoom_y: zoom_y + zoom_height]
            create_tensor_in_hd5(hd5, mri_group, f'{series_zoom}{HD5_GROUP_CHAR}{slicer.InstanceNumber}', zoom_slice, stats, mri_date)
            create_tensor_in_hd5(hd5, mri_group, f'{series_zoom_segmented}{HD5_GROUP_CHAR}{slicer.InstanceNumber}', zoom_mask, stats, mri_date)

            if (slicer.InstanceNumber - 1) % MRI_FRAMES == 0:  # Diastole frame is always the first
                diastoles[cur_angle] = slicer
                diastoles_masks[cur_angle] = mask
            if cur_angle not in systoles:
                systoles[cur_angle] = slicer
                systoles_pix[cur_angle] = ventricle_pixels
                systoles_masks[cur_angle] = mask
            else:
                if ventricle_pixels < systoles_pix[cur_angle]:
                    systoles[cur_angle] = slicer
                    systoles_pix[cur_angle] = ventricle_pixels
                    systoles_masks[cur_angle] = mask

    for angle in diastoles:
        logging.info(f'Found systole, instance:{systoles[angle].InstanceNumber} ventricle pixels:{systoles_pix[angle]}')
        full_slice = diastoles[angle].pixel_array.astype(np.float32)
        create_tensor_in_hd5(hd5, mri_group, f'diastole_frame_b{angle}', full_slice, stats, mri_date)
        create_tensor_in_hd5(hd5, mri_group, f'diastole_mask_b{angle}', diastoles_masks[angle], stats, mri_date)
        if write_pngs:
            plt.imsave(tensors + 'diastole_frame_b' + str(angle) + IMAGE_EXT, full_slice)
            plt.imsave(tensors + 'diastole_mask_b' + str(angle) + IMAGE_EXT, full_mask)

        full_slice = systoles[angle].pixel_array.astype(np.float32)
        create_tensor_in_hd5(hd5, mri_group, f'systole_frame_b{angle}', full_slice, stats, mri_date)
        create_tensor_in_hd5(hd5, mri_group, f'systole_mask_b{angle}', systoles_masks[angle], stats, mri_date)
        if write_pngs:
            plt.imsave(tensors + 'systole_frame_b' + str(angle) + IMAGE_EXT, full_slice)
            plt.imsave(tensors + 'systole_mask_b' + str(angle) + IMAGE_EXT, full_mask)


def _tensorize_brain_mri(slices: List[pydicom.Dataset], series: str, mri_date: datetime.datetime, mri_group: str, hd5: h5py.File) -> None:
    mri_data1 = np.zeros((slices[0].Rows, slices[0].Columns, len(slices) // 2), dtype=np.float32)
    mri_data2 = np.zeros((slices[0].Rows, slices[0].Columns, len(slices) // 2), dtype=np.float32)
    for slicer in slices:
        _save_pixel_dimensions_if_missing(slicer, series, hd5)
        _save_slice_thickness_if_missing(slicer, series, hd5)
        _save_series_orientation_and_position_if_missing(slicer, series, hd5)
        slice_index = slicer.InstanceNumber - 1
        if slicer.SeriesNumber in [5, 11]:
            mri_data1[..., slice_index] = slicer.pixel_array.astype(np.float32)
        elif slicer.SeriesNumber in [6, 12]:
            mri_data2[..., slice_index] = slicer.pixel_array.astype(np.float32)
    create_tensor_in_hd5(hd5, mri_group, series + '_1', mri_data1, date=mri_date)
    create_tensor_in_hd5(hd5, mri_group, series + '_2', mri_data2, date=mri_date)


def _save_pixel_dimensions_if_missing(slicer, series, hd5):
    if MRI_PIXEL_WIDTH + '_' + series not in hd5 and series in MRI_BRAIN_SERIES + MRI_CARDIAC_SERIES + MRI_CARDIAC_SERIES_SEGMENTED + MRI_LIVER_SERIES + MRI_LIVER_SERIES_12BIT:
        hd5.create_dataset(MRI_PIXEL_WIDTH + '_' + series, data=float(slicer.PixelSpacing[0]))
    if MRI_PIXEL_HEIGHT + '_' + series not in hd5 and series in MRI_BRAIN_SERIES + MRI_CARDIAC_SERIES + MRI_CARDIAC_SERIES_SEGMENTED + MRI_LIVER_SERIES + MRI_LIVER_SERIES_12BIT:
        hd5.create_dataset(MRI_PIXEL_HEIGHT + '_' + series, data=float(slicer.PixelSpacing[1]))


def _save_slice_thickness_if_missing(slicer, series, hd5):
    if MRI_SLICE_THICKNESS + '_' + series not in hd5 and series in MRI_BRAIN_SERIES + MRI_CARDIAC_SERIES + MRI_CARDIAC_SERIES_SEGMENTED + MRI_LIVER_SERIES + MRI_LIVER_SERIES_12BIT:
        hd5.create_dataset(MRI_SLICE_THICKNESS + '_' + series, data=float(slicer.SliceThickness))


def _save_series_orientation_and_position_if_missing(slicer, series, hd5, instance=None):
    orientation_ds_name = MRI_PATIENT_ORIENTATION + '_' + series
    position_ds_name = MRI_PATIENT_POSITION + '_' + series
    if instance:
        orientation_ds_name += HD5_GROUP_CHAR + instance
        position_ds_name += HD5_GROUP_CHAR + instance
    if orientation_ds_name not in hd5 and series in MRI_BRAIN_SERIES + MRI_CARDIAC_SERIES + MRI_CARDIAC_SERIES_SEGMENTED + MRI_LIVER_SERIES + MRI_LIVER_SERIES_12BIT:
        hd5.create_dataset(orientation_ds_name, data=[float(x) for x in slicer.ImageOrientationPatient])
    if position_ds_name not in hd5 and series in MRI_BRAIN_SERIES + MRI_CARDIAC_SERIES + MRI_CARDIAC_SERIES_SEGMENTED + MRI_LIVER_SERIES + MRI_LIVER_SERIES_12BIT:
        hd5.create_dataset(position_ds_name, data=[float(x) for x in slicer.ImagePositionPatient])


def _has_overlay(d) -> bool:
    try:
        _ = d[0x6000, 0x3000].value
        return True
    except KeyError:
        return False


def _is_mitral_valve_segmentation(d) -> bool:
    return d.SliceThickness == 6


def _slice_index_from_ideal_protocol(d, min_ideal_series):
    return 6*(d.InstanceNumber-1) + ((d.SeriesNumber-min_ideal_series)//2)


def _get_overlay_from_dicom(d, debug=False) -> Tuple[np.ndarray, np.ndarray]:
    """Get an overlay from a DICOM file

    Morphological operators are used to transform the pixel outline of the myocardium
    to the labeled pixel masks for myocardium and left ventricle

    Arguments
        d: the dicom file
        stats: Counter to keep track of summary statistics

    Returns
        Tuple of two numpy arrays.
        The first is the raw overlay array with myocardium outline,
        The second is a pixel mask with 0 for background 1 for myocardium and 2 for ventricle
    """
    i_overlay = 0
    dicom_tag = 0x6000 + 2 * i_overlay
    overlay_raw = d[dicom_tag, 0x3000].value
    rows = d[dicom_tag, 0x0010].value  # rows = 512
    cols = d[dicom_tag, 0x0011].value  # cols = 512
    overlay_frames = d[dicom_tag, 0x0015].value
    bits_allocated = d[dicom_tag, 0x0100].value

    np_dtype = np.dtype('uint8')
    length_of_pixel_array = len(overlay_raw)
    expected_length = rows * cols
    if bits_allocated == 1:
        expected_bit_length = expected_length
        bit = 0
        overlay = np.ndarray(shape=(length_of_pixel_array * 8), dtype=np_dtype)
        for byte in overlay_raw:
            for bit in range(bit, bit + 8):
                overlay[bit] = byte & 0b1
                byte >>= 1
            bit += 1
        overlay = overlay[:expected_bit_length]
    if overlay_frames == 1:
        overlay = overlay.reshape(rows, cols)
        idx = np.where(overlay == 1)
        min_pos = (np.min(idx[0]), np.min(idx[1]))
        max_pos = (np.max(idx[0]), np.max(idx[1]))
        short_side = min((max_pos[0] - min_pos[0]), (max_pos[1] - min_pos[1]))
        small_radius = max(MRI_MIN_RADIUS, short_side * MRI_SMALL_RADIUS_FACTOR)
        big_radius = max(MRI_MIN_RADIUS+1, short_side * MRI_BIG_RADIUS_FACTOR)
        small_structure = _unit_disk(small_radius)
        m1 = binary_closing(overlay, small_structure).astype(np.int)
        big_structure = _unit_disk(big_radius)
        m2 = binary_closing(overlay, big_structure).astype(np.int)
        anatomical_mask = m1 + m2
        ventricle_pixels = np.count_nonzero(anatomical_mask == MRI_SEGMENTED_CHANNEL_MAP['ventricle'])
        myocardium_pixels = np.count_nonzero(anatomical_mask == MRI_SEGMENTED_CHANNEL_MAP['myocardium'])
        if ventricle_pixels == 0 and myocardium_pixels > MRI_MAX_MYOCARDIUM:  # try to rescue small ventricles
            erode_structure = _unit_disk(small_radius*1.5)
            anatomical_mask = anatomical_mask - binary_erosion(m1, erode_structure).astype(np.int)
            ventricle_pixels = np.count_nonzero(anatomical_mask == MRI_SEGMENTED_CHANNEL_MAP['ventricle'])
            myocardium_pixels = np.count_nonzero(anatomical_mask == MRI_SEGMENTED_CHANNEL_MAP['myocardium'])
        return overlay, anatomical_mask, ventricle_pixels, myocardium_pixels


def _unit_disk(r) -> np.ndarray:
    y, x = np.ogrid[-r: r + 1, -r: r + 1]
    return (x ** 2 + y ** 2 <= r ** 2).astype(np.int)


def _outline_to_mask(labeled_outline, idx) -> np.ndarray:
    idx = np.where(labeled_outline == idx)
    poly = list(zip(idx[1].tolist(), idx[0].tolist()))
    img = Image.new("L", [labeled_outline.shape[1], labeled_outline.shape[0]], 0)
    ImageDraw.Draw(img).polygon(poly, outline=1, fill=1)
    return np.array(img)


def _write_tensors_from_xml(xml_field_ids, xml_folder, hd5, sample_id, write_pngs, stats, continuous_stats) -> None:
    for xml_field in xml_field_ids:
        xmlp = xml_folder + str(sample_id) + '_' + xml_field + '*.xml'
        ecgs = glob.glob(xmlp)
        if xml_field == ECG_REST_FIELD:
            _write_ecg_rest_tensors(ecgs, xml_field, hd5, sample_id, write_pngs, stats, continuous_stats)
        elif xml_field == ECG_BIKE_FIELD:
            _write_ecg_bike_tensors(ecgs, xml_field, hd5, sample_id, stats)
        else:
            raise ValueError('Unknown ECG field ID:', xml_field)


def _write_ecg_rest_tensors(ecgs, xml_field, hd5, sample_id, write_pngs, stats, continuous_stats) -> None:
    rest_group = 'ukb_ecg_rest'
    for ecg in ecgs:
        logging.info('Got ECG for sample:{} XML field:{}'.format(sample_id, xml_field))
        root = et.parse(ecg).getroot()
        ecg_date = _str2date(_date_str_from_ecg(root))
        diagnosis_text = []
        for d in root.findall("./Interpretation/Diagnosis/DiagnosisText"):
            if 'QRS Complexes:' in d.text:
                qrs = float(d.text.replace('QRS Complexes:', '').split(',')[0].strip())
                create_tensor_in_hd5(hd5, rest_group, 'QRSComplexes', qrs, stats, date=ecg_date)
            elif '---' in d.text or 'Arrhythmia results of the full-disclosure ECG' in d.text:
                continue
            else:
                diagnosis_text.append(d.text.replace(',', '').replace('*', '').replace('&', 'and').replace('  ', ' '))

        diagnosis_str = ' '.join(diagnosis_text)
        create_tensor_in_hd5(hd5, rest_group, 'ecg_rest_text', diagnosis_str, stats, date=ecg_date, storage_type=StorageType.STRING)

        for c in root.findall("./StripData/WaveformData"):
            lead_data = list(map(float, c.text.strip().split(',')))
            dataset_name = 'strip_' + str(c.attrib['lead'])
            create_tensor_in_hd5(hd5, rest_group, dataset_name, lead_data, stats, date=ecg_date)
            stats[dataset_name] += 1

        for c in root.findall("./RestingECGMeasurements"):
            for child in c:
                if child.text is not None and child.tag in ECG_TAGS_TO_WRITE:
                    create_tensor_in_hd5(hd5, rest_group, child.tag, float(child.text), stats, date=ecg_date)
                    stats[child.tag] += 1
                    if write_pngs:
                        continuous_stats[child.tag].append(float(child.text))
                if child.tag == 'MedianSamples':
                    for median_c in child:
                        if median_c.tag == 'WaveformData':
                            median_wave = list(map(float, median_c.text.strip().split(',')))
                            dataset_name = 'median_' + str(median_c.attrib['lead'])
                            create_tensor_in_hd5(hd5, 'ukb_ecg_rest', dataset_name, median_wave, stats, date=ecg_date)
        for c in root.findall("./RestingECGMeasurements/MeasurementTable"):
            for child in c:
                if child.tag not in ECG_TABLE_TAGS:
                    continue
                values = list(map(_to_float_or_nan, child.text.strip().split(',')))
                create_tensor_in_hd5(hd5, 'ukb_ecg_rest', child.tag.lower(), values, stats, date=ecg_date)


def create_tensor_in_hd5(
    hd5: h5py.File, path_prefix: str, name: str, value, stats: Counter = None, date: datetime.datetime = None,
    storage_type: StorageType = None, attributes: Dict[str, Any] = None,
):
    hd5_path = tensor_path(path_prefix, name)
    if hd5_path in hd5:
        hd5_path = f'{hd5_path}instance_{len(hd5[hd5_path])}'
    else:
        hd5_path = f'{hd5_path}instance_0'
    if stats is not None:
        stats[hd5_path] += 1
    if storage_type == StorageType.STRING:
        d = hd5.create_dataset(hd5_path, data=value, dtype=h5py.special_dtype(vlen=str))
    elif isinstance(value, int) or isinstance(value, float):
        d = hd5.create_dataset(hd5_path, data=[value])
    elif isinstance(value, np.ndarray) or isinstance(value, list):
        d = hd5.create_dataset(hd5_path, data=value, compression='gzip')
    else:
        raise NotImplementedError(f'{storage_type} cannot be automatically written yet')  # TODO: Add categorical, etc.
    if date is not None:
        d.attrs['date'] = time.mktime(date.timetuple())
    if attributes is not None:
        for k in attributes:
            d.attrs[k] = attributes[k]


def tensor_path(path_prefix: str, name: str) -> str:
    """
    In the future, TMAPs should be generated using this same function
    """
    return f'/{path_prefix}/{name}/'


def first_dataset_at_path(hd5, path, gather_fxn=min):
    if path not in hd5:
        raise ValueError(f'Could not find key:{path} in hd5.')
    data = hd5[path]
    if isinstance(data, h5py.Dataset):
        return data
    deeper_key_prefix = f'{path}{gather_fxn(hd5[path])}/'
    return first_dataset_at_path(hd5, deeper_key_prefix)


def _datetime_to_str(dt: datetime.datetime) -> str:
    return dt.strftime(DATE_FORMAT)


def path_date_to_datetime(date: str) -> datetime.datetime:
    return datetime.datetime.strptime(date, DATE_FORMAT)


def _write_ecg_bike_tensors(ecgs, xml_field, hd5, sample_id, stats):
    for ecg in ecgs:
        root = et.parse(ecg).getroot()
        date = datetime.datetime.strptime(_date_str_from_ecg(root), '%Y-%m-%d')
        write_to_hd5 = partial(create_tensor_in_hd5, hd5=hd5, path_prefix='ukb_ecg_bike', stats=stats, date=date)
        logging.info('Got ECG for sample:{} XML field:{}'.format(sample_id, xml_field))

        instance = ecg.split(JOIN_CHAR)[-2]
        write_to_hd5(storage_type=StorageType.STRING, name='instance', value=instance)

        protocol = root.findall('./Protocol/Phase')[0].find('ProtocolName').text
        write_to_hd5(storage_type=StorageType.STRING, name='protocol', value=protocol)

        median_ecgs = defaultdict(list)
        for median_waves in root.findall('./MedianData/Median/WaveformData'):
            median_ecgs[median_waves.attrib['lead']].extend(list(map(float, median_waves.text.strip().split(','))))
        if len(median_ecgs) > 0:
            median_np = np.zeros(ECG_BIKE_MEDIAN_SIZE)
            for lead in median_ecgs:
                median_idx = min(ECG_BIKE_MEDIAN_SIZE[0], len(median_ecgs[lead]))
                median_np[:median_idx, ECG_BIKE_LEADS[lead]] = median_ecgs[lead][:median_idx]
            write_to_hd5(name='median', value=median_np)
        else:
            stats['missing median bike ECG'] += 1

        counter = 0
        strip_np = np.zeros(ECG_BIKE_STRIP_SIZE)
        for strip_waves in root.findall("./StripData/Strip/WaveformData"):
            counter += 1
            strip_list = list(map(float, strip_waves.text.strip().split(',')))
            strip_idx = min(ECG_BIKE_MEDIAN_SIZE[0], len(strip_list))
            strip_np[:strip_idx, ECG_BIKE_LEADS[strip_waves.attrib['lead']]] = strip_list[:strip_idx]
        if counter > 0:
            write_to_hd5(name='strip', value=strip_np)
        else:
            stats['missing strip bike ECG'] += 1

        full_ekgs = [[] for _ in range(ECG_BIKE_NUM_LEADS)]
        count = 0
        for full_d in root.findall("./FullDisclosure/FullDisclosureData"):
            for full_line in re.split('\n|\t', full_d.text):
                for sample in re.split(',', full_line):
                    if sample == '':
                        continue
                    lead = (count % (ECG_BIKE_NUM_LEADS * ECG_BIKE_SAMPLE_RATE)) // ECG_BIKE_SAMPLE_RATE
                    full_ekgs[lead].append(float(sample))
                    count += 1

        if all(full_ekgs):  # Does each lead have data?
            full_np = np.zeros(ECG_BIKE_FULL_SIZE)
            for i, lead in enumerate(full_ekgs):
                full_idx = min(ECG_BIKE_FULL_SIZE[0], len(lead))
                full_np[:full_idx, i] = lead[:full_idx]
            write_to_hd5(name='full', value=full_np)
        else:
            stats['missing full disclosure bike ECG'] += 1

        # Patient info
        patient_fields = ('Age', 'Height', 'Weight')
        for field in patient_fields:
            val = [_xml_path_to_float(root, f'./PatientInfo/{field}')]
            write_to_hd5(name=str.lower(field), value=val)

        # Trend measurements
        trend_entry_fields = ['HeartRate', 'Load', 'Grade', 'Mets', 'VECount', 'PaceCount']
        trend_lead_measurements = [
            'JPointAmplitude',
            'STAmplitude20ms',
            'STAmplitude',
            'RAmplitude',
            'R1Amplitude',
            'STSlope',
            'STIntegral',
        ]
        phase_to_int = {'Pretest': 0, 'Exercise': 1, 'Rest': 2}
        lead_to_int = {'I': 0, '2': 1, '3': 2}
        trend_entries = root.findall("./TrendData/TrendEntry")

        trends = {trend: np.full(len(trend_entries), np.nan) for trend in trend_entry_fields + ['time', 'PhaseTime', 'PhaseTime', 'PhaseName', 'Artifact']}
        trends.update({trend: np.full((len(trend_entries), 3), np.nan) for trend in trend_lead_measurements})

        for i, trend_entry in enumerate(trend_entries):
            for field in trend_entry_fields:
                field_val = trend_entry.find(field)
                if field_val is None:
                    continue
                field_val = _to_float_or_false(field_val.text)
                if field_val is False:
                    continue
                trends[field][i] = field_val
            for lead_field, lead in product(trend_lead_measurements, trend_entry.findall('LeadMeasurements')):
                lead_num = lead.attrib['lead']
                field_val = lead.find(lead_field)
                if field_val is None:
                    continue
                field_val = _to_float_or_false(field_val.text.strip('?'))
                if field_val is False:
                    continue
                trends[lead_field][i, lead_to_int[lead_num]] = field_val
            trends['time'][i] = SECONDS_PER_MINUTE * int(trend_entry.find("EntryTime/Minute").text) + int(trend_entry.find("EntryTime/Second").text)
            trends['PhaseTime'][i] = SECONDS_PER_MINUTE * int(trend_entry.find("PhaseTime/Minute").text) + int(trend_entry.find("PhaseTime/Second").text)
            trends['PhaseName'][i] = phase_to_int[trend_entry.find('PhaseName').text]
            trends['Artifact'][i] = float(trend_entry.find('Artifact').text.strip('%')) / 100  # Artifact is reported as a percentage

        for field, trend_list in trends.items():
            write_to_hd5(name=f'trend_{str.lower(field)}', value=trend_list)

        # Last 60 seconds of raw given that the rest phase is 60s
        phase_durations = {}
        for protocol in root.findall("./Protocol/Phase"):
            phase_name = protocol.find("PhaseName").text
            phase_duration = SECONDS_PER_MINUTE * int(protocol.find("PhaseDuration/Minute").text) + int(
                protocol.find("PhaseDuration/Second").text,
            )
            phase_durations[phase_name] = phase_duration
            write_to_hd5(name=f'{str.lower(phase_name)}_duration', value=[phase_duration])

        # HR stats
        max_hr = _xml_path_to_float(root, './ExerciseMeasurements/MaxHeartRate')
        resting_hr = _xml_path_to_float(root, './ExerciseMeasurements/RestingStats/RestHR')
        max_pred_hr = _xml_path_to_float(root, './ExerciseMeasurements/MaxPredictedHR')
        write_to_hd5(name='max_hr', value=[max_hr])
        write_to_hd5(name='resting_hr', value=[resting_hr])
        write_to_hd5(name='max_pred_hr', value=[max_pred_hr])


def _write_tensors_from_niftis(folder: str, hd5: h5py.File, field_id: str, stats: Counter):
    niftis = glob.glob(os.path.join(folder, MRI_NIFTI_FIELD_ID_TO_ROOT[field_id], '**/*nii.gz'), recursive=True)
    logging.info(f'Found {len(niftis)} NIFTI files at {os.path.join(folder, MRI_NIFTI_FIELD_ID_TO_ROOT[field_id])} ')
    for nifti in niftis:  # iterate through all nii.gz files and add them to the hd5
        nifti_mri = nib.load(nifti)
        nifti_array = nifti_mri.get_fdata()
        nii_name = os.path.basename(nifti).replace('.nii.gz', '')  # removes .nii.gz
        parent = nifti.split('/')[-2]
        if parent not in nii_name:
            nii_name = f'{parent}_{nii_name}'
        stats[nii_name] += 1
        stats[f'{nii_name} shape:{nifti_array.shape}'] += 1
        if MRI_NIFTI_FIELD_ID_TO_ROOT[field_id] == 'SWI' and nifti_array.shape[-1] == MRI_SWI_SLICES_TO_AXIS_SHIFT and nifti_array.shape[0] > nifti_array.shape[1]:
            nifti_array = np.moveaxis(nifti_array, 0, 1)
            stats[f'{nii_name} shape post SWI shift:{nifti_array.shape}'] += 1
        create_tensor_in_hd5(hd5=hd5, path_prefix='ukb_brain_mri', name=nii_name, value=nifti_array, stats=stats)


def _xml_path_to_float(root: et, path: str) -> float:
    return float(root.find(path).text)


def _date_str_from_ecg(root):
    date_str = ''
    for d in root.findall('./ObservationDateTime/Year'):
        date_str += d.text + '-'
    for d in root.findall('./ObservationDateTime/Month'):
        date_str += d.text + '-'
    for d in root.findall('./ObservationDateTime/Day'):
        date_str += d.text
    return date_str


def _get_disease_censor_dates(disease2tsv, min_sample_id=0, max_sample_id=9e9) -> Dict[str, Dict[int, List[datetime.datetime]]]:
    censor_dates = defaultdict(dict)
    for d in disease2tsv:
        with open(disease2tsv[d], 'r') as my_tsv:
            lol = list(csv.reader(my_tsv, delimiter='\t'))
            for row in lol[1:]:
                sample_id = int(row[0])
                if min_sample_id <= sample_id <= max_sample_id:
                    censor_dates[d][sample_id] = _str2date(row[4])
    return censor_dates


def disease_censor_status(disease2tsv, min_sample_id=0, max_sample_id=9e9) -> Dict[str, Dict[int, List[int]]]:
    censor_status = defaultdict(dict)
    for d in disease2tsv:
        with open(disease2tsv[d], 'r') as my_tsv:
            lol = list(csv.reader(my_tsv, delimiter='\t'))
            for row in lol[1:]:
                sample_id = int(row[0])
                if min_sample_id <= sample_id <= max_sample_id:
                    censor_status[d][sample_id] = 0 if row[1] == 'NA' else int(row[1])
    return censor_status


def disease_prevalence_status(disease2tsv, min_sample_id=0, max_sample_id=9e9) -> Dict[str, Dict[int, List[int]]]:
    censor_status = defaultdict(dict)
    for d in disease2tsv:
        with open(disease2tsv[d], 'r') as my_tsv:
            lol = list(csv.reader(my_tsv, delimiter='\t'))
            for row in lol[1:]:
                sample_id = int(row[0])
                if min_sample_id <= sample_id <= max_sample_id:
                    censor_status[d][sample_id] = 0 if row[2] == 'NA' else int(row[2])
    return censor_status


def disease_incidence_status(disease2tsv, min_sample_id=0, max_sample_id=9e9) -> Dict[str, Dict[int, List[int]]]:
    censor_status = defaultdict(dict)
    for d in disease2tsv:
        with open(disease2tsv[d], 'r') as my_tsv:
            lol = list(csv.reader(my_tsv, delimiter='\t'))
            for row in lol[1:]:
                sample_id = int(row[0])
                if min_sample_id <= sample_id <= max_sample_id:
                    censor_status[d][sample_id] = 0 if row[3] == 'NA' else int(row[3])
    return censor_status


def get_disease2tsv(tsv_folder) -> Dict[str, str]:
    ukb_prefix = 'ukb9222_'
    ukb_postfix = '_phenoV1.tsv'
    disease2tsv = {}
    tsvs = os.listdir(tsv_folder)
    for tsv in tsvs:
        disease_name = tsv.replace(ukb_prefix, '').replace(ukb_postfix, '').lower()
        disease2tsv[disease_name] = tsv_folder + tsv
    return disease2tsv


def append_fields_from_csv(tensors, csv_file, group, delimiter):
    stats = Counter()
    data_maps = defaultdict(dict)
    categorical_channel_maps = {MRI_ANNOTATION_NAME: MRI_ANNOTATION_CHANNEL_MAP}  # str: dict
    with open(csv_file, 'r') as volumes:
        lol = list(csv.reader(volumes, delimiter=delimiter))
        fields = lol[0][1:]  # Assumes sample id is the first field
        logging.info(f"CSV has {len(fields)} fields:{fields}")
        for row in lol[1:]:
            sample_id = row[0]
            data_maps[sample_id] = {fields[i]: row[i+1] for i in range(len(fields))}

    logging.info(f"Data maps:{len(data_maps)} for group:{group}")
    for tp in os.listdir(tensors):
        if os.path.splitext(tp)[-1].lower() != TENSOR_EXT:
            continue
        try:
            with h5py.File(tensors + tp, 'a') as hd5:
                sample_id = tp.replace(TENSOR_EXT, '')
                if sample_id in data_maps:
                    stats['Samples found'] += 1
                    if stats['Samples found'] % 250 == 0:
                        for k in stats:
                            logging.info(f'{k} has: {stats[k]}')
                    for field in data_maps[sample_id]:
                        value = data_maps[sample_id][field]
                        if group == 'continuous':
                            try:
                                value = float(value.strip())
                            except ValueError:
                                stats[f'could not cast field: {field} with value: {value} to float'] += 1
                                continue
                        elif group == 'categorical':
                            is_channel_mapped = False
                            for cm_name in categorical_channel_maps:
                                if value in categorical_channel_maps[cm_name]:
                                    hd5_key = group + HD5_GROUP_CHAR + cm_name
                                    if cm_name in hd5[group]:
                                        data = hd5[hd5_key]
                                        data[0] = categorical_channel_maps[cm_name][value]
                                        stats['updated'] += 1
                                    else:
                                        hd5.create_dataset(hd5_key, data=[categorical_channel_maps[cm_name][value]])
                                        stats['created'] += 1
                                    is_channel_mapped = True
                            if is_channel_mapped:
                                continue

                            if value.lower() in ['false', 'f', '0']:
                                value = 0
                            elif value.lower() in ['true', 't', '1']:
                                value = 1
                            else:
                                stats[f'Could not parse categorical field: {field} with value: {value}'] += 1
                                continue

                        hd5_key = group + HD5_GROUP_CHAR + field
                        if field in hd5[group]:
                            data = hd5[hd5_key]
                            data[0] = value
                            stats['updated'] += 1
                        else:
                            hd5.create_dataset(hd5_key, data=[value])
                            stats['created'] += 1
                else:
                    stats['sample id missing']
        except:
            print('could not open', tp, traceback.format_exc())
            stats['failed'] += 1

    logging.info(f'Finished appending data!')
    for k in stats:
        logging.info(f'{k} has: {stats[k]}')


def append_gene_csv(tensors, csv_file, delimiter):
    stats = Counter()
    data_maps = defaultdict(dict)
    with open(csv_file, 'r') as volumes:
        lol = list(csv.reader(volumes, delimiter=delimiter))
        fields = lol[0][1:]  # Assumes sample id is the first field
        logging.info(f"CSV of flag data header:{fields}")
        for row in lol[1:]:
            sample_id = row[0]
            data_maps[sample_id] = {fields[i]: row[i+1] for i in range(len(fields))}

    logging.info(f"Data maps:{len(data_maps)}")
    hd5_prefix = 'categorical' + HD5_GROUP_CHAR
    for tp in os.listdir(tensors):
        if os.path.splitext(tp)[-1].lower() != TENSOR_EXT:
            continue
        try:
            with h5py.File(tensors + tp, 'a') as hd5:
                sample_id = tp.replace(TENSOR_EXT, '')
                if sample_id in data_maps:
                    for field in data_maps[sample_id]:
                        if field in hd5[hd5_prefix]:
                            data = hd5[data_maps[sample_id][field]]
                            data[0] = 1.0
                            stats['updated'] += 1
                        else:
                            hd5.create_dataset(hd5_prefix + field, data=[1.0])
                            stats['created'] += 1
                else:
                    stats['sample id missing']
        except:
            print('couldnt open', tp, traceback.format_exc())
            stats['failed'] += 1

    for k in stats:
        logging.info("{}: {}".format(k, stats[k]))


# TODO Use 'with' or explicitly close files opened in this method
def _ukbb_stats(run_id, output_folder, phenos_folder, volume_csv, icd_csv, app_csv, zip_folder) -> None:
    stats = Counter()

    steve2sek = {}
    lol = list(csv.reader(open(app_csv, 'r'), delimiter=','))
    for row in lol[1:]:
        steve2sek[row[0]] = row[1]

    disease2tsv = get_disease2tsv(phenos_folder)
    logging.info('got disease tsvs:{}'.format(disease2tsv))
    dates = _get_disease_censor_dates(disease2tsv, 1000000, 2000000)
    status = disease_censor_status(disease2tsv, 1000000, 2000000)
    lol = list(csv.reader(open(volume_csv, 'r'), delimiter='\t'))
    logging.info(list(enumerate(lol[0])))
    lvesv = {}
    lvedv = {}
    lvef = {}
    for row in lol[1:]:
        sample_id = int(row[0])
        lvesv[sample_id] = float(row[2])
        lvedv[sample_id] = float(row[4])
        lvef[sample_id] = float(row[6])

    lol = list(csv.reader(open(icd_csv, 'r'), delimiter='\t'))
    logging.info('CSV of ICDs header:{}'.format(list(enumerate(lol[0]))))

    icd_indexes = {}
    icds = defaultdict(dict)
    for i, header in enumerate(lol[0]):
        if header.lower() in disease2tsv:
            icd_indexes[header.lower()] = i
    logging.info('got icd indexes:{}'.format(icd_indexes))

    disease_dates = defaultdict(list)
    for row in lol[1:]:
        sample_id = int(row[0])
        if os.path.exists(zip_folder + row[0] + '_20209_2_0.zip') \
                or os.path.exists(zip_folder + row[0] + '_20209_0_0.zip') \
                or os.path.exists(zip_folder + row[0] + '_20209_1_0.zip'):
            stats['sample_with_sax_mri'] += 1
        for disease in disease2tsv:
            if sample_id not in dates[disease]:
                stats['sample_id but no status or dates '] += 1
                continue
            if disease not in icd_indexes:
                continue
            icds[disease][sample_id] = row[icd_indexes[disease]]
            if icds[disease][sample_id] != 'NA' and sample_id in status[disease]:
                if status[disease][sample_id] != int(icds[disease][sample_id]):
                    stats[
                        disease + '_icd_disparity_' + icds[disease][sample_id] + '_censor_' + str(
                        status[disease][sample_id],
                        )
                    ] += 1
                else:
                    stats[disease + '_icd_agree'] += 1

            icds[disease][sample_id] = status[disease][sample_id]
            if icds[disease][sample_id] != 'NA' and int(icds[disease][sample_id]) == 1:
                stats[disease + '_total'] += 1
                if sample_id in dates[disease]:
                    stats[disease + '_with_date'] += 1
                    if dates[disease][sample_id] > _str2date('1940-11-11'):
                        disease_dates[disease].append(dates[disease][sample_id])
                if sample_id in lvef:
                    stats[disease + '_with_mri'] += 1
                if sample_id in lvef and sample_id in dates[disease]:
                    stats[disease + '_with_mri_and_date'] += 1

    for k in sorted(list(stats.keys())):
        logging.info('{} has: {}'.format(k, stats[k]))

    logging.info('Plot dates for diseases:{}'.format(list(disease_dates.keys())))
    for i, d in enumerate(disease_dates):
        figure_path = os.path.join(output_folder, run_id, d + '_dates' + IMAGE_EXT)
        if not os.path.exists(os.path.dirname(figure_path)):
            os.makedirs(os.path.dirname(figure_path))
        logging.info('Try to plot figure at:{} with {} dates.'.format(figure_path, len(disease_dates[d])))
        plt.figure(figsize=(16, 16))
        plt.title(d)
        plt.hist(disease_dates[d], bins=60)
        plt.savefig(figure_path)
        plt.close()

    # rows = max(1, int(math.ceil(math.sqrt(len(disease_dates)))))
    # cols = math.ceil(len(disease_dates)/rows)
    # fig, axes = plt.subplots(rows, cols, figsize=(28,28))
    # for i,d in enumerate(disease_dates):
    #     ax = plt.subplot(rows, cols, i+1)
    #     ax.set_title(d)
    #     ax.hist(disease_dates[d], bins=60)
    # plt.savefig(output_folder + 'disease_dates.png')

    _log_extreme_n(lvef, 3)
    _log_extreme_n(lvedv, 3)
    _log_extreme_n(lvesv, 3)


def _print_disease_tensor_maps(phenos_folder) -> None:
    disease2tsv = get_disease2tsv(phenos_folder)
    status = disease_censor_status(disease2tsv)
    for d in sorted(list(disease2tsv.keys())):
        total = len(status[d])
        diseased = np.sum(list(status[d].values()))
        factor = int(total / (diseased * 2))
        print(
            "'{}': TensorMap('{}', group='categorical_index', channel_map={{'no_{}':0, '{}':1}}, loss=weighted_crossentropy([1.0, {}], '{}')),".format(
            d, d, d, d, factor, d,
            ),
        )


def _print_disease_tensor_maps_incident_prevalent(phenos_folder) -> None:
    disease2tsv = get_disease2tsv(phenos_folder)
    status_p = disease_prevalence_status(disease2tsv, 1000000, 2000000)
    status_i = disease_incidence_status(disease2tsv, 1000000, 2000000)
    for d in sorted(list(disease2tsv.keys())):
        total = len(status_p[d])
        diseased_p = np.sum(list(status_p[d].values()))
        factor_p = int(total / (1 + (diseased_p * 3)))
        diseased_i = np.sum(list(status_i[d].values()))
        factor_i = int(total / (1 + (diseased_i * 3)))
        print(
            "'{}_prevalent_incident': TensorMap('{}', group='categorical_date', channel_map={{'no_{}':0, 'prevalent_{}':1, 'incident_{}':2}}, loss=weighted_crossentropy([1.0, {}, {}], '{}_prevalent_incident')),".format(
                d, d, d, d, d, factor_p, factor_i, d,
            ),
        )


def _print_disease_tensor_maps_time(phenos_folder) -> None:
    disease2tsv = get_disease2tsv(phenos_folder)
    disease_tm_str = "'{}_time': TensorMap('{}', group='diagnosis_time', channel_map={{'{}_time':0}}, loss='mse'),"
    for d in sorted(list(disease2tsv.keys())):
        print(disease_tm_str.format(d, d, d))


def _plot_mi_hospital_only(db, run_id, output_folder) -> None:
    conn = sqlite3.connect(db)
    sql_cursor = conn.cursor()
    q = "SELECT datething.value FROM phenotype datething where datething.fieldid=42000 and sample_id in  (select sample_id from phenotype where FieldID=42001 and value=0);"
    dates = []
    for data_row in sql_cursor.execute(q):
        dates.append(_str2date(data_row[0]))
    plt.figure(figsize=(12, 12))
    plt.xlabel('MI Date')
    plt.hist(dates, bins=60)
    plt.savefig(os.path.join(output_folder, run_id, 'mi_self_report_dates' + IMAGE_EXT))


def _str2date(d) -> datetime.date:
    parts = d.split(CONCAT_CHAR)
    if len(parts) < 2:
        return datetime.datetime.now().date()
    return datetime.date(int(parts[0]), int(parts[1]), int(parts[2]))


def _date_from_dicom(d) -> str:
    return d.AcquisitionDate[0:4] + CONCAT_CHAR + d.AcquisitionDate[4:6] + CONCAT_CHAR + d.AcquisitionDate[6:]


def _datetime_from_dicom(d) -> datetime.date:
    return _str2date(_date_from_dicom(d))


def _log_extreme_n(stats, n) -> None:
    logging.info('Min values:')
    i = 0
    ordered = sorted(stats.items(), key=operator.itemgetter(1))
    for k, v in ordered:
        logging.info('{} has: {}'.format(k, v))
        i += 1
        if i > n:
            break
    logging.info('Max values:')
    i = 0
    for k, v in ordered[-n:-1]:
        logging.info('{} has: {}'.format(k, v))
        i += 1
        if i > n:
            break
    logging.info('\n\n')


def _prune_sample(
    sample_id: int, min_sample_id: int, max_sample_id: int, mri_field_ids: List[int],
    xml_field_ids: List[int], zip_folder: str, xml_folder: str,
):
    """Return True if the sample ID is missing associated MRI, EKG, or GT data.  Or if the sample_id is below the given minimum."""

    if sample_id < min_sample_id:
        return True
    if sample_id > max_sample_id:
        return True
    if len(mri_field_ids) > 0 and not (_sample_has_dicom_mris(zip_folder, sample_id) or _sample_has_nifti_mris(zip_folder, sample_id)):
        return True
    if len(xml_field_ids) > 0 and not _sample_has_ecgs(xml_folder, xml_field_ids, sample_id):
        return True

    return False
