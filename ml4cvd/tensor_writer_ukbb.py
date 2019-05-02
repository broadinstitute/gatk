# tensor_writer_ukbb.py
#
# UK Biobank-specific tensor writing, SQL querying, data munging goes here
#

# Imports
import os
import re
import csv
import glob
import h5py
import shutil
import logging
import sqlite3
import zipfile
import pydicom
import datetime
import operator
import numpy as np
from typing import Dict, List, Tuple
from timeit import default_timer as timer
from collections import Counter, defaultdict

import matplotlib

matplotlib.use('Agg')  # Need this to write images from the GSA servers.  Order matters:
import matplotlib.pyplot as plt  # First import matplotlib, then use Agg, then import plt
from PIL import Image, ImageDraw  # Polygon to mask
import xml.etree.ElementTree as et
from scipy.ndimage.morphology import binary_closing  # Morphological operator

from plots import plot_value_counter, plot_histograms
from defines import IMAGE_EXT, TENSOR_EXT, DICOM_EXT, JOIN_CHAR, CONCAT_CHAR, HD5_GROUP_CHAR
from defines import ECG_BIKE_LEADS, ECG_BIKE_MEDIAN_SIZE, ECG_BIKE_STRIP_SIZE, ECG_BIKE_FULL_SIZE
from defines import MRI_DATE, MRI_FRAMES, MRI_SEGMENTED, MRI_TO_SEGMENT, MRI_ZOOM_INPUT, MRI_ZOOM_MASK

MRI_SERIES_TO_WRITE = ['cine_segmented_lax_2ch', 'cine_segmented_lax_3ch', 'cine_segmented_lax_4ch', 'cine_segmented_sax_b1', 'cine_segmented_sax_b2',
                       'cine_segmented_sax_b3', 'cine_segmented_sax_b4', 'cine_segmented_sax_b5', 'cine_segmented_sax_b6', 'cine_segmented_sax_b7',
                       'cine_segmented_sax_b8', 'cine_segmented_sax_b9', 'cine_segmented_sax_b10', 'cine_segmented_sax_b11',
                       'cine_segmented_sax_inlinevf']

ECG_BIKE_FIELD = '6025'
ECG_REST_FIELD = '20205'
ECG_SINUS = ['Normal_sinus_rhythm', 'Sinus_bradycardia', 'Marked_sinus_bradycardia', 'Atrial_fibrillation']
ECG_NORMALITY = ['Normal_ECG', 'Abnormal_ECG', 'Borderline_ECG', 'Otherwise_normal_ECG']
ECG_TAGS_TO_WRITE = ['VentricularRate', 'PQInterval', 'PDuration', 'QRSDuration', 'QTInterval', 'QTCInterval', 'RRInterval', 'PPInterval',
                     'SokolovLVHIndex', 'PAxis', 'RAxis', 'TAxis', 'QTDispersion', 'QTDispersionBazett', 'QRSNum', 'POnset', 'POffset', 'QOnset',
                     'QOffset', 'TOffset']


def write_tensors(a_id: str,
                  db: str,
                  xml_folder: str,
                  zip_folder: str,
                  phenos_folder: str,
                  output_folder: str,
                  tensors: str,
                  dicoms: str,
                  volume_csv: str,
                  lv_mass_csv: str,
                  icd_csv: str,
                  categorical_field_ids: List[int],
                  continuous_field_ids: List[int],
                  mri_field_ids: List[int],
                  xml_field_ids: List[int],
                  x: int,
                  y: int,
                  z: int,
                  include_heart_zoom: bool,
                  zoom_x: int,
                  zoom_y: int,
                  zoom_width: int,
                  zoom_height: int,
                  write_pngs: bool,
                  min_sample_id: int,
                  max_sample_id: int,
                  min_values_to_print: int) -> None:
    """Write tensors as HD5 files containing any kind of data from UK BioBank

    One HD5 file is generated per sample.  Each file may contain many tensor encodings of data including:
     survey responses, MRI, and ECG.

    :param a_id: User chosen string to identify this run
    :param db: Path to SQL database for UKBB
    :param xml_folder: Path to folder containing ECG XML files
    :param zip_folder: Path to folder containing zipped DICOM files
    :param phenos_folder: Path to folder containing phenotype CSVs
    :param output_folder: Folder to write outputs to (mostly for debugging)
    :param tensors: Folder to populate with HD5 tensors
    :param dicoms: Folder where zipped DICOM will be decompressed
    :param volume_csv: CSV containing systole and diastole volumes and ejection fraction for samples with MRI
    :param lv_mass_csv: CSV containing LV mass and other data from a subset of samples with MRI
    :param icd_csv: CSV of diagnoses
    :param categorical_field_ids: List of categorical field IDs from UKBB
    :param continuous_field_ids: List of continuous field IDs from UKBB
    :param mri_field_ids: List of MRI field IDs from UKBB
    :param xml_field_ids: List of ECG field IDs from UKBB
    :param x: Maximum x dimension of MRIs
    :param y: Maximum y dimension of MRIs
    :param z: Maximum z dimension of MRIs
    :param include_heart_zoom: If True, also extract zoomed MRI of the heart
    :param zoom_x: x coordinate of the zoom
    :param zoom_y: y coordinate of the zoom
    :param zoom_width: width of the zoom
    :param zoom_height: height of the zoom
    :param write_pngs: write MRIs as PNG images for debugging
    :param min_sample_id: Minimum sample id to generate, for parallelization
    :param max_sample_id: Maximum sample id to generate, for parallelization
    :param min_values_to_print: Minimum number of samples that have responded to question for it to be included in the
            categorical or continuous dictionaries printed after tensor generation
    :param ukbb7089_sample_id_to_hail_pkl_path: path to the pkl file containing all available genotypes.
    :param filtered_genotypes_array_path: path to the filtered genotypes.

    :return: None
    """
    stats = Counter()
    continuous_stats = defaultdict(list)
    conn = sqlite3.connect(db)
    sql_cursor = conn.cursor()

    field_meanings, dates, status, lvef, lvesv, lvedv, lv_mass, icds, sample_ids = _load_meta_data_for_tensor_writing(sql_cursor, phenos_folder,
                                                                                                                      volume_csv, lv_mass_csv,
                                                                                                                      icd_csv, categorical_field_ids,
                                                                                                                      continuous_field_ids,
                                                                                                                      min_sample_id, max_sample_id)

    for sample_id in sorted(sample_ids):

        start_time = timer()  # Keep track of elapsed execution time

        tensor_path = os.path.join(tensors, str(sample_id) + TENSOR_EXT)
        if not os.path.exists(os.path.dirname(tensor_path)):
            os.makedirs(os.path.dirname(tensor_path))
        if _prune_sample(sample_id, min_sample_id, max_sample_id, mri_field_ids, xml_field_ids, zip_folder, xml_folder):
            continue
        try:
            with h5py.File(tensor_path, 'w') as hd5:
                _write_tensors_from_zipped_dicoms(x, y, z, include_heart_zoom, zoom_x, zoom_y, zoom_width, zoom_height, write_pngs, tensors, dicoms,
                                                  mri_field_ids, zip_folder, hd5, sample_id, stats)
                _write_tensors_from_xml(xml_field_ids, xml_folder, hd5, sample_id, write_pngs, stats, continuous_stats)
                _write_tensors_from_sql(sql_cursor, hd5, sample_id, categorical_field_ids, continuous_field_ids, field_meanings, write_pngs,
                                        continuous_stats, stats)
                _write_tensors_from_icds(hd5, sample_id, icds, dates, stats)
                _write_tensors_from_volumes(hd5, sample_id, lv_mass, lvef, lvesv, lvedv, continuous_stats)
                stats['Tensors written'] += 1
        except AttributeError:
            logging.exception('Encountered AttributeError trying to write a UKBB tensor at path:{}'.format(tensor_path))
            logging.info('Deleting attempted tensor at path:{}'.format(tensor_path))
            os.remove(tensor_path)
        except ValueError:
            logging.exception('Encountered ValueError trying to write a UKBB tensor at path:{}'.format(tensor_path))
            logging.info('Deleting attempted tensor at path:{}'.format(tensor_path))
            os.remove(tensor_path)
        except RuntimeError:
            logging.exception('Encountered RuntimeError trying to write a UKBB tensor at path:{}'.format(tensor_path))
            logging.info('Deleting attempted tensor at path:{}'.format(tensor_path))
            os.remove(tensor_path)
        except OSError:
            logging.exception('Encountered OSError trying to write a UKBB tensor at path:{}'.format(tensor_path))
            logging.info('Deleting attempted tensor at path:{}'.format(tensor_path))
            os.remove(tensor_path)

        end_time = timer()
        elapsed_time = end_time - start_time
        logging.info("Populated {} in {} seconds.".format(tensor_path, elapsed_time))

    _dicts_and_plots_from_tensorization(a_id, output_folder, min_values_to_print, write_pngs, continuous_stats, stats)


def _load_meta_data_for_tensor_writing(sql_cursor: sqlite3.Cursor,
                                       phenos_folder: str,
                                       volume_csv: str,
                                       lv_mass_csv: str,
                                       icd_csv: str,
                                       categorical_field_ids: List[int],
                                       continuous_field_ids: List[int],
                                       min_sample_id: int,
                                       max_sample_id: int) -> Tuple[Dict[str, str],
                                                                    Dict[str, Dict[int, List[datetime.date]]],
                                                                    Dict[str, Dict[int, List[int]]],
                                                                    Dict[int, float],
                                                                    Dict[int, float],
                                                                    Dict[int, float],
                                                                    Dict[str, Dict[int, List[int]]],
                                                                    List[int]]:
    """ Gather metadata necessary to write tensors from UK biobank

    Loads the field IDs of survey data, dates of assessment, diagnosis of diseases,
    ejection fractions, diastolic volumes, systolic volumes, and sample IDs to make tensors from.

    :param sql_cursor: SQL Cursor object used to query meanings of the field IDs
    :param phenos_folder: Path to folder containing phenotype CSVs
    :param volume_csv: CSV containing systole and diastole volumes and ejection fraction for samples with MRI
    :param lv_mass_csv: TSV containing left ventricular mass and other cardiac MRI readouts on ~5000 people returned from app 2964
    :param icd_csv: CSV of diagnoses
    :param categorical_field_ids: List of categorical field IDs from UKBB
    :param continuous_field_ids: List of continuous field IDs from UKBB
    :param min_sample_id: Minimum sample id to generate, for parallelization
    :param max_sample_id: Maximum sample id to generate, for parallelization
    :return: Tuple of metadata
        field_meanings: Dict mapping field IDs (as strings) to their meaning as a string
        dates: Nested Dictionary mapping each disease string to a dictionary
                mapping sample ids to the date of their diagnosis or last update
        status: Nested Dictionary mapping each disease string to a dictionary
                mapping sample ids to an int (1 if they have the disease, 0 if not)
        lvef: Dictionary mapping sample IDs (as ints) to left ventricle ejection fraction as a float
        lvesv: Dictionary mapping sample IDs (as ints) to left ventricle end systole volume as a float
        lvedv: Dictionary mapping sample IDs (as ints) to left ventricle end diastole volume as a float
        icds: Nested Dictionary mapping each disease string to a dictionary
                mapping sample ids to an int (1 if they have the disease, 0 if not)
        sample_ids: List of sample IDs (as ints) to generate tensors for
    """
    field_meanings = {}
    for fid in categorical_field_ids + continuous_field_ids:
        for field_row in sql_cursor.execute("SELECT field FROM dictionary WHERE fieldid=%d" % (fid)):
            field_meanings[fid] = field_row[0]

    disease2tsv = get_disease2tsv(phenos_folder)
    logging.info('Got disease TSVs:{}'.format(disease2tsv))
    dates = _get_disease_censor_dates(disease2tsv, min_sample_id, max_sample_id)
    logging.info('Got dates for diseases:{}'.format(list(dates.keys())))
    status = disease_censor_status(disease2tsv, min_sample_id, max_sample_id)
    logging.info('Got statuses and dates for diseases:{}'.format(list(status.keys())))

    lvef = {}
    lvesv = {}
    lvedv = {}
    with open(volume_csv, 'r') as volumes:
        lol = list(csv.reader(volumes, delimiter='\t'))
        logging.info('CSV of MRI volumes header:{}'.format(list(enumerate(lol[0]))))
        for row in lol[1:]:
            sample_id = int(row[0])
            if min_sample_id <= sample_id <= max_sample_id:
                lvesv[sample_id] = float(row[2])
                lvedv[sample_id] = float(row[4])
                lvef[sample_id] = float(row[6])

    lv_mass = {}
    with open(lv_mass_csv, 'r') as lvm:
        lol = list(csv.reader(lvm, delimiter='\t'))
        logging.info('CSV of returned MRI mass, etc, header:{}'.format(list(enumerate(lol[0]))))
        for row in lol[1:]:
            # column 0 is the original app's sample ID. Column 1 is app 7089's sample ID.
            sample_id = int(row[1])
            if min_sample_id <= sample_id <= max_sample_id and len(row[13]) > 0:
                # Zero-based column #13 is the LV mass
                lv_mass[sample_id] = float(row[13])

    sample_ids = []
    icds = defaultdict(dict)
    lol = list(csv.reader(open(icd_csv, 'r'), delimiter='\t'))
    logging.info('CSV of ICDs header:{}'.format(list(enumerate(lol[0]))))
    for row in lol[1:]:
        sample_id = int(row[0])
        if min_sample_id <= sample_id <= max_sample_id:
            sample_ids.append(sample_id)
            for disease in disease2tsv:
                if sample_id not in status[disease]:
                    continue
                icds[disease][sample_id] = status[disease][sample_id]

    return field_meanings, dates, status, lvef, lvesv, lvedv, lv_mass, icds, sample_ids


def _sample_has_mris(zip_folder, sample_id) -> bool:
    sample_str = str(sample_id)
    return os.path.exists(zip_folder + sample_str + '_20209_2_0.zip') or os.path.exists(
        zip_folder + sample_str + '_20209_0_0.zip') or os.path.exists(
        zip_folder + sample_str + '_20209_1_0.zip')


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


def _dicts_and_plots_from_tensorization(a_id: str,
                                        output_folder: str,
                                        min_values_to_print: int,
                                        write_pngs: bool,
                                        continuous_stats: Dict[str, List[float]],
                                        stats: Dict[str, int]) -> None:
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


def _write_tensors_from_sql(sql_cursor: sqlite3.Cursor,
                            hd5: h5py.File,
                            sample_id: int,
                            categorical_field_ids: List[int],
                            continuous_field_ids: List[int],
                            field_meanings: Dict[str, str],
                            write_pngs: bool,
                            continuous_stats: Dict[str, List[float]],
                            stats: Dict[str, int]) -> None:
    """Tensorize and write UK biobank survey data from continuous and categorical fields
    into the HD5 file of a given sample

    :param sql_cursor: SQL Cursor object used to query meanings of the field IDs
    :param hd5: HD5 File where all tensors for this sample ID are saved
    :param sample_id: The sample ID whose tensor we are currently writing
    :param categorical_field_ids: List of categorical field IDs from UKBB
    :param continuous_field_ids: List of continuous field IDs from UKBB
    :param field_meanings: Dict mapping field IDs (as strings) to their meaning as a string
    :param write_pngs: write MRIs as PNG images for debugging
    :param continuous_stats: Dictionary mapping field meanings to the list of continuous values found for them
    :param stats: Dictionary mapping strings to ints keeps track of categorical responses and other info.
    :return: None
    """
    data_query = "SELECT m.meaning, p.instance, p.array_idx, p.value FROM phenotype p JOIN Coding m"
    data_query += " ON m.coding=p.value AND m.coding_file_id=p.coding_file_id"
    data_query += " WHERE p.fieldid=%d and p.sample_id=%d;"
    for fid in categorical_field_ids:
        try:
            for data_row in sql_cursor.execute(data_query % (fid, sample_id)):
                dataset_name = _dataset_name_from_meaning('categorical',
                                                          [field_meanings[fid], str(data_row[0]),
                                                           str(data_row[1]), str(data_row[2])])
                float_category = _to_float_or_false(data_row[3])
                if float_category is not False:
                    hd5.create_dataset(dataset_name, data=[float_category])
                    stats[dataset_name] += 1
                else:
                    logging.warning('Cannot cast float from:{} categorical field:{} means:{} sample id:{}'.format(
                        data_row[3], fid, data_row[0], sample_id))
        except:
            logging.exception(
                'problem with data query field id:{} data_query is:{} for sample id:{}'.format(
                    fid, data_query, sample_id))

    continuous_query = "SELECT p.instance, p.array_idx, p.value FROM phenotype p "
    continuous_query += "WHERE p.fieldid=%d and p.sample_id=%d;"
    for fid in continuous_field_ids:
        for data_row in sql_cursor.execute(continuous_query % (fid, sample_id)):
            dataset_name = _dataset_name_from_meaning('continuous', [str(fid), field_meanings[fid],
                                                                     str(data_row[0]), str(data_row[1])])
            float_continuous = _to_float_or_false(data_row[2])
            if float_continuous is not False:
                hd5.create_dataset(dataset_name, data=[float_continuous])
                stats[dataset_name] += 1
                if write_pngs:
                    continuous_stats[field_meanings[fid]].append(float_continuous)
            else:
                logging.warning('Cannot cast float:{} continuous:{} sample:{}'.format(data_row[2], fid, sample_id))

    date_query = "SELECT p.instance, p.array_idx, p.value FROM phenotype p "
    date_query += "WHERE p.fieldid=53 and p.sample_id=%d;"
    for data_row in sql_cursor.execute(date_query % sample_id):
        hd5.create_dataset('assessment-date_{}_{}'.format(data_row[0], data_row[1]), (1,), data=data_row[2],
                           dtype=h5py.special_dtype(vlen=str))
    _job_title_from_sql_to_tensor(sql_cursor, hd5, sample_id, stats)
    _icd10_from_sql_to_tensor(sql_cursor, hd5, sample_id, stats)


def _job_title_from_sql_to_tensor(sql_cursor: sqlite3.Cursor,
                                  hd5: h5py.File,
                                  sample_id: int,
                                  stats: Dict[str, int]) -> None:
    """If found, tensorize the free text job titles as entered by user, removing duplicates

    :param sql_cursor: SQL Cursor object used to query meanings of the field IDs
    :param hd5: HD5 File where all tensors for this sample ID are saved
    :param sample_id: The sample ID whose tensor we are currently writing
    :param stats: Dictionary mapping strings to ints keeps track of categorical responses and other info.
    :return: None
    """
    jobs = []
    job_text_query = "SELECT instance, array_idx, value FROM phenotype WHERE fieldid=22600 AND sample_id=%d;"
    for data_row in sql_cursor.execute(job_text_query % sample_id):
        jobs.append(data_row[2])
        stats['Job_' + data_row[2]] += 1
    if len(jobs) > 0:
        jobs = list(set(jobs))
        hd5.create_dataset('jobs', (1,), data=JOIN_CHAR.join(jobs), dtype=h5py.special_dtype(vlen=str))


def _icd10_from_sql_to_tensor(sql_cursor: sqlite3.Cursor,
                              hd5: h5py.File,
                              sample_id: int,
                              stats: Dict[str, int]) -> None:
    """If found, tensorize the ICD10 codes

    :param sql_cursor: SQL Cursor object used to query meanings of the field IDs
    :param hd5: HD5 File where all tensors for this sample ID are saved
    :param sample_id: The sample ID whose tensor we are currently writing
    :param stats: Dictionary mapping strings to ints keeps track of categorical responses and other info.
    :return: None
    """
    icds = []
    icd10_fields = [41202, 41204, 40001, 40002, 40006]
    icd_query = "SELECT value FROM phenotype WHERE fieldid=%d AND sample_id=%d;"
    for icd10_field in icd10_fields:
        for data_row in sql_cursor.execute(icd_query % (icd10_field, sample_id)):
            icds.append(data_row[0])
            # stats['ICD10_' + data_row[0]] += 1
    if len(icds) > 0:
        icds = list(set(icds))
        hd5.create_dataset('icd10', (1,), data=JOIN_CHAR.join(icds), dtype=h5py.special_dtype(vlen=str))


def _write_tensors_from_volumes(hd5: h5py.File,
                                sample_id: int,
                                lv_mass: Dict[int, float],
                                lvef: Dict[int, float],
                                lvesv: Dict[int, float],
                                lvedv: Dict[int, float],
                                continuous_stats: Dict[str, List[float]]) -> None:
    """Write volumes of the left ventricle into tensor file of a particular sample if we have them

    :param hd5: HD5 File where all tensors for this sample ID are saved
    :param sample_id: The sample ID whose tensor we are currently writing
    :param lv_mass: Dictionary mapping sample IDs (as ints) to left ventricular mass as a float
    :param lvef: Dictionary mapping sample IDs (as ints) to left ventricle ejection fraction as a float
    :param lvesv: Dictionary mapping sample IDs (as ints) to left ventricle end systole volume as a float
    :param lvedv: Dictionary mapping sample IDs (as ints) to left ventricle end diastole volume as a float
    :param continuous_stats: Dictionary mapping field meanings to the list of continuous values found for them
    :return: None
    """
    if sample_id in lvef and sample_id in lvesv and sample_id in lvedv:
        hd5.create_dataset('continuous' + HD5_GROUP_CHAR + 'end_systole_volume', data=[lvesv[sample_id]])
        continuous_stats['end_systole_volume'].append(lvesv[sample_id])
        hd5.create_dataset('continuous' + HD5_GROUP_CHAR + 'end_diastole_volume', data=[lvedv[sample_id]])
        continuous_stats['end_diastole_volume'].append(lvedv[sample_id])
        hd5.create_dataset('continuous' + HD5_GROUP_CHAR + 'ejection_fraction', data=[lvef[sample_id]])
        continuous_stats['ejection_fraction'].append(lvef[sample_id])

    if sample_id in lv_mass:
        hd5.create_dataset('lv_mass', data=[lv_mass[sample_id]])
        continuous_stats['lv_mass'].append(lv_mass[sample_id])


def _write_tensors_from_icds(hd5: h5py.File,
                             sample_id: int,
                             icds: Dict[str, Dict[int, List[int]]],
                             dates: Dict[str, Dict[int, List[datetime.date]]],
                             stats: Dict[str, int]) -> None:
    """Write disease status and dates into the HD5 file of a given sample

    :param hd5: HD5 File where all tensors for this sample ID are saved
    :param sample_id: The sample ID whose tensor we are currently writing
    :param dates: Nested Dictionary mapping each disease string to a dictionary
                mapping sample ids to the date of their diagnosis or last update
    :param icds: Nested Dictionary mapping each disease string to a dictionary
                mapping sample ids to an int (1 if they have the disease, 0 if not)
    :param stats: Dictionary mapping strings to ints keeps track of categorical responses and other info.
    :return: None
    """
    for disease in icds:
        if sample_id in icds[disease] and icds[disease][sample_id] != 'NA':
            hd5.create_dataset(disease, data=[float(icds[disease][sample_id])])
            disease_str = disease.replace(JOIN_CHAR, CONCAT_CHAR)
            stats['categorical/' + disease_str + '-' + str(icds[disease][sample_id])] += 1
            if sample_id in dates[disease]:
                disease_date = dates[disease][sample_id].strftime('%Y-%m-%d')
                hd5.create_dataset(disease + '_date', (1,), data=disease_date, dtype=h5py.special_dtype(vlen=str))


def _dataset_name_from_meaning(group: str, fields: List[str]) -> str:
    clean_fields = []
    for f in fields:
        clean_fields.append(''.join(e for e in f if e.isalnum() or e == ' '))
    joined = JOIN_CHAR.join(clean_fields).replace('  ', CONCAT_CHAR).replace(' ', CONCAT_CHAR)
    return group + HD5_GROUP_CHAR + joined


def _to_float_or_false(s):
    try:
        return float(s)
    except ValueError:
        return False


def _write_tensors_from_zipped_dicoms(x,
                                      y,
                                      z,
                                      include_heart_zoom,
                                      zoom_x,
                                      zoom_y,
                                      zoom_width,
                                      zoom_height,
                                      write_pngs,
                                      tensors,
                                      dicoms,
                                      mri_field_ids,
                                      zip_folder,
                                      hd5,
                                      sample_id,
                                      stats) -> None:
    sample_str = str(sample_id)
    for mri_field in mri_field_ids:
        mris = glob.glob(zip_folder + sample_str + '_' + mri_field + '*.zip')
        for zipped in mris:
            logging.info("Got zipped dicoms for sample: {} with MRI field: {}".format(sample_id, mri_field))
            dicom_folder = os.path.join(dicoms, sample_str, mri_field)
            if not os.path.exists(dicom_folder):
                os.makedirs(dicom_folder)
            with zipfile.ZipFile(zipped, "r") as zip_ref:
                zip_ref.extractall(dicom_folder)
                _write_tensors_from_dicoms(x, y, z, include_heart_zoom, zoom_x, zoom_y, zoom_width, zoom_height,
                                           write_pngs, tensors, dicom_folder, hd5, stats)
                stats['MRI fields written'] += 1
            shutil.rmtree(dicom_folder)


def _write_tensors_from_dicoms(x,
                               y,
                               z,
                               include_heart_zoom,
                               zoom_x,
                               zoom_y,
                               zoom_width,
                               zoom_height,
                               write_pngs,
                               tensors,
                               dicom_folder,
                               hd5,
                               stats) -> None:
    """Convert a folder of DICOMs from a sample into tensors for each series

    Segmented dicoms require special processing and are written to tensor per-slice

    Arguments
        x: Width of the tensors (actual MRI width will be padded with 0s or cropped to this number)
        y: Height of the tensors (actual MRI width will be padded with 0s or cropped to this number)
        dicom_folder: Folder with all dicoms associated with one sample.
        hd5: Tensor file in which to create datasets for each series and each segmented slice
        stats: Counter to keep track of summary statistics

    """
    views = defaultdict(list)
    for dicom in os.listdir(dicom_folder):
        if os.path.splitext(dicom)[-1] != DICOM_EXT:
            continue
        d = pydicom.read_file(os.path.join(dicom_folder, dicom))
        if d.SeriesDescription.lower() in MRI_SERIES_TO_WRITE:
            views[d.SeriesDescription.lower()].append(d)
            stats[d.SeriesDescription.lower()] += 1

    diastoles = {}
    diastoles_pix = {}
    systoles = {}
    systoles_pix = {}

    for v in views:
        mri_shape = (views[v][0].Rows, views[v][0].Columns, len(views[v]))
        stats['mri shape:' + str(mri_shape)] += 1

        if v != MRI_TO_SEGMENT:
            mri_data = np.zeros((x, y, max(z, len(views[v]))), dtype=np.float32)
        else:
            full_slice = np.zeros((x, y), dtype=np.float32)
            full_mask = np.zeros((x, y), dtype=np.float32)

        for slicer in views[v]:
            if slicer.pixel_array.shape[0] == mri_shape[0] and slicer.pixel_array.shape[1] == mri_shape[1]:
                sx = min(slicer.Rows, x)
                sy = min(slicer.Columns, y)
                slice_data = slicer.pixel_array.astype(np.float32)[:sx, :sy]
                if v != MRI_TO_SEGMENT:
                    mri_data[:sx, :sy, slicer.InstanceNumber - 1] = slice_data
                elif v == MRI_TO_SEGMENT and _has_overlay(slicer):
                    overlay, mask = _get_overlay_from_dicom(slicer)
                    ventricle_pixels = np.count_nonzero(mask == 1)
                    if ventricle_pixels == 0:
                        continue
                    cur_angle = (slicer.InstanceNumber - 1) // MRI_FRAMES  # dicom InstanceNumber is 1-based
                    if not cur_angle in diastoles:
                        diastoles[cur_angle] = slicer
                        diastoles_pix[cur_angle] = ventricle_pixels
                        systoles[cur_angle] = slicer
                        systoles_pix[cur_angle] = ventricle_pixels
                    else:
                        if ventricle_pixels > diastoles_pix[cur_angle]:
                            diastoles[cur_angle] = slicer
                            diastoles_pix[cur_angle] = ventricle_pixels
                        if ventricle_pixels < systoles_pix[cur_angle]:
                            systoles[cur_angle] = slicer
                            systoles_pix[cur_angle] = ventricle_pixels

                    full_slice[:sx, :sy] = slice_data
                    full_mask[:sx, :sy] = mask
                    hd5.create_dataset(MRI_TO_SEGMENT + HD5_GROUP_CHAR + str(slicer.InstanceNumber), data=full_slice,
                                       compression='gzip')
                    hd5.create_dataset(MRI_SEGMENTED + HD5_GROUP_CHAR + str(slicer.InstanceNumber), data=full_mask,
                                       compression='gzip')
                    if MRI_DATE not in hd5:
                        hd5.create_dataset(MRI_DATE, (1,), data=_date_from_dicom(slicer),
                                           dtype=h5py.special_dtype(vlen=str))
                    if include_heart_zoom:
                        zoom_slice = full_slice[zoom_x: zoom_x + zoom_width, zoom_y: zoom_y + zoom_height]
                        zoom_mask = full_mask[zoom_x: zoom_x + zoom_width, zoom_y: zoom_y + zoom_height]
                        hd5.create_dataset(MRI_ZOOM_INPUT + HD5_GROUP_CHAR + str(slicer.InstanceNumber),
                                           data=zoom_slice, compression='gzip')
                        hd5.create_dataset(MRI_ZOOM_MASK + HD5_GROUP_CHAR + str(slicer.InstanceNumber), data=zoom_mask,
                                           compression='gzip')
                    if write_pngs:
                        overlay = np.ma.masked_where(overlay != 0, slicer.pixel_array)
                        # Note that plt.imsave renders the first dimension (our x) as vertical and our y as horizontal
                        plt.imsave(tensors + v + '_{0:3d}'.format(slicer.InstanceNumber) + '_mask' + IMAGE_EXT, mask)
                        plt.imsave(tensors + v + '_{0:3d}'.format(slicer.InstanceNumber) + '_overlay' + IMAGE_EXT, overlay)
                        if include_heart_zoom:
                            plt.imsave(tensors + v + '_{}'.format(slicer.InstanceNumber) + '_zslice' + IMAGE_EXT, zoom_slice)
                            plt.imsave(tensors + v + '_{}'.format(slicer.InstanceNumber) + '_zmask' + IMAGE_EXT, zoom_mask)
                if write_pngs:
                    plt.imsave(tensors + v + '_' + str(slicer.InstanceNumber) + IMAGE_EXT, slicer.pixel_array)

        if v == MRI_TO_SEGMENT:
            # if len(diastoles) == 0:
            #     raise ValueError('Error not writing tensors if diastole could not be found.')
            for angle in diastoles:
                sx = min(diastoles[angle].Rows, x)
                sy = min(diastoles[angle].Columns, y)
                full_slice[:sx, :sy] = diastoles[angle].pixel_array.astype(np.float32)[:sx, :sy]
                overlay, full_mask[:sx, :sy] = _get_overlay_from_dicom(diastoles[angle])
                hd5.create_dataset('diastole_frame_b' + str(angle), data=full_slice, compression='gzip')
                hd5.create_dataset('diastole_mask_b' + str(angle), data=full_mask, compression='gzip')
                if write_pngs:
                    plt.imsave(tensors + 'diastole_frame_b' + str(angle) + IMAGE_EXT, full_slice)
                    plt.imsave(tensors + 'diastole_mask_b' + str(angle) + IMAGE_EXT, full_mask)

                sx = min(systoles[angle].Rows, x)
                sy = min(systoles[angle].Columns, y)
                full_slice[:sx, :sy] = systoles[angle].pixel_array.astype(np.float32)[:sx, :sy]
                overlay, full_mask[:sx, :sy] = _get_overlay_from_dicom(systoles[angle])
                hd5.create_dataset('systole_frame_b' + str(angle), data=full_slice, compression='gzip')
                hd5.create_dataset('systole_mask_b' + str(angle), data=full_mask, compression='gzip')
                if write_pngs:
                    plt.imsave(tensors + 'systole_frame_b' + str(angle) + IMAGE_EXT, full_slice)
                    plt.imsave(tensors + 'systole_mask_b' + str(angle) + IMAGE_EXT, full_mask)
        else:
            hd5.create_dataset(v, data=mri_data, compression='gzip')


def _has_overlay(d) -> bool:
    try:
        _ = d[0x6000, 0x3000].value
        return True
    except KeyError:
        return False


def _get_overlay_from_dicom(d) -> Tuple[np.ndarray, np.ndarray]:
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
        arr = np.ndarray(shape=(length_of_pixel_array * 8), dtype=np_dtype)
        for byte in overlay_raw:
            for bit in range(bit, bit + 8):
                arr[bit] = byte & 0b1
                byte >>= 1
            bit += 1
        arr = arr[:expected_bit_length]
    if overlay_frames == 1:
        arr = arr.reshape(rows, cols)
        idx = np.where(arr == 1)
        min_pos = (np.min(idx[0]), np.min(idx[1]))
        max_pos = (np.max(idx[0]), np.max(idx[1]))
        short_side = min((max_pos[0] - min_pos[0]), (max_pos[1] - min_pos[1]))
        small_radius = max(2, short_side * 0.185)
        big_radius = max(2, short_side * 0.9)
        myocardium_structure = _unit_disk(small_radius)
        m1 = binary_closing(arr, myocardium_structure).astype(np.int)
        ventricle_structure = _unit_disk(big_radius)
        m2 = binary_closing(arr, ventricle_structure).astype(np.int)
        return arr, m1 + m2


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
    rest_group = 'ecg_rest' + HD5_GROUP_CHAR
    categorical_group = 'categorical' + HD5_GROUP_CHAR
    for ecg in ecgs:
        logging.info('Got ECG for sample:{} XML field:{}'.format(sample_id, xml_field))
        root = et.parse(ecg).getroot()
        hd5.create_dataset('ecg_rest_date', (1,), data=_date_str_from_ecg(root), dtype=h5py.special_dtype(vlen=str))

        diagnosis_text = []
        found_infarct = False
        for d in root.findall("./Interpretation/Diagnosis/DiagnosisText"):
            if 'QRS Complexes:' in d.text:
                qrs = float(d.text.replace('QRS Complexes:', '').split(',')[0].strip())
                hd5.create_dataset('continuous' + HD5_GROUP_CHAR + 'QRSComplexes', data=[qrs])
            elif '---' in d.text or 'Arrhythmia results of the full-disclosure ECG' in d.text:
                continue
            elif d.text.replace(' ', '_') in ECG_NORMALITY:
                hd5.create_dataset(categorical_group + d.text.replace(' ', '_'), data=[1])
                stats[d.text] += 1
            else:
                if not found_infarct and 'infarct' in d.text:
                    hd5.create_dataset(categorical_group + 'infarct', data=[1])
                    found_infarct = True
                    stats['infarct'] += 1
                for sinus in ECG_SINUS:
                    if sinus in d.text.replace(' ', '_'):
                        hd5.create_dataset(categorical_group + sinus, data=[1])
                        stats[sinus] += 1
                diagnosis_text.append(d.text.replace(',', '').replace('*', '').replace('&', 'and').replace('  ', ' '))
        if not found_infarct:
            hd5.create_dataset(categorical_group + 'no_infarct', data=[1])
        diagnosis_str = ' '.join(diagnosis_text)
        hd5.create_dataset('ecg_rest_text', (1,), data=diagnosis_str, dtype=h5py.special_dtype(vlen=str))

        for c in root.findall("./StripData/WaveformData"):
            lead_data = list(map(float, c.text.strip().split(',')))
            dataset_name = 'strip_' + str(c.attrib['lead'])
            hd5.create_dataset(rest_group + dataset_name, data=lead_data, compression='gzip')
            stats[dataset_name] += 1

        for c in root.findall("./RestingECGMeasurements"):
            for child in c:
                if child.text is not None and child.tag in ECG_TAGS_TO_WRITE:
                    hd5.create_dataset('continuous' + HD5_GROUP_CHAR + child.tag, data=[float(child.text)])
                    stats[child.tag] += 1
                    if write_pngs:
                        continuous_stats[child.tag].append(float(child.text))
                if child.tag == 'MedianSamples':
                    for median_c in child:
                        if median_c.tag == 'WaveformData':
                            median_wave = list(map(float, median_c.text.strip().split(',')))
                            dataset_name = 'median_' + str(median_c.attrib['lead'])
                            hd5.create_dataset(rest_group + dataset_name, data=median_wave, compression='gzip')


def _write_ecg_bike_tensors(ecgs, xml_field, hd5, sample_id, stats):
    for ecg in ecgs:
        instance = ecg.split(JOIN_CHAR)[-2]
        logging.info('Got ECG for sample:{} XML field:{}'.format(sample_id, xml_field))
        root = et.parse(ecg).getroot()
        hd5.create_dataset('ecg_bike_date_{}'.format(instance), (1,),
                           data=_date_str_from_ecg(root), dtype=h5py.special_dtype(vlen=str))

        hd5_group = 'ecg_bike' + HD5_GROUP_CHAR
        median_ecgs = defaultdict(list)
        for median_waves in root.findall('./MedianData/Median/WaveformData'):
            median_ecgs[median_waves.attrib['lead']].extend(list(map(float, median_waves.text.strip().split(','))))
        if len(median_ecgs) > 0:
            median_np = np.zeros(ECG_BIKE_MEDIAN_SIZE)
            for lead in median_ecgs:
                median_idx = min(ECG_BIKE_MEDIAN_SIZE[0], len(median_ecgs[lead]))
                median_np[:median_idx, ECG_BIKE_LEADS[lead]] = median_ecgs[lead][:median_idx]
            hd5.create_dataset(hd5_group + 'median_{}'.format(instance), data=median_np, compression='gzip')
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
            hd5.create_dataset(hd5_group + 'strip_{}'.format(instance), data=strip_np, compression='gzip')
        else:
            stats['missing strip bike ECG'] += 1

        counter = 0
        full_ekgs = defaultdict(list)
        for lead_order in root.findall("./FullDisclosure/LeadOrder"):
            full_leads = {i: lead for i, lead in enumerate(lead_order.text.split(','))}
        for full_d in root.findall("./FullDisclosure/FullDisclosureData"):
            for full_line in re.split('\n|\t', full_d.text):
                for sample in re.split(',', full_line):
                    if sample == '':
                        continue
                    full_ekgs[full_leads[counter % 3]].append(float(sample))
                    counter += 1

        if len(full_ekgs) != 0:
            full_np = np.zeros(ECG_BIKE_FULL_SIZE)
            for lead in full_ekgs:
                full_idx = min(ECG_BIKE_FULL_SIZE[0], len(full_ekgs[lead]))
                full_np[:full_idx, ECG_BIKE_LEADS[lead]] = full_ekgs[lead][:full_idx]
            hd5.create_dataset(hd5_group + 'full_{}'.format(instance), data=full_np, compression='gzip')
        else:
            stats['missing full disclosure bike ECG'] += 1


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
                    stats[disease + '_icd_disparity_' + icds[disease][sample_id] + '_censor_' + str(
                        status[disease][sample_id])] += 1
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
        print("'{}': TensorMap('{}', group='categorical_index', channel_map={{'no_{}':0, '{}':1}}, loss=weighted_crossentropy([1.0, {}], '{}')),".format(
            d, d, d, d, factor, d))



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
                d, d, d, d, d, factor_p, factor_i, d))


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
    parts = d.split('-')
    if len(parts) < 2:
        return datetime.datetime.now().date()
    return datetime.date(int(parts[0]), int(parts[1]), int(parts[2]))


def _date_from_dicom(d) -> str:
    return d.AcquisitionDate[0:4] + CONCAT_CHAR + d.AcquisitionDate[4:6] + CONCAT_CHAR + d.AcquisitionDate[6:]


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


def _prune_sample(sample_id: int, min_sample_id: int, max_sample_id: int, mri_field_ids: List[int],
                  xml_field_ids: List[int], zip_folder: str, xml_folder: str):
    """Return True if the sample ID is missing associated MRI, EKG, or GT data.  Or if the sample_id is below the given minimum."""

    if sample_id < min_sample_id:
        return True
    if sample_id > max_sample_id:
        return True
    if len(mri_field_ids) > 0 and not _sample_has_mris(zip_folder, sample_id):
        return True
    if len(xml_field_ids) > 0 and not _sample_has_ecgs(xml_folder, xml_field_ids, sample_id):
        return True

    return False


if __name__ == '__main__':
    from arguments import parse_args

    args = parse_args()
    _print_disease_tensor_maps(args.phenos_folder)
    _print_disease_tensor_maps_time(args.phenos_folder)
    print('\n\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n\n\n')
    _print_disease_tensor_maps_incident_prevalent(args.phenos_folder)
    _plot_mi_hospital_only(args.db, args.id, args.output_folder)
    _ukbb_stats(args.id, args.output_folder, args.phenos_folder, args.volume_csv,
                args.icd_csv, args.app_csv, args.zip_folder)
