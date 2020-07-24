from typing import List
from enum import Enum, auto


class StorageType(Enum):
    CONTINUOUS = auto()
    CATEGORICAL_INDEX = auto()
    CATEGORICAL_FLAG = auto()
    ONE_HOT = auto()
    STRING = auto()
    BYTE_STRING = auto()

    def __str__(self):
        """StorageType.FLOAT_ARRAY becomes float_array"""
        return str.lower(super().__str__().split('.')[1])


DATE_FORMAT = '%Y-%m-%d_%H-%M-%S'
PARTNERS_DATE_FORMAT = '%m-%d-%Y'
PARTNERS_TIME_FORMAT = '%H:%M:%S'
PARTNERS_DATETIME_FORMAT = '%Y-%m-%dT%H:%M:%S'
CARDIAC_SURGERY_DATE_FORMAT = "%Y-%m-%d"

EPS = 1e-7

DICOM_EXT = '.dcm'
IMAGE_EXT = '.png'
PDF_EXT = '.pdf'
TENSOR_EXT = '.hd5'
MODEL_EXT = '.h5'
XML_EXT = '.xml'

STOP_CHAR = '!'
JOIN_CHAR = '_'
CONCAT_CHAR = '-'
HD5_GROUP_CHAR = '/'

MRI_FRAMES = 50
MRI_TO_SEGMENT = 'cine_segmented_sax_inlinevf'
MRI_LAX_TO_SEGMENT = 'cine_segmented_lax_inlinevf'
MRI_SEGMENTED = 'cine_segmented_sax_inlinevf_segmented'
MRI_LAX_SEGMENTED = 'cine_segmented_lax_inlinevf_segmented'
MRI_ZOOM_INPUT = 'sax_inlinevf_zoom'
MRI_ZOOM_MASK = 'sax_inlinevf_zoom_mask'
MRI_DATE = 'mridate'
MRI_ANNOTATION_NAME = 'mri_critic_annotation'
MRI_PIXEL_WIDTH = 'mri_pixel_width'
MRI_PIXEL_HEIGHT = 'mri_pixel_height'
MRI_SLICE_THICKNESS = 'mri_slice_thickness'
MRI_PATIENT_POSITION = 'mri_patient_position'
MRI_PATIENT_ORIENTATION = 'mri_patient_orientation'
MRI_SEGMENTED_CHANNEL_MAP = {'background': 0, 'ventricle': 1, 'myocardium': 2}
MRI_ANNOTATION_CHANNEL_MAP = {'good': 0, 'included-lvot': 1, 'mistrace': 2, 'phantom-apex': 3, 'hardclip': 4}
MRI_LAX_3CH_SEGMENTED_CHANNEL_MAP = {'background': 0, 'LV_anteroseptum': 1, 'left_atrium': 2, 'LV_inferior_wall': 3, 'LV_Papillary': 4, 'LV_Cavity': 5}
MRI_LAX_4CH_SEGMENTED_CHANNEL_MAP = {
    'background': 0, 'RV_free_wall': 1, 'RA_free_wall': 2, 'LA_free_wall': 3, 'LV_anterolateral_wall': 4,
    'interventricular_septum': 5, 'interatrial_septum': 6, 'crista_terminalis': 7, 'RA_cavity': 8, 'RV_cavity': 9,
    'LA_cavity': 10, 'LV_cavity': 11, 'descending_aorta': 12, 'thoracic_cavity': 13,
}
MRI_SAX_SEGMENTED_CHANNEL_MAP = {
    'background': 0, 'RV_free_wall': 1, 'interventricular_septum': 2, 'LV_free_wall': 3, 'LV_pap': 4, 'LV_cavity': 5,
    'RV_cavity': 6, 'thoracic_cavity': 7, 'liver': 8, 'stomach': 9, 'spleen': 10,
}
MRI_AO_SEGMENTED_CHANNEL_MAP = {
    'background': 0, 'svc': 1, 'pulmonary_artery': 2, 'ascending_aortic_wall': 3, 'ascending_aorta': 4,
    'descending_aortic_wall': 5, 'descending_aorta': 6, 'thorax': 7, 'body': 8, 'breast_implant': 9,
}
MRI_LIVER_SEGMENTED_CHANNEL_MAP = {'background': 0, 'liver': 1, 'inferior_vena_cava': 2, 'abdominal_aorta': 3, 'body': 4}

CAD_ICDS = [
    'K401', 'K402', 'K403', 'K404', 'K411', 'K412', 'K413', 'K414', 'K451', 'K452', 'K453', 'K454', 'K455',
    'K491', 'K492', 'K498', 'K499', 'K502', 'K751', 'K752', 'K753', 'K754', 'K758', 'K759',
]

# TODO: These values should ultimately come from the coding table
CODING_VALUES_LESS_THAN_ONE = [-10, -1001]
CODING_VALUES_MISSING = [-3, -1, -2, -11, -818, -121, -313, -906]

ECG_REST_LEADS = {
    'strip_I': 0, 'strip_II': 1, 'strip_III': 2, 'strip_V1': 3, 'strip_V2': 4, 'strip_V3': 5,
    'strip_V4': 6, 'strip_V5': 7, 'strip_V6': 8, 'strip_aVF': 9, 'strip_aVL': 10, 'strip_aVR': 11,
}
ECG_REST_MEDIAN_LEADS = {
    'median_I': 0, 'median_II': 1, 'median_III': 2, 'median_V1': 3, 'median_V2': 4, 'median_V3': 5,
    'median_V4': 6, 'median_V5': 7, 'median_V6': 8, 'median_aVF': 9, 'median_aVL': 10, 'median_aVR': 11,
}
ECG_REST_AMP_LEADS = {
    'I': 0, 'II': 1, 'III': 2, 'aVR': 3, 'aVL': 4, 'aVF': 5,
    'V1': 6, 'V2': 7, 'V3': 8, 'V4': 9, 'V5': 10, 'V6': 11,
}
ECG_SEGMENTED_CHANNEL_MAP = {'unknown': 0, 'TP_segment': 1, 'P_wave': 2, 'PQ_segment': 3, 'QRS_complex': 4, 'ST_segment': 5, 'T_wave': 6, 'U_wave': 7}

ECG_BIKE_LEADS = {"I": 0, "2": 1, "3": 2}
ECG_BIKE_MEDIAN_SIZE = (5500, len(ECG_BIKE_LEADS))
ECG_BIKE_STRIP_SIZE = (5000, len(ECG_BIKE_LEADS))
ECG_BIKE_FULL_SIZE = (216500, len(ECG_BIKE_LEADS))
ECG_BIKE_RECOVERY_SIZE = (500 * 60, len(ECG_BIKE_LEADS))

ECG_CHAR_2_IDX = {
    '!': 0, ' ': 1, "'": 2, ')': 4, '(': 3, '-': 5, '/': 6, '1': 7, '3': 9, '2': 8, '4': 10, ':': 11,
    '?': 12, 'A': 13, 'C': 15, 'B': 14, 'E': 17, 'D': 16, 'G': 18, 'I': 20, 'H': 19, 'J': 21, 'M': 23,
    'L': 22, 'N': 24, 'Q': 26, 'P': 25, 'S': 28, 'R': 27, 'U': 30, 'T': 29, 'W': 32, 'V': 31, 'a': 33,
    'c': 35, 'b': 34, 'e': 37, 'd': 36, 'g': 39, 'f': 38, 'i': 41, 'h': 40, 'k': 43, 'j': 42, 'm': 45,
    'l': 44, 'o': 47, 'n': 46, 'q': 49, 'p': 48, 's': 51, 'r': 50, 'u': 53, 't': 52, 'w': 55, 'v': 54,
    'y': 57, 'x': 56, 'z': 58, 'O': 59, '5': 60,
}
ECG_IDX_2_CHAR = {
    0: '!', 1: ' ', 2: "'", 3: '(', 4: ')', 5: '-', 6: '/', 7: '1', 8: '2', 9: '3', 10: '4', 11: ':',
    12: '?', 13: 'A', 14: 'B', 15: 'C', 16: 'D', 17: 'E', 18: 'G', 19: 'H', 20: 'I', 21: 'J', 22: 'L',
    23: 'M', 24: 'N', 25: 'P', 26: 'Q', 27: 'R', 28: 'S', 29: 'T', 30: 'U', 31: 'V', 32: 'W', 33: 'a',
    34: 'b', 35: 'c', 36: 'd', 37: 'e', 38: 'f', 39: 'g', 40: 'h', 41: 'i', 42: 'j', 43: 'k', 44: 'l',
    45: 'm', 46: 'n', 47: 'o', 48: 'p', 49: 'q', 50: 'r', 51: 's', 52: 't', 53: 'u', 54: 'v', 55: 'w',
    56: 'x', 57: 'y', 58: 'z', 59: 'O', 60: '5',
}

PARTNERS_READ_TEXT = 'read_'
PARTNERS_CHAR_2_IDX = {
    ' ': 0, '0': 1, '1': 2, '2': 3, '3': 4, '4': 5, '5': 6, '6': 7, '7': 8, '8': 9, '9': 10, 'a': 11, 'b': 12, 'c': 13, 'd': 14, 'e': 15, 'f': 16, 'g': 17,
    'h': 18, 'i': 19, 'j': 20, 'k': 21, 'l': 22, 'm': 23, 'n': 24, 'o': 25, 'p': 26, 'q': 27, 'r': 28, 's': 29, 't': 30, 'u': 31, 'v': 32, 'w': 33, 'x': 34,
    'y': 35, 'z': 36, '\\': 37, '.': 38, STOP_CHAR: 39,
}
PARTNERS_IDX_2_CHAR = {
    0: ' ', 1: '0', 2: '1', 3: '2', 4: '3', 5: '4', 6: '5', 7: '6', 8: '7', 9: '8', 10: '9', 11: 'a', 12: 'b', 13: 'c', 14: 'd', 15: 'e', 16: 'f', 17: 'g',
    18: 'h', 19: 'i', 20: 'j', 21: 'k', 22: 'l', 23: 'm', 24: 'n', 25: 'o', 26: 'p', 27: 'q', 28: 'r', 29: 's', 30: 't', 31: 'u', 32: 'v', 33: 'w', 34: 'x',
    35: 'y', 36: 'z', 37: '\\', 38: '.', 39: STOP_CHAR,
}

TENSOR_MAPS_FILE_NAME = 'tensor_maps_by_script'

# BigQuery tables
SQL_DATASET = "ukbb7089_201904"
DICTIONARY_TABLE = SQL_DATASET+".dictionary"
CODING_TABLE = SQL_DATASET+".coding"
PHENOTYPE_TABLE = SQL_DATASET+".phenotype"

GCS_BUCKET = 'ml4cvd'

TENSOR_MAP_GROUP_CONTINUOUS = 'multi_field_no_missing_channel_continuous'
TENSOR_MAP_GROUP_MISSING_CONTINUOUS = 'multi_field_continuous'
IMPUTATION_RANDOM = 'random'
IMPUTATION_MEAN = 'mean'


def dataset_name_from_meaning(group: str, fields: List[str]) -> str:
    clean_fields = []
    for f in fields:
        clean_fields.append(''.join(e for e in f if e.isalnum() or e == ' '))
    joined = JOIN_CHAR.join(clean_fields).replace('  ', CONCAT_CHAR).replace(' ', CONCAT_CHAR)
    if group is None:
        return joined
    return group + HD5_GROUP_CHAR + joined
