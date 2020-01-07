import logging
import datetime
from typing import Any
from dateutil import relativedelta

import numpy as np
from scipy.ndimage import zoom
from keras.utils import to_categorical

from ml4cvd.metrics import sentinel_logcosh_loss, survival_likelihood_loss, pearson
from ml4cvd.metrics import per_class_recall, per_class_recall_3d, per_class_recall_4d, per_class_recall_5d
from ml4cvd.metrics import per_class_precision, per_class_precision_3d, per_class_precision_4d, per_class_precision_5d
from ml4cvd.defines import DataSetType, CODING_VALUES_MISSING, TENSOR_MAP_GROUP_MISSING_CONTINUOUS, TENSOR_MAP_GROUP_CONTINUOUS
from ml4cvd.defines import EPS, JOIN_CHAR, IMPUTATION_RANDOM, IMPUTATION_MEAN, CODING_VALUES_LESS_THAN_ONE, MRI_SEGMENTED_CHANNEL_MAP
from ml4cvd.defines import MRI_FRAMES, MRI_SEGMENTED, MRI_TO_SEGMENT, MRI_ZOOM_INPUT, MRI_ZOOM_MASK, MRI_ANNOTATION_NAME, MRI_ANNOTATION_CHANNEL_MAP

np.set_printoptions(threshold=np.inf)

CONTINUOUS_NEVER_ZERO = ['bike_max_hr', 'bike_resting_hr', 'ecg-bike-max-pred-hr-no0',
                         '25006_Volume-of-grey-matter_2', '25021_Volume-of-amygdala-left_2',
                         '25737_Discrepancy-between-dMRI-brain-image-and-T1-brain-image_2', '25738_Discrepancy-between-SWI-brain-image-and-T1-brain-image_2',
                         '25739_Discrepancy-between-rfMRI-brain-image-and-T1-brain-image_2', '25740_Discrepancy-between-tfMRI-brain-image-and-T1-brain-image_2',
                         '25736_Discrepancy-between-T2-FLAIR-brain-image-and-T1-brain-image_2',
                         ]

CONTINUOUS_WITH_CATEGORICAL_ANSWERS = ['92_Operation-yearage-first-occurred_0_0', '1807_Fathers-age-at-death_0_0',
                                       '130_Place-of-birth-in-UK--east-coordinate_0_0',
                                       '87_Noncancer-illness-yearage-first-occurred_0_0',
                                       '1883_Number-of-full-sisters_0_0', '2966_Age-high-blood-pressure-diagnosed_0_0',
                                       '129_Place-of-birth-in-UK--north-coordinate_0_0',
                                       '1070_Time-spent-watching-television-TV_0_0', '1438_Bread-intake_0_0',
                                       '3526_Mothers-age-at-death_0_0',
                                       '2217_Age-started-wearing-glasses-or-contact-lenses_0_0', '1488_Tea-intake_0_0',
                                       '1060_Time-spent-outdoors-in-winter_0_0', '1528_Water-intake_0_0',
                                       '874_Duration-of-walks_0_0', '894_Duration-of-moderate-activity_0_0',
                                       '1458_Cereal-intake_0_0',
                                       '884_Number-of-daysweek-of-moderate-physical-activity-10-minutes_0_0',
                                       '1873_Number-of-full-brothers_0_0', '1845_Mothers-age_0_0',
                                       '1090_Time-spent-driving_0_0', '1289_Cooked-vegetable-intake_0_0',
                                       '3809_Time-since-last-prostate-specific-antigen-PSA-test_0_0',
                                       '1568_Average-weekly-red-wine-intake_0_0', '2897_Age-stopped-smoking_0_0',
                                       '864_Number-of-daysweek-walked-10-minutes_0_0',
                                       '1588_Average-weekly-beer-plus-cider-intake_0_0',
                                       '2355_Most-recent-bowel-cancer-screening_0_0', '2976_Age-diabetes-diagnosed_0_0',
                                       '3761_Age-hay-fever-rhinitis-or-eczema-diagnosed_0_0',
                                       '3786_Age-asthma-diagnosed_0_0',
                                       '1578_Average-weekly-champagne-plus-white-wine-intake_0_0',
                                       '1598_Average-weekly-spirits-intake_0_0',
                                       '1608_Average-weekly-fortified-wine-intake_0_0',
                                       '1299_Salad--raw-vegetable-intake_0_0', '1309_Fresh-fruit-intake_0_0',
                                       '1319_Dried-fruit-intake_0_0', '3680_Age-when-last-ate-meat_0_0',
                                       '914_Duration-of-vigorous-activity_0_0',
                                       '1050_Time-spend-outdoors-in-summer_0_0', '1737_Childhood-sunburn-occasions_0_0',
                                       '1269_Exposure-to-tobacco-smoke-at-home_0_0',
                                       '2867_Age-started-smoking-in-former-smokers_0_0',
                                       '2887_Number-of-cigarettes-previously-smoked-daily_0_0',
                                       '2926_Number-of-unsuccessful-stopsmoking-attempts_0_0',
                                       '2684_Years-since-last-breast-cancer-screening--mammogram_0_0',
                                       '2734_Number-of-live-births_0_0',
                                       '2804_Age-when-last-used-oral-contraceptive-pill_0_0',
                                       '2824_Age-at-hysterectomy_0_0',
                                       '3536_Age-started-hormonereplacement-therapy-HRT_0_0',
                                       '3546_Age-last-used-hormonereplacement-therapy-HRT_0_0',
                                       '3581_Age-at-menopause-last-menstrual-period_0_0',
                                       '3839_Number-of-spontaneous-miscarriages_0_0',
                                       '2405_Number-of-children-fathered_0_0',
                                       '3992_Age-emphysemachronic-bronchitis-diagnosed_0_0',
                                       '4022_Age-pulmonary-embolism-blood-clot-in-lung-diagnosed_0_0',
                                       '4429_Average-monthly-beer-plus-cider-intake_0_0'
                                       ]

MRI_ANNOTATION_GOOD_NEEDED = ['corrected_extracted_lvesv', 'corrected_extracted_lvedv', 'corrected_extracted_lvef', 'sax_inlinevf_zoom_mask',
                              'cine_segmented_sax_inlinevf_segmented', 'mri_systole_diastole_8_segmented', 'mri_systole_diastole_segmented',
                              'mri_slice_segmented'] # , 'sax_all_diastole_segmented'

MERGED_MAPS = ['mothers_age_0', 'fathers_age_0',]
NOT_MISSING = 'not-missing'

MEAN_IDX = 0
STD_IDX = 1


class TensorMap(object):
    """Tensor maps encode the semantics, shapes and types of tensors available

        The mapping can be to numpy nd arrays, categorical labels, or continuous values.
        The tensor shapes can be inferred for categorical TensorMaps which provide a channel mapping dictionary.
        The channel map is a dict mapping a description string to an index into a numpy array i.e the tensor.
        For categorical data the resulting tensor is a one hot vector with a 1 at the channel index and zeros elsewhere.
        In general, new data sources require new TensorMaps and new tensor writers.
        Input and output names are treated differently to allow self mappings, for example autoencoders
    """
    def __init__(self,
                 name,
                 shape=None,
                 group=None,
                 loss=None,
                 model=None,
                 metrics=None,
                 parents=None,
                 sentinel=None,
                 activation=None,
                 loss_weight=1.0,
                 channel_map=None,
                 hd5_override=None,
                 dependent_map=None,
                 required_inputs=None,
                 normalization=None,
                 annotation_units=32,
                 imputation=None,
                 tensor_from_file=None,
                 dtype=None,
                 validator=None,
                 cacheable=True,):
        """TensorMap constructor


        :param name: String name of the tensor mapping
        :param shape: Tuple of integers specifying tensor shape
        :param group: String group of the tensor mapping
        :param loss: Loss function or str specifying pre-defined loss function
        :param model: Model for hidden layer tensor maps
        :param metrics: List of metric functions of strings
        :param parents: List of tensorMaps which must be attached to the graph before this one
        :param sentinel: If set, this value should never naturally occur in this TensorMap, it will be used for masking loss function
        :param activation: String specifying activation function
        :param loss_weight: Relative weight of the loss from this tensormap
        :param channel_map: Dictionary mapping strings indicating channel meaning to channel index integers
        :param hd5_override: Override default behavior of tensor_from_file
        :param dependent_map: TensorMap that depends on or is determined by this one
        :param required_inputs: List of TensorMaps that are required by this one, used by hidden layer TensorMaps
        :param normalization: Dictionary specifying normalization values
        :param annotation_units: Size of embedding dimension for unstructured input tensor maps.
        :param imputation: Method of imputation for missing values. Options are mean or random.
        :param tensor_from_file: Function that returns numpy array from hd5 file for this TensorMap
        :param dtype: DataSetType of tensor map
        """
        self.name = name
        self.loss = loss
        self.model = model
        self.shape = shape
        self.group = group
        self.metrics = metrics
        self.parents = parents
        self.sentinel = sentinel
        self.activation = activation
        self.loss_weight = loss_weight
        self.channel_map = channel_map
        self.hd5_override = hd5_override
        self.normalization = normalization
        self.dependent_map = dependent_map
        self.required_inputs = required_inputs
        self.annotation_units = annotation_units
        self.imputation = imputation
        self.tensor_from_file = tensor_from_file
        self.dtype = dtype
        self.validator = validator
        self.cacheable = cacheable

        if self.shape is None:
            if self.is_multi_field_continuous_with_missing_channel():
                self.shape = (len(channel_map) * 2,)
            else:
                self.shape = (len(channel_map),)

        if self.activation is None and self.is_categorical_any():
            self.activation = 'softmax'

        if self.activation is None and self.is_continuous_any():
            self.activation = 'linear'

        if self.loss is None and self.is_categorical_any():
            self.loss = 'categorical_crossentropy'
        elif self.loss is None and self.is_continuous() and self.sentinel is not None:
            self.loss = sentinel_logcosh_loss(self.sentinel)
        elif self.loss is None and self.is_continuous():
            self.loss = 'mse'
        elif self.loss is None and self.is_proportional_hazard():
            self.loss = survival_likelihood_loss(self.shape[0]//2)
            self.activation = 'sigmoid'
        elif self.loss is None:
            self.loss = 'mse'

        if self.metrics is None and self.is_categorical_any():
            self.metrics = ['categorical_accuracy']
            if len(self.shape) == 1:
                self.metrics += per_class_precision(self.channel_map)
                self.metrics += per_class_recall(self.channel_map)
            elif len(self.shape) == 2:
                self.metrics += per_class_precision_3d(self.channel_map)
                self.metrics += per_class_recall_3d(self.channel_map)
            elif len(self.shape) == 3:
                self.metrics += per_class_precision_4d(self.channel_map)
                self.metrics += per_class_recall_4d(self.channel_map)
            elif len(self.shape) == 4:
                self.metrics += per_class_precision_5d(self.channel_map)
                self.metrics += per_class_recall_5d(self.channel_map)
        elif self.metrics is None and self.is_continuous_any():
            self.metrics = [pearson]
        elif self.metrics is None:
            self.metrics = []

        if self.tensor_from_file is None:
            self.tensor_from_file = _default_tensor_from_file

        if self.validator is None:
            self.validator = lambda tm, x: x

    def __hash__(self):
        return hash((self.name, self.shape, self.group))

    def __eq__(self, other):
        if not isinstance(other, TensorMap):
            return NotImplemented
        else:
            self_attributes = self.__dict__.items()
            other_attributes = other.__dict__.items()

            for (self_field, self_value), (other_field, other_value) in zip(self_attributes, other_attributes):
                if self_field != other_field:
                    return False
                if not _is_equal_field(self_value, other_value):
                    logging.debug(f"Comparing two '{self.name}' tensor maps: "
                                  f"'{self_field}' values '{self_value}' and '{other_value}' are not equal.")
                    return False

            return True

    def output_name(self):
        if self.group is None:
            return JOIN_CHAR.join(['output', self.name])
        else:
            return JOIN_CHAR.join(['output', self.name, self.group])

    def input_name(self):
        if self.group is None:
            return JOIN_CHAR.join(['input', self.name])
        else:
            return JOIN_CHAR.join(['input', self.name, self.group])

    def is_categorical(self):
        return self.group == 'categorical'

    def is_categorical_index(self):
        return self.group == 'categorical_index'

    def is_categorical_date(self):
        return self.group == 'categorical_date'

    def is_categorical_flag(self):
        return self.group == 'categorical_flag'

    def is_categorical_any(self):
        return self.is_categorical_index() or self.is_categorical() or self.is_categorical_date() or self.is_categorical_flag() or self.is_ecg_categorical_interpretation() or self.dtype == DataSetType.CATEGORICAL

    def is_continuous(self):
        return self.group == 'continuous' or self.dtype == DataSetType.CONTINUOUS

    def is_multi_field_continuous(self):
        return self.group == TENSOR_MAP_GROUP_MISSING_CONTINUOUS or self.group == TENSOR_MAP_GROUP_CONTINUOUS

    def is_root_array(self):
        return self.group == 'root_array'

    def is_multi_field_continuous_with_missing_channel(self):
        return self.group == TENSOR_MAP_GROUP_MISSING_CONTINUOUS

    def is_diagnosis_time(self):
        return self.group == 'diagnosis_time'

    def is_continuous_any(self):
        return self.is_continuous() or self.is_diagnosis_time()

    def is_ecg_rest(self):
        return self.group == 'ecg_rest'

    def is_ecg_categorical_interpretation(self):
        return self.group == 'ecg_categorical_interpretation'

    def is_ecg_bike(self):
        return self.group == 'ecg_bike'

    def is_ecg_bike_recovery(self):
        return self.group == 'ecg_bike_recovery'

    def is_ecg_text(self):
        return self.group == 'ecg_text'

    def is_hidden_layer(self):
        return self.group == 'hidden_layer' and self.model is not None

    def is_categorical_any_with_shape_len(self, length):
        return self.is_categorical_any() and len(self.shape) == length

    def is_imputation_random(self):
        return self.is_multi_field_continuous() and self.imputation == IMPUTATION_RANDOM

    def is_imputation_mean(self):
        return self.is_multi_field_continuous() and self.imputation == IMPUTATION_MEAN

    def is_proportional_hazard(self):
        return self.group == 'proportional_hazard'

    def zero_mean_std1(self, np_tensor):
        np_tensor -= np.mean(np_tensor)
        np_tensor /= np.std(np_tensor) + EPS
        np_tensor = np.nan_to_num(np_tensor)
        return np_tensor

    def normalize_and_validate(self, np_tensor):
        self.validator(self, np_tensor)
        if self.normalization is None:
            return np_tensor
        if 'zero_mean_std1' in self.normalization:
            return self.zero_mean_std1(np_tensor)
        if 'mean' in self.normalization and 'std' in self.normalization:
            not_missing_in_channel_map = False
            if self.channel_map is not None:
                not_missing_in_channel_map = NOT_MISSING in self.channel_map
            if self.is_continuous() and not_missing_in_channel_map:
                for i in range(0, len(np_tensor)):
                    if self.channel_map[NOT_MISSING] == i:
                        continue
                    # If the not-missing channel exists in the channel_map and it is marked as "missing" (value of 0)
                    # and the data itself is 0, then overwrite the value with a draw from a N(0,1)
                    if np_tensor[self.channel_map[NOT_MISSING]] == 0 and np_tensor[i] == 0:
                        np_tensor[i] = np.random.normal(1)
                    elif np_tensor[i] == 0:
                        np_tensor[i] -= self.normalization['mean']
                        np_tensor[i] /= (self.normalization['std'] + EPS)
            else:
                np_tensor -= self.normalization['mean']
                np_tensor /= (self.normalization['std'] + EPS)
            return np_tensor

    def normalize_multi_field_continuous(self, np_tensor):
        if self.normalization is None:
            return np_tensor

        for k in self.channel_map:
            idx = self.channel_map[k] * 2
            # If both the value (at idx) and not-missing channel (at idx + 1) are 0 then impute the value.
            if np_tensor[idx + 1] == 0 and np_tensor[idx] == 0:
                np_tensor[idx] = self.impute()
            else:
                np_tensor[idx] -= self.normalization[k][MEAN_IDX]
                np_tensor[idx] /= (self.normalization[k][STD_IDX] + EPS)

        return np_tensor

    def normalize_multi_field_continuous_no_missing_channels(self, np_tensor, missing_array):
        if self.normalization is None:
            return np_tensor

        for k in self.channel_map:
            idx = self.channel_map[k]
            if missing_array[idx]:
                np_tensor[idx] = self.impute()
            else:
                np_tensor[idx] -= self.normalization[k][MEAN_IDX]
                np_tensor[idx] /= (self.normalization[k][STD_IDX] + EPS)

        return np_tensor

    def impute(self):
        if self.normalization is None:
            return ValueError('Imputation requires normalization.')
        if self.is_imputation_random():
            return np.random.normal(1)
        elif self.is_imputation_mean():
            return 0
        else:
            return ValueError('Imputation method unknown.')

    def rescale(self, np_tensor):
        if self.normalization is None:
            return np_tensor
        elif 'mean' in self.normalization and 'std' in self.normalization:
            np_tensor = np.array(np_tensor) * self.normalization['std']
            np_tensor = np.array(np_tensor) + self.normalization['mean']
            return np_tensor
        elif 'zero_mean_std1' in self.normalization:
            return self.zero_mean_std1(np_tensor)
        else:
            return np_tensor

    # Special cases for tensor maps that merge multiple continuous fields (ie combine age of mother with mother's age
    # at death into one channel)
    def _merged_tensor_from_file(self, hd5):
        if self.name == 'mothers_age_0':
            data = np.zeros(self.shape, dtype=np.float32)
            if 'Mother-still-alive_Yes_0_0' in hd5['categorical']:
                if 'mother_alive' in self.channel_map.keys():
                    data[self.channel_map['mother_alive']] = 1
                if '1845_Mothers-age_0_0' in hd5[self.group]:
                    value = hd5[self.group].get('1845_Mothers-age_0_0')[0]
                    if value > 0:
                        data[self.channel_map['mother_age']] = value
                        data[self.channel_map[NOT_MISSING]] = 1
            elif 'Mother-still-alive_No_0_0' in hd5['categorical']:
                if 'mother_dead' in self.channel_map.keys():
                    data[self.channel_map['mother_dead']] = 1
                if '3526_Mothers-age-at-death_0_0' in hd5[self.group]:
                    value = hd5[self.group].get('3526_Mothers-age-at-death_0_0')[0]
                    if value > 0:
                        data[self.channel_map['mother_age']] = value
                        data[self.channel_map[NOT_MISSING]] = 1
            return self.normalize_and_validate(data)
        elif self.name == 'fathers_age_0':
            data = np.zeros(self.shape, dtype=np.float32)
            if 'Father-still-alive_Yes_0_0' in hd5['categorical']:
                if 'father_alive' in self.channel_map.keys():
                    data[self.channel_map['father_alive']] = 1
                if '2946_Fathers-age_0_0' in hd5[self.group]:
                    value = hd5[self.group].get('2946_Fathers-age_0_0')[0]
                    if value > 0:
                        data[self.channel_map['father_age']] = value
                        data[self.channel_map[NOT_MISSING]] = 1
            elif 'Father-still-alive_No_0_0' in hd5['categorical']:
                if 'father_dead' in self.channel_map.keys():
                    data[self.channel_map['father_dead']] = 1
                if '1807_Fathers-age-at-death_0_0' in hd5[self.group]:
                    value = hd5[self.group].get('1807_Fathers-age-at-death_0_0')[0]
                    if value > 0:
                        data[self.channel_map['father_age']] = value
                        data[self.channel_map[NOT_MISSING]] = 1
            return self.normalize_and_validate(data)
        raise ValueError('No Merged Tensor Map handling found for ' + self.name + ".")


def make_range_validator(minimum: float, maximum: float):
    def _range_validator(tm: TensorMap, tensor: np.ndarray):
        if not ((tensor > minimum).all() and (tensor < maximum).all()):
            raise ValueError(f'TensorMap {tm.name} failed range check.')
    return _range_validator


def no_nans(tm: TensorMap, tensor: np.ndarray):
    if np.isnan(tensor).any():
        raise ValueError(f'Skipping TensorMap {tm.name} with NaNs.')


def _translate(val, cur_min, cur_max, new_min, new_max):
    val -= cur_min
    val /= (cur_max - cur_min)
    val *= (new_max - new_min)
    val += new_min
    return val


def str2date(d):
    parts = d.split('-')
    if len(parts) < 2:
        return datetime.datetime.now().date()
    return datetime.date(int(parts[0]), int(parts[1]), int(parts[2]))


def _is_equal_field(field1: Any, field2: Any) -> bool:
    """We consider two fields equal if
            a. they are not functions and they are equal, or
            b. one or both are functions and their names match
        If the fields are lists, we check for the above equality for corresponding
        elements from the list.
    """
    if isinstance(field1, list) and isinstance(field2, list):
        if len(field1) != len(field2):
            return False
        elif len(field1) == 0:
            return True

        fields1 = map(_get_name_if_function, field1)
        fields2 = map(_get_name_if_function, field2)

        return all([f1 == f2] for f1, f2 in zip(sorted(fields1), sorted(fields2)))
    else:
        return _get_name_if_function(field1) == _get_name_if_function(field2)


def _get_name_if_function(field: Any) -> Any:
    """We assume 'field' is a function if it's 'callable()'"""
    if callable(field):
        return field.__name__
    else:
        return field


def _default_tensor_from_file(tm, hd5, dependents={}):
    """Reconstruct a tensor from an hd5 file

    Arguments
        tm: The TensorMap that describes the type of tensor to make
        hd5: The file where the tensor was saved
        dependents: A dict that maps dependent TensorMaps to numpy arrays
            if self has a dependent TensorMap it will be constructed and added here

    Returns
        A numpy array whose dimension and type is dictated by tm
    """
    if tm.is_categorical_index():
        categorical_data = np.zeros(tm.shape, dtype=np.float32)
        if tm.name in hd5:
            index = int(hd5[tm.name][0])
            categorical_data[index] = 1.0
        elif tm.name in hd5['categorical']:
            index = int(hd5['categorical'][tm.name][0])
            categorical_data[index] = 1.0
        else:
            raise ValueError(f"No categorical index found for tensor map: {tm.name}.")
        return categorical_data
    elif tm.is_categorical_flag():
        categorical_data = np.zeros(tm.shape, dtype=np.float32)
        index = 0
        if tm.name in hd5 and int(hd5[tm.name][0]) != 0:
            index = 1
        elif tm.name in hd5['categorical'] and int(hd5['categorical'][tm.name][0]) != 0:
            index = 1
        categorical_data[index] = 1.0
        return categorical_data
    elif tm.is_categorical_date():
        categorical_data = np.zeros(tm.shape, dtype=np.float32)
        if tm.name in hd5:
            index = int(hd5[tm.name][0])
        elif tm.name in hd5['categorical']:
            index = int(hd5['categorical'][tm.name][0])
        else:
            index = 0  # Assume no disease if the tensor does not have the dataset
        if index != 0:
            if tm.name + '_date' in hd5:
                disease_date = str2date(str(hd5[tm.name + '_date'][0]))
                assess_date = str2date(str(hd5['assessment-date_0_0'][0]))
            elif tm.name + '_date' in hd5['dates']:
                disease_date = str2date(str(hd5['dates'][tm.name + '_date'][0]))
                assess_date = str2date(str(hd5['dates']['enroll_date'][0]))
            else:
                raise ValueError(f"No date found for tensor map: {tm.name}.")
            index = 1 if disease_date < assess_date else 2
        categorical_data[index] = 1.0
        return categorical_data
    elif tm.is_diagnosis_time():
        time_data = np.zeros((1,), dtype=np.float32)
        disease_status = int(hd5[tm.name][0])
        assess_date = str2date(str(hd5['assessment-date_0_0'][0]))
        disease_date = str2date(str(hd5[tm.name + '_date'][0]))
        delta = relativedelta.relativedelta(disease_date, assess_date)
        difference = (delta.years * 12) + delta.months
        if disease_status == 0 or difference == 0:
            raise ValueError('Ignoring healthy people in diagnosis time.')
        time_data[0] = difference
        return time_data
    elif tm.name in [MRI_TO_SEGMENT, MRI_ZOOM_INPUT]:
        mask_group = MRI_SEGMENTED if tm.name == MRI_TO_SEGMENT else MRI_ZOOM_MASK
        slice_idx = np.random.choice(list(hd5[tm.name].keys()))
        angle_idx = int(slice_idx) // MRI_FRAMES
        tensor = np.zeros(tm.shape, dtype=np.float32)
        dependents[tm.dependent_map] = np.zeros(tm.dependent_map.shape, dtype=np.float32)
        for i in range(tm.shape[-2]):
            cur_slice = str((angle_idx * MRI_FRAMES) + i + 1)  # Instance Number off by 1
            tensor[:, :, i, 0] = np.array(hd5[tm.name].get(cur_slice), dtype=np.float32)
            label_tensor = np.array(hd5[mask_group].get(cur_slice), dtype=np.float32)
            dependents[tm.dependent_map][:, :, i, :] = to_categorical(label_tensor, tm.dependent_map.shape[-1])
        return tm.normalize_and_validate(tensor)
    elif tm.name == 'aligned_distance':
        return np.zeros((1,), dtype=np.float32)
    elif tm.name == 'lms_ideal_optimised_low_flip_6dyn_4slice':
        whole_liver = np.array(hd5['lms_ideal_optimised_low_flip_6dyn'])
        cur_index = np.random.randint(whole_liver.shape[2] - tm.shape[2])
        tensor = whole_liver[:, :, cur_index:cur_index + 4]
        return tm.normalize_and_validate(tensor)
    elif tm.name == 'mri_slice':
        cur_slice = np.random.choice(list(hd5[MRI_TO_SEGMENT].keys()))
        tensor = np.zeros(tm.shape, dtype=np.float32)
        dependents[tm.dependent_map] = np.zeros(tm.dependent_map.shape, dtype=np.float32)
        tensor[:, :, 0] = np.array(hd5[MRI_TO_SEGMENT].get(cur_slice), dtype=np.float32)
        label_tensor = np.array(hd5[MRI_SEGMENTED].get(cur_slice), dtype=np.float32)
        dependents[tm.dependent_map][:, :, :] = to_categorical(label_tensor, tm.dependent_map.shape[-1])
        return tm.normalize_and_validate(tensor)
    elif tm.name == 'sax_inlinevf_zoom_blackout':
        mask_group = MRI_ZOOM_MASK
        slice_idx = np.random.choice(list(hd5[MRI_ZOOM_INPUT].keys()))
        angle_idx = int(slice_idx) // MRI_FRAMES
        tensor = np.zeros(tm.shape, dtype=np.float32)
        dependents[tm.dependent_map] = np.zeros(tm.dependent_map.shape, dtype=np.float32)
        for i in range(tm.shape[-2]):
            cur_slice = str((angle_idx * MRI_FRAMES) + i + 1)  # Instance Number off by 1
            tensor[:, :, i, 0] = np.array(hd5[MRI_ZOOM_INPUT].get(cur_slice), dtype=np.float32)
            label_tensor = np.array(hd5[mask_group].get(cur_slice), dtype=np.float32)
            dependents[tm.dependent_map][:, :, i, :] = to_categorical(label_tensor, tm.dependent_map.shape[-1])
            tensor[:, :, i, 0] *= np.not_equal(label_tensor, 0, dtype=np.float32)
        return tm.normalize_and_validate(tensor)
    elif tm.name == 'cine_segmented_sax_inlinevf_blackout':
        mask_group = MRI_SEGMENTED
        slice_idx = np.random.choice(list(hd5[MRI_TO_SEGMENT].keys()))
        angle_idx = int(slice_idx) // MRI_FRAMES
        tensor = np.zeros(tm.shape, dtype=np.float32)
        dependents[tm.dependent_map] = np.zeros(tm.dependent_map.shape, dtype=np.float32)
        for i in range(tm.shape[-2]):
            cur_slice = str((angle_idx * MRI_FRAMES) + i + 1)  # Instance Number off by 1
            tensor[:, :, i, 0] = np.array(hd5[MRI_TO_SEGMENT].get(cur_slice), dtype=np.float32)
            label_tensor = np.array(hd5[mask_group].get(cur_slice), dtype=np.float32)
            dependents[tm.dependent_map][:, :, i, :] = to_categorical(label_tensor, tm.dependent_map.shape[-1])
            tensor[:, :, i, 0] *= np.not_equal(label_tensor, 0, dtype=np.float32)
        return tm.normalize_and_validate(tensor)
    elif tm.name == 'mri_systole_diastole':
        if tm.hd5_override is not None:
            b_number = 'b' + str(np.random.choice(tm.hd5_override))
            diastole_slice = 'diastole_frame_' + b_number
        else:
            frames = [str(frame) for frame in hd5.keys() if 'diastole_frame_' in str(frame)]
            if len(frames) == 0:
                raise ValueError('No diastole frames found.')
            diastole_slice = np.random.choice(frames)
            b_number = diastole_slice.split('_')[-1]  # (e.g b1, b2, b3, ... b12 ish)
        tensor = np.zeros(tm.shape, dtype=np.float32)
        dependents[tm.dependent_map] = np.zeros(tm.dependent_map.shape, dtype=np.float32)
        tensor[:, :, 0, 0] = np.array(hd5[diastole_slice], dtype=np.float32)
        dependents[tm.dependent_map][:, :, 0, :] = to_categorical(np.array(hd5['diastole_mask_' + b_number]), tm.dependent_map.shape[-1])
        tensor[:, :, 1, 0] = np.array(hd5['systole_frame_' + b_number], dtype=np.float32)
        dependents[tm.dependent_map][:, :, 1, :] = to_categorical(np.array(hd5['systole_mask_' + b_number]), tm.dependent_map.shape[-1])
        return tm.normalize_and_validate(tensor)
    elif tm.name == 'mri_systole_diastole_8':
        tensor = np.zeros(tm.shape, dtype=np.float32)
        dependents[tm.dependent_map] = np.zeros(tm.dependent_map.shape, dtype=np.float32)
        tensor[:, :, 0, 0] = np.array(hd5['diastole_frame_b2'], dtype=np.float32)
        dependents[tm.dependent_map][:, :, 0, :] = to_categorical(np.array(hd5['diastole_mask_b2']), tm.dependent_map.shape[-1])
        tensor[:, :, 1, 0] = np.array(hd5['systole_frame_b2'], dtype=np.float32)
        dependents[tm.dependent_map][:, :, 1, :] = to_categorical(np.array(hd5['systole_mask_b2']), tm.dependent_map.shape[-1])
        tensor[:, :, 2, 0] = np.array(hd5['diastole_frame_b4'], dtype=np.float32)
        dependents[tm.dependent_map][:, :, 2, :] = to_categorical(np.array(hd5['diastole_mask_b4']), tm.dependent_map.shape[-1])
        tensor[:, :, 3, 0] = np.array(hd5['systole_frame_b4'], dtype=np.float32)
        dependents[tm.dependent_map][:, :, 3, :] = to_categorical(np.array(hd5['systole_mask_b4']), tm.dependent_map.shape[-1])
        tensor[:, :, 4, 0] = np.array(hd5['diastole_frame_b6'], dtype=np.float32)
        dependents[tm.dependent_map][:, :, 4, :] = to_categorical(np.array(hd5['diastole_mask_b6']), tm.dependent_map.shape[-1])
        tensor[:, :, 5, 0] = np.array(hd5['systole_frame_b6'], dtype=np.float32)
        dependents[tm.dependent_map][:, :, 5, :] = to_categorical(np.array(hd5['systole_mask_b6']), tm.dependent_map.shape[-1])
        tensor[:, :, 6, 0] = np.array(hd5['diastole_frame_b8'], dtype=np.float32)
        dependents[tm.dependent_map][:, :, 6, :] = to_categorical(np.array(hd5['diastole_mask_b8']), tm.dependent_map.shape[-1])
        tensor[:, :, 7, 0] = np.array(hd5['systole_frame_b8'], dtype=np.float32)
        dependents[tm.dependent_map][:, :, 7, :] = to_categorical(np.array(hd5['systole_mask_b8']), tm.dependent_map.shape[-1])
        return tm.normalize_and_validate(tensor)
    elif tm.name == 'sax_all_diastole':
        missing = 0
        tensor = np.zeros(tm.shape, dtype=np.float32)
        dependents[tm.dependent_map] = np.zeros(tm.dependent_map.shape, dtype=np.float32)
        for b in range(tm.shape[-2]):
            try:
                tensor[:, :, b, 0] = np.array(hd5[f'diastole_frame_b{b}'], dtype=np.float32)
                dependents[tm.dependent_map][:, :, b, :] = to_categorical(np.array(hd5[f'diastole_mask_b{b}']), tm.dependent_map.shape[-1])
            except KeyError:
                missing += 1
                tensor[:, :, b, 0] = 0
                dependents[tm.dependent_map][:, :, b, MRI_SEGMENTED_CHANNEL_MAP['background']] = 1
        if missing == tm.shape[-2]:
            raise ValueError(f'Could not find any slices in {tm.name} was hoping for {tm.shape[-2]}')
        return tm.normalize_and_validate(tensor)
    elif tm.is_root_array():
        tensor = np.zeros(tm.shape, dtype=np.float32)
        tensor[:] = np.array(hd5[tm.name], dtype=np.float32)
        return tm.normalize_and_validate(tensor)
    elif tm.name in MRI_ANNOTATION_GOOD_NEEDED:
        continuous_data = np.zeros(tm.shape, dtype=np.float32)  # Automatic left ventricular analysis with InlineVF
        if MRI_ANNOTATION_NAME in hd5['categorical'] and hd5['categorical'][MRI_ANNOTATION_NAME][0] != MRI_ANNOTATION_CHANNEL_MAP['good']:
            raise ValueError('MRI Critic annotation not good or unreviewed.')
        continuous_data[0] = float(hd5['continuous'][tm.name][0])
        if continuous_data[0] == 0 and (tm.sentinel == None and tm.name in CONTINUOUS_NEVER_ZERO):
            raise ValueError(tm.name + ' is a continuous value that cannot be set to 0, but no value was found.')
        if continuous_data[0] == 0 and tm.sentinel is not None:
            continuous_data[:] = tm.sentinel
            return continuous_data
        return tm.normalize_and_validate(continuous_data)
    elif tm.name == 'ecg_coarse':
        categorical_data = np.zeros(tm.shape, dtype=np.float32)
        if 'poor_data_quality' in hd5['categorical']:
            raise ValueError('Poor data skipped by ecg_coarse.')
        ecg_interpretation = str(hd5['ecg_rest_text'][0])
        for afib in ['Atrial fibrillation']:
            if afib in ecg_interpretation:
                categorical_data[tm.channel_map['Atrial_fibrillation']] = 1.0
                return categorical_data
        for rhythm in ['sinus', 'Sinus']:
            if rhythm in ecg_interpretation:
                categorical_data[tm.channel_map['Sinus_rhythm']] = 1.0
                return categorical_data
        categorical_data[tm.channel_map['Other_rhythm']] = 1.0
        return categorical_data
    elif tm.name == 'ecg_semi_coarse':
        categorical_data = np.zeros(tm.shape, dtype=np.float32)
        if 'poor_data_quality' in hd5['categorical']:
            raise ValueError('Poor data skipped by ecg_coarse.')
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
    elif tm.name == 'ecg_semi_coarse_with_poor':
        categorical_data = np.zeros(tm.shape, dtype=np.float32)
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
    elif tm.is_ecg_categorical_interpretation():
        categorical_data = np.zeros(tm.shape, dtype=np.float32)
        for channel in tm.channel_map:
            if channel in str(hd5['ecg_rest_text'][0]):
                categorical_data[tm.channel_map[channel]] = 1.0
                return categorical_data
        if 'no_' + tm.name in tm.channel_map:
            categorical_data[tm.channel_map['no_' + tm.name]] = 1.0
            return categorical_data
        else:
            raise ValueError(f"ECG categorical interpretation could not find any of these keys: {tm.channel_map.keys()}")
    elif tm.is_categorical() and tm.channel_map is not None:
        categorical_data = np.zeros(tm.shape, dtype=np.float32)
        for channel in tm.channel_map:
            if channel in hd5['categorical']:
                categorical_data[tm.channel_map[channel]] = 1.0
        return categorical_data
    elif tm.name in MERGED_MAPS:
        return tm._merged_tensor_from_file(hd5)
    elif tm.is_continuous():
        continuous_data = np.zeros(tm.shape, dtype=np.float32)
        if tm.name in hd5:
            if hasattr(hd5[tm.name], "__shape__"):
                continuous_data[0] = hd5[tm.name][0]
            else:
                continuous_data[0] = hd5[tm.name][()]
        missing = True
        for k in tm.channel_map:
            if k in hd5[tm.group]:
                value = hd5[tm.group][k][0]
                missing = False
                if k in CONTINUOUS_WITH_CATEGORICAL_ANSWERS:
                    if value in CODING_VALUES_LESS_THAN_ONE:
                        value = .5
                    if value in CODING_VALUES_MISSING:
                        # need to set missing values to 0 so normalization works
                        value = 0
                        missing = True
                continuous_data[tm.channel_map[k]] = value
        if NOT_MISSING in tm.channel_map and not missing:
            continuous_data[tm.channel_map[NOT_MISSING]] = 1
        if continuous_data[0] == 0 and (tm.sentinel is None and tm.name in CONTINUOUS_NEVER_ZERO):
            raise ValueError(tm.name + ' is a continuous value that cannot be set to 0, but no value was found.')
        return tm.normalize_and_validate(continuous_data)
    elif tm.is_multi_field_continuous():
        if tm.is_multi_field_continuous_with_missing_channel():
            multiplier = 2
        else:
            multiplier = 1
        missing_array = [False] * len(tm.channel_map)
        continuous_data = np.zeros(tm.shape, dtype=np.float32)
        for k in tm.channel_map:
            missing = True
            if k in hd5['continuous']:
                value = hd5['continuous'][k][0]
                missing = False
                if tm.name in CONTINUOUS_WITH_CATEGORICAL_ANSWERS:
                    if value in CODING_VALUES_LESS_THAN_ONE:
                        value = .5
                    if value in CODING_VALUES_MISSING:
                        # need to set missing values to 0 so normalization works
                        value = 0
                        missing = True
                # Put value at index k (times 2 to make space for the not-missing channels), and put whether or not
                # this value is not missing in the following element.
                continuous_data[tm.channel_map[k] * multiplier] = value
            if tm.is_multi_field_continuous_with_missing_channel():
                continuous_data[tm.channel_map[k] * multiplier + 1] = not missing
            else:
                missing_array[tm.channel_map[k]] = missing
        if tm.is_multi_field_continuous_with_missing_channel():
            return tm.normalize_multi_field_continuous(continuous_data)
        else:
            return tm.normalize_multi_field_continuous_no_missing_channels(continuous_data, missing_array)
    elif tm.is_ecg_bike():
        tensor = np.array(hd5[tm.group][tm.name], dtype=np.float32)
        return tm.normalize_and_validate(tensor)
    elif tm.is_ecg_bike_recovery():
        tensor = np.zeros(tm.shape)
        for channel, idx in tm.channel_map.items():
            tensor[:, idx] = hd5[tm.group][channel]
        return tm.normalize_and_validate(tensor)
    elif tm.is_ecg_text():
        tensor = np.zeros(tm.shape, dtype=np.float32)
        dependents[tm.dependent_map] = np.zeros(tm.dependent_map.shape, dtype=np.float32)
        caption = str(hd5[tm.name][0]).strip()
        char_idx = np.random.randint(len(caption) + 1)
        if char_idx == len(caption):
            next_char = '!'
        else:
            next_char = caption[char_idx]
        dependents[tm.dependent_map][tm.dependent_map.channel_map[next_char]] = 1.0
        window_offset = max(0, tm.shape[0] - char_idx)

        for k in range(max(0, char_idx - tm.shape[0]), char_idx):
            tensor[window_offset, tm.dependent_map.channel_map[caption[k]]] = 1.0
            window_offset += 1
        return tensor
    elif tm.is_hidden_layer():
        input_dict = {}
        for input_tm in tm.required_inputs:
            input_dict[input_tm.input_name()] = np.expand_dims(input_tm.tensor_from_file(input_tm, hd5), axis=0)
        return tm.model.predict(input_dict)
    elif tm.dependent_map is not None:  # Assumes dependent maps are 1-hot categoricals
        dataset_key = np.random.choice(list(tm.dependent_map.channel_map.keys()))
        one_hot = np.zeros(tm.dependent_map.shape)
        one_hot[tm.dependent_map.channel_map[dataset_key]] = 1.0
        dependents[tm.dependent_map] = one_hot
        tensor = np.array(hd5.get(dataset_key), dtype=np.float32)
        return tm.normalize_and_validate(tensor)
    elif tm.channel_map is None and tm.dependent_map is None:
        return np.array(hd5.get(tm.name), dtype=np.float32)
