import h5py
import numpy as np
import logging
from typing import List, Tuple
from ml4h.tensormap.general import tensor_path
from ml4h.TensorMap import TensorMap, Interpretation, str2date, make_range_validator
from ml4h.defines import StorageType


def is_genetic_man(hd5):
    return 'Genetic-sex_Male_0_0' in hd5['categorical']


def is_genetic_woman(hd5):
    return 'Genetic-sex_Female_0_0' in hd5['categorical']


def age_in_years_tensor(
    date_key,
    birth_key='continuous/34_Year-of-birth_0_0',
    population_normalize=False,
):
    def age_at_tensor_from_file(
        tm: TensorMap,
        hd5: h5py.File,
        dependents=None,
    ):
        try:
            age = np.array([hd5['ecg/latest/patient_info/Age'][()]])
        except:
            logging.info('could not get age')
            raise KeyError('cold not')
        # age = age.astype("float")

        return age
        # return tm.normalize_and_validate(np.array([assess_date.year-birth_year]))

    return age_at_tensor_from_file


def prevalent_incident_tensor(start_date_key, event_date_key):
    def _prevalent_incident_tensor_from_file(
        tm: TensorMap,
        hd5: h5py.File,
        dependents=None,
    ):
        index = 0
        categorical_data = np.zeros(tm.shape, dtype=np.float32)
        if tm.hd5_key_guess() in hd5:
            data = tm.hd5_first_dataset_in_group(hd5, tm.hd5_key_guess())
            if tm.storage_type == StorageType.CATEGORICAL_INDEX or tm.storage_type == StorageType.CATEGORICAL_FLAG:
                index = int(data[0])
                categorical_data[index] = 1.0
            else:
                categorical_data = np.array(data)
        elif tm.storage_type == StorageType.CATEGORICAL_FLAG:
            categorical_data[index] = 1.0
        else:
            raise ValueError(
                f"No HD5 Key at prefix {tm.path_prefix} found for tensor map: {tm.name}.",
            )

        if index != 0:
            if event_date_key in hd5 and start_date_key in hd5:
                disease_date = str2date(str(hd5[event_date_key][0]))
                assess_date = str2date(str(hd5[start_date_key][0]))
            else:
                raise ValueError(f"No date found for tensor map: {tm.name}.")
            index = 1 if disease_date < assess_date else 2
        categorical_data[index] = 1.0
        return categorical_data

    return _prevalent_incident_tensor_from_file


def preprocess_with_function(fxn, hd5_key=None):
    def preprocess_tensor_from_file(tm, hd5, dependents={}):
        missing = True
        continuous_data = np.zeros(tm.shape, dtype=np.float32)
        my_key = tm.hd5_key_guess() if hd5_key is None else hd5_key
        if my_key in hd5:
            missing = False
            continuous_data[0] = tm.hd5_first_dataset_in_group(hd5, my_key)[0]
        if missing and tm.sentinel is None:
            raise ValueError(
                f'No value found for {tm.name}, a continuous TensorMap with no sentinel value, and channel keys:{list(tm.channel_map.keys())}.',
            )
        elif missing:
            continuous_data[:] = tm.sentinel
        return fxn(continuous_data)

    return preprocess_tensor_from_file


def _weekly_alcohol(instance):
    alcohol_keys = [
        f'1568_Average-weekly-red-wine-intake_{instance}_0',
        f'1578_Average-weekly-champagne-plus-white-wine-intake_{instance}_0',
        f'1588_Average-weekly-beer-plus-cider-intake_{instance}_0',
        f'1598_Average-weekly-spirits-intake_{instance}_0',
        f'1608_Average-weekly-fortified-wine-intake_{instance}_0',
    ]

    def alcohol_from_file(tm, hd5, dependents={}):
        drinks = 0
        for k in alcohol_keys:
            data = tm.hd5_first_dataset_in_group(
                hd5, key_prefix=f'{tm.path_prefix}/{k}',
            )
            drinks += float(data[0])
        return np.array([drinks], dtype=np.float32)

    return alcohol_from_file


log_25781_2 = TensorMap(
    '25781_Total-volume-of-white-matter-hyperintensities-from-T1-and-T2FLAIR-images_2_0',
    loss='logcosh',
    path_prefix='continuous',
    normalization={
        'mean': 7,
        'std': 8,
    },
    tensor_from_file=preprocess_with_function(np.log),
    channel_map={'white-matter-hyper-intensities': 0},
)

weight_lbs_2 = TensorMap(
    'weight_lbs',
    Interpretation.CONTINUOUS,
    normalization={
        'mean': 168.74,
        'std': 34.1,
    },
    loss='logcosh',
    channel_map={'weight_lbs': 0},
    tensor_from_file=preprocess_with_function(
        lambda x: x * 2.20462,
        'continuous/21002_Weight_2_0',
    ),
)

weekly_alcohol_0 = TensorMap(
    'weekly_alcohol_0',
    loss='logcosh',
    path_prefix='continuous',
    channel_map={'weekly_alcohol_0': 0},
    tensor_from_file=_weekly_alcohol(0),
)
weekly_alcohol_1 = TensorMap(
    'weekly_alcohol_1',
    loss='logcosh',
    path_prefix='continuous',
    channel_map={'weekly_alcohol_1': 0},
    tensor_from_file=_weekly_alcohol(1),
)
weekly_alcohol_2 = TensorMap(
    'weekly_alcohol_2',
    loss='logcosh',
    path_prefix='continuous',
    channel_map={'weekly_alcohol_2': 0},
    tensor_from_file=_weekly_alcohol(2),
)

###
weight_kg = TensorMap('weight_kg',  Interpretation.CONTINUOUS, normalization={'mean': 76.54286701805927, 'std': 15.467605416933122}, loss='logcosh', channel_map={'weight_kg': 0})
height_cm = TensorMap('height_cm',  Interpretation.CONTINUOUS, normalization={'mean': 169.18064748408653, 'std': 9.265265197273026}, loss='logcosh', channel_map={'height_cm': 0})
bmi_bsa = TensorMap('bmi',  Interpretation.CONTINUOUS, normalization={'mean': 26.65499238706321, 'std': 4.512077188749083}, loss='logcosh', channel_map={'bmi': 0})

mothers_age = TensorMap(
    'mothers_age_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'mother_age': 0, 'mother_alive': 2, 'mother_dead': 3, 'not-missing': 1},
    normalization={'mean': 75.555, 'std': 11.977}, annotation_units = 4,
)

fathers_age = TensorMap(
    'fathers_age_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'father_age': 0, 'father_alive': 2, 'father_dead': 3, 'not-missing': 1},
    normalization={'mean':70.928, 'std': 12.746}, annotation_units = 4,
)

genetic_sex = TensorMap(
    'Genetic-sex_Male_0_0', Interpretation.CATEGORICAL, storage_type=StorageType.CATEGORICAL_FLAG, path_prefix='categorical', annotation_units=2,
    channel_map={'Genetic-sex_Female_0_0': 0, 'Genetic-sex_Male_0_0': 1}, loss='categorical_crossentropy',
)
sex = TensorMap(
    'Sex_Male_0_0', Interpretation.CATEGORICAL, storage_type=StorageType.CATEGORICAL_FLAG, path_prefix='categorical', annotation_units=2,
    channel_map={'Sex_Female_0_0': 0, 'Sex_Male_0_0': 1}, loss='categorical_crossentropy',
)
bmi = TensorMap(
    '23104_Body-mass-index-BMI_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'23104_Body-mass-index-BMI_0_0': 0}, annotation_units=1,
    validator=make_range_validator(0, 300), normalization={'mean': 27.432061533712652, 'std': 4.785244772462738}, loss='logcosh',
)
bmi_ukb = TensorMap(
    'bmi', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'23104_Body-mass-index-BMI_0_0': 0}, annotation_units=1,
    validator=make_range_validator(0, 300), normalization={'mean': 27.432061533712652, 'std': 4.785244772462738}, loss='logcosh',
)
bmi_21 = TensorMap(
    '21001_Body-mass-index-BMI_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'21001_Body-mass-index-BMI_0_0': 0}, annotation_units=1,
    validator=make_range_validator(0, 300), normalization={'mean': 27.3397, 'std': 4.7721}, loss='logcosh',
)
birth_year = TensorMap(
    '22200_Year-of-birth_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'22200_Year-of-birth_0_0': 0}, annotation_units=1, loss='logcosh',
    validator=make_range_validator(1901, 2025), normalization={'mean': 1952.0639129359386, 'std': 7.656326148519739},
)
birth_year_34 = TensorMap(
    '34_Year-of-birth_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'34_Year-of-birth_0_0': 0}, annotation_units=1, loss='logcosh',
    validator=make_range_validator(1901, 2025), normalization = {'mean': 1952.0639129359386, 'std': 7.656326148519739},
)
age_0 = TensorMap(
    '21003_Age-when-attended-assessment-centre_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh', validator=make_range_validator(1, 120),
    normalization={'mean': 56.52847159208494, 'std': 8.095287610193827}, channel_map={'21003_Age-when-attended-assessment-centre_0_0': 0},
)
age_1 = TensorMap(
    '21003_Age-when-attended-assessment-centre_1_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh', validator=make_range_validator(1, 120),
    normalization={'mean': 61.4476555588322, 'std': 7.3992113757847005}, channel_map={'21003_Age-when-attended-assessment-centre_1_0': 0},
)
age_2 = TensorMap(
    '21003_Age-when-attended-assessment-centre_2_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh', validator=make_range_validator(1, 120),
    normalization={'mean': 63.35798891483556, 'std': 7.554638350423902}, channel_map={'21003_Age-when-attended-assessment-centre_2_0': 0},
)

brain_volume = TensorMap(
    '25010_Volume-of-brain-greywhite-matter_2_0', Interpretation.CONTINUOUS, path_prefix='continuous', normalization={'mean': 1165940.0, 'std': 111511.0},
    channel_map={'25010_Volume-of-brain-greywhite-matter_2_0': 0}, loss='logcosh', loss_weight=0.1,
)

sodium = TensorMap(
    '30530_Sodium-in-urine_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'30530_Sodium-in-urine_0_0': 0},
    normalization={'mean': 77.45323967267045, 'std': 44.441236848463774}, annotation_units=1, loss='logcosh',
)
potassium = TensorMap(
    '30520_Potassium-in-urine_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'30520_Potassium-in-urine_0_0': 0},
    normalization={'mean': 63.06182700345117, 'std': 33.84208704773539}, annotation_units=1, loss='logcosh',
)
cholesterol_hdl = TensorMap(
    '30760_HDL-cholesterol_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'30760_HDL-cholesterol_0_0': 0},
    normalization={'mean': 1.4480129055069355, 'std': 0.3823115953478376}, annotation_units=1, loss='logcosh',
)
cholesterol = TensorMap(
    '30690_Cholesterol_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'30690_Cholesterol_0_0': 0},
    normalization={'mean': 5.692381214399044, 'std': 1.1449409331668705}, annotation_units=1, loss='logcosh',
)

cigarettes = TensorMap('2887_Number-of-cigarettes-previously-smoked-daily_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'2887_Number-of-cigarettes-previously-smoked-daily_0_0': 0}, normalization = {'mean': 18.92662147068755, 'std':10.590930376362259 }, annotation_units=1)
alcohol = TensorMap('5364_Average-weekly-intake-of-other-alcoholic-drinks_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'5364_Average-weekly-intake-of-other-alcoholic-drinks_0_0': 0}, normalization = {'mean': 0.03852570253005904, 'std':0.512608370266108 }, annotation_units=1)


def alcohol_channel_map(instance=0, array_idx=0):
    return {
        f'Alcohol-intake-frequency_Never_{instance}_{array_idx}': 0,
        f'Alcohol-intake-frequency_Special-occasions-only_{instance}_{array_idx}': 1,
        f'Alcohol-intake-frequency_One-to-three-times-a-month_{instance}_{array_idx}': 2,
        f'Alcohol-intake-frequency_Once-or-twice-a-week_{instance}_{array_idx}': 3,
        f'Alcohol-intake-frequency_Three-or-four-times-a-week_{instance}_{array_idx}': 4,
        f'Alcohol-intake-frequency_Daily-or-almost-daily_{instance}_{array_idx}': 5,
    }


alcohol_0 = TensorMap('alcohol_0', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_channel_map(instance=0))
alcohol_1 = TensorMap('alcohol_1', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_channel_map(instance=1))
alcohol_2 = TensorMap('alcohol_2', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_channel_map(instance=2))


def alcohol_status_map(instance=0, array_idx=0):
    return {
        f'Alcohol-drinker-status_Never_{instance}_{array_idx}': 0,
        f'Alcohol-drinker-status_Previous_{instance}_{array_idx}': 1,
        f'Alcohol-drinker-status_Current_{instance}_{array_idx}': 2,
    }


alcohol_status_0 = TensorMap('alcohol_status_0', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_status_map(instance=0))
alcohol_status_1 = TensorMap('alcohol_status_1', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_status_map(instance=1))
alcohol_status_2 = TensorMap('alcohol_status_2', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_status_map(instance=2))


def alcohol_meals_map(instance=0, array_idx=0):
    return {
        f'Alcohol-usually-taken-with-meals_No_{instance}_{array_idx}': 0,
        f'Alcohol-usually-taken-with-meals_It-varies_{instance}_{array_idx}': 1,
        f'Alcohol-usually-taken-with-meals_Yes_{instance}_{array_idx}': 2,
    }


alcohol_meals_0 = TensorMap('alcohol_meals_0', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_meals_map(instance=0))
alcohol_meals_1 = TensorMap('alcohol_meals_1', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_meals_map(instance=1))
alcohol_meals_2 = TensorMap('alcohol_meals_2', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_meals_map(instance=2))

coffee = TensorMap(
    '1498_Coffee-intake_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'1498_Coffee-intake_0_0': 0},
    normalization={'mean': 2.015086529948216, 'std': 2.0914960998390497}, annotation_units=1,
)
water = TensorMap(
    '1528_Water-intake_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'1528_Water-intake_0_0': 0},
    normalization={'mean': 2.7322977785723324, 'std': 2.261996814128837}, annotation_units=1,
)
meat = TensorMap(
    '3680_Age-when-last-ate-meat_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'3680_Age-when-last-ate-meat_0_0': 0},
    normalization={'mean': 29.74062983480561, 'std': 14.417292213873964}, annotation_units=1,
)
walks = TensorMap(
    '864_Number-of-daysweek-walked-10-minutes_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'864_Number-of-daysweek-walked-10-minutes_0_0': 0},
    normalization={'mean': 5.369732285440756, 'std': 1.9564911925721618}, annotation_units=1,
)
walk_duration = TensorMap(
    '874_Duration-of-walks_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'874_Duration-of-walks_0_0': 0},
    normalization={'mean': 61.64092215093373, 'std': 78.79522990818906}, annotation_units=1,
)
physical_activities = TensorMap(
    '884_Number-of-daysweek-of-moderate-physical-activity-10-minutes_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'884_Number-of-daysweek-of-moderate-physical-activity-10-minutes_0_0': 0 },
    normalization={'mean': 3.6258833281089258, 'std': 2.3343738999823676}, annotation_units=1,
)
physical_activity = TensorMap(
    '894_Duration-of-moderate-activity_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'894_Duration-of-moderate-activity_0_0': 0 },
    normalization={'mean': 66.2862593866103, 'std': 77.28681218835422}, annotation_units=1,
)
physical_activity_vigorous = TensorMap(
    '904_Number-of-daysweek-of-vigorous-physical-activity-10-minutes_0_0', Interpretation.CONTINUOUS,
    channel_map={'904_Number-of-daysweek-of-vigorous-physical-activity-10-minutes_0_0': 0}, path_prefix='continuous',
    normalization={'mean': 1.838718301735063, 'std': 1.9593505421480895}, annotation_units=1,
)
physical_activity_vigorous_duration = TensorMap(
    '914_Duration-of-vigorous-activity_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'914_Duration-of-vigorous-activity_0_0': 0},
    normalization={'mean': 44.854488382965144, 'std': 48.159967071781466}, annotation_units=1,
)
tv = TensorMap(
    '1070_Time-spent-watching-television-TV_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'1070_Time-spent-watching-television-TV_0_0': 0},
    normalization={'mean': 2.7753595642790914, 'std': 1.7135478462887321}, annotation_units=1,
)
computer = TensorMap(
    '1080_Time-spent-using-computer_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'1080_Time-spent-using-computer_0_0': 0},
    normalization={'mean': 0.9781465855433753, 'std': 1.4444414103121512}, annotation_units=1,
)
car = TensorMap(
    '1090_Time-spent-driving_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'1090_Time-spent-driving_0_0': 0},
    normalization={'mean': 0.8219851505445748, 'std': 1.304094814200189}, annotation_units=1,
)
summer = TensorMap(
    '1050_Time-spend-outdoors-in-summer_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'1050_Time-spend-outdoors-in-summer_0_0': 0},
    normalization={'mean': 3.774492304870845, 'std': 2.430483731404539}, annotation_units=1,
)
winter = TensorMap(
    '1060_Time-spent-outdoors-in-winter_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'1060_Time-spent-outdoors-in-winter_0_0': 0},
    normalization={'mean': 1.8629686916635555, 'std': 1.88916218603397}, annotation_units=1,
)

systolic_blood_pressure_0 = TensorMap(
    '4080_Systolic-blood-pressure-automated-reading_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh',
    channel_map={'4080_Systolic-blood-pressure-automated-reading_0_0': 0}, validator=make_range_validator(40, 400),
    normalization={'mean': 137.79964191990328, 'std': 19.292863700283757},
)
diastolic_blood_pressure_0 = TensorMap(
    '4079_Diastolic-blood-pressure-automated-reading_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh',
    channel_map={'4079_Diastolic-blood-pressure-automated-reading_0_0': 0}, validator=make_range_validator(20, 300),
    normalization={'mean': 82.20657551284782, 'std': 10.496040770224475},
)

systolic_blood_pressure_1 = TensorMap(
    '4080_Systolic-blood-pressure-automated-reading_1_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh',
    channel_map={'4080_Systolic-blood-pressure-automated-reading_1_0': 0}, validator=make_range_validator(40, 400),
    normalization={'mean': 137.79964191990328, 'std': 19.292863700283757},
)
diastolic_blood_pressure_1 = TensorMap(
    '4079_Diastolic-blood-pressure-automated-reading_1_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh',
    channel_map={'4079_Diastolic-blood-pressure-automated-reading_1_0': 0}, validator=make_range_validator(20, 300),
    normalization={'mean': 82.20657551284782, 'std': 10.496040770224475},
)

systolic_blood_pressure_2 = TensorMap(
    '4080_Systolic-blood-pressure-automated-reading_2_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh',
    channel_map={'4080_Systolic-blood-pressure-automated-reading_2_0': 0}, validator=make_range_validator(40, 400),
    normalization={'mean': 137.79964191990328, 'std': 19.292863700283757},
)
diastolic_blood_pressure_2 = TensorMap(
    '4079_Diastolic-blood-pressure-automated-reading_2_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh',
    channel_map={'4079_Diastolic-blood-pressure-automated-reading_2_0': 0}, validator=make_range_validator(20, 300),
    normalization={'mean': 82.20657551284782, 'std': 10.496040770224475},
)


categorical_phenotypes_25 = TensorMap(
    'categorical-phenotypes-25', Interpretation.CATEGORICAL, path_prefix='categorical',
    channel_map={
        'Adopted-as-a-child_No_0_0': 0,
        'Beef-intake_Less-than-once-a-week_0_0': 1,
        'Breastfed-as-a-baby_Yes_0_0': 2,
        'Country-of-birth-UKelsewhere_England_0_0': 3,
        'Current-tobacco-smoking_No_0_0': 4,
        'Drive-faster-than-motorway-speed-limit_Neverrarely_0_0': 5,
        'Fracturedbroken-bones-in-last-5-years_No_0_0': 6,
        'Genetic-sex_Female_0_0': 7,
        'Hearing-difficultyproblems_No_0_0': 8,
        'Major-dietary-changes-in-the-last-5-years_No_0_0': 9,
        'Processed-meat-intake_Less-than-once-a-week_0_0': 10,
        'Processed-meat-intake_Once-a-week_0_0': 11,
        'Salt-added-to-food_Neverrarely_0_0': 12,
        'Salt-added-to-food_Sometimes_0_0': 13,
        'Shortness-of-breath-walking-on-level-ground_No_0_0': 14,
        'Smoking-status_Never_0_0': 15,
        'Smokingsmokers-in-household_No_0_0': 16,
        'Smoking-status_Previous_0_0': 17,
        'Usual-walking-pace_Brisk-pace_0_0': 18,
        'Usual-walking-pace_Steady-average-pace_0_0': 19,
        'Variation-in-diet_Sometimes_0_0': 20,
        'Vitamin-and-mineral-supplements_None-of-the-above_0_0': 21,
        'Wears-glasses-or-contact-lenses_Yes_0_0': 22,
        'Weight-change-compared-with-1-year-ago_No-weigh-about-the-same_0_0': 23,
        'Weight-change-compared-with-1-year-ago_Yes-gained-weight_0_0': 24,

    },
)

categorical_phenotypes_36 = TensorMap(
    'categorical-phenotypes-36', Interpretation.CATEGORICAL, path_prefix='categorical',
    channel_map={
        'Adopted-as-a-child_No_0_0': 0,
        'Breastfed-as-a-baby_Yes_0_0': 1,
        'Country-of-birth-UKelsewhere_England_0_0': 2,
        'Current-tobacco-smoking_No_0_0': 3,
        'Drive-faster-than-motorway-speed-limit_Neverrarely_0_0': 4,
        'Facial-ageing_Younger-than-you-are_0_0': 5,
        'Father-still-alive_No_0_0': 6,
        'Fracturedbroken-bones-in-last-5-years_No_0_0': 7,
        'Genetic-sex_Female_0_0': 8,
        'Handedness-chiralitylaterality_Righthanded_0_0': 9,
        'Hearing-difficultyproblems_No_0_0': 10,
        'Longstanding-illness-disability-or-infirmity_No_0_0': 11,
        'Maternal-smoking-around-birth_No_0_0': 12,
        'Milk-type-used_Semiskimmed_0_0': 13,
        'Mineral-and-other-dietary-supplements_None-of-the-above_0_0': 14,
        'Mother-still-alive_No_0_0': 15,
        'Pacemaker_No_0_0': 16,
        'Part-of-a-multiple-birth_No_0_0': 17,
        'Past-tobacco-smoking_I-have-never-smoked_0_0': 18,
        'Past-tobacco-smoking_Smoked-on-most-or-all-days_0_0': 19,
        'Pork-intake_Less-than-once-a-week_0_0': 20,
        'Pork-intake_Never_0_0': 21,
        'Poultry-intake_24-times-a-week_0_0': 22,
        'Poultry-intake_Once-a-week_0_0': 23,
        'Processed-meat-intake_Less-than-once-a-week_0_0': 24,
        'Processed-meat-intake_Once-a-week_0_0': 25,
        'Salt-added-to-food_Neverrarely_0_0': 26,
        'Salt-added-to-food_Sometimes_0_0': 27,
        'Smoking-status_Previous_0_0': 28,
        'Smokingsmokers-in-household_No_0_0': 29,
        'Types-of-physical-activity-in-last-4-weeks_Walking-for-pleasure-not-as-a-means-of-transport_0_0': 30,
        'Types-of-transport-used-excluding-work_Carmotor-vehicle_0_0': 31,
        'Vitamin-and-mineral-supplements_None-of-the-above_0_0': 32,
        'Wears-glasses-or-contact-lenses_Yes_0_0': 33,
        'Wheeze-or-whistling-in-the-chest-in-last-year_No_0_0': 35,
    },
)

categorical_phenotypes_78 = TensorMap(
    'categorical-phenotypes-78', Interpretation.CATEGORICAL, path_prefix='categorical', annotation_units=64,
    channel_map={
        'Adopted-as-a-child_No_0_0': 0,
        'Alcohol-intake-versus-10-years-previously_Less-nowadays_0_0': 1,
        'Alcohol-usually-taken-with-meals_It-varies_0_0': 2,
        'Alcohol-usually-taken-with-meals_Yes_0_0': 3,
        'Beef-intake_Less-than-once-a-week_0_0': 4,
        'Bread-type_Wholemeal-or-wholegrain_0_0': 5,
        'Breastfed-as-a-baby_Yes_0_0': 6,
        'Cereal-type_Other-eg-Cornflakes-Frosties_0_0': 7,
        'Cheese-intake_24-times-a-week_0_0': 8,
        'Coffee-type_Instant-coffee_0_0': 9,
        'Comparative-body-size-at-age-10_About-average_0_0': 10,
        'Comparative-body-size-at-age-10_Thinner_0_0': 11,
        'Comparative-height-size-at-age-10_About-average_0_0': 12,
        'Country-of-birth-UKelsewhere_England_0_0': 13,
        'Current-tobacco-smoking_No_0_0': 14,
        'Drive-faster-than-motorway-speed-limit_Neverrarely_0_0': 15,
        'Duration-walking-for-pleasure_Between-30-minutes-and-1-hour_0_0': 16,
        'Ease-of-skin-tanning_Get-moderately-tanned_0_0': 17,
        'FI1-numeric-addition-test_15_0_0': 18,
        'FI3-word-interpolation_Adult_0_0': 19,
        'FI4-positional-arithmetic_6_0_0': 20,
        'FI7-synonym_Cease_0_0': 21,
        'Facial-ageing_Younger-than-you-are_0_0': 21,
        'Father-still-alive_No_0_0': 22,
        'Father-still-alive_Yes_0_0': 23,
        'Fracturedbroken-bones-in-last-5-years_No_0_0': 24,
        'Frequency-of-stair-climbing-in-last-4-weeks_610-times-a-day_0_0': 25,
        'Genetic-sex_Female_0_0': 26,
        'Had-menopause_Yes_0_0': 27,
        'Hair-colour-natural-before-greying_Dark-brown_0_0': 28,
        'Hair-colour-natural-before-greying_Light-brown_0_0': 29,
        'Handedness-chiralitylaterality_Righthanded_0_0': 30,
        'Hearing-difficultyproblems_No_0_0': 31,
        'Hot-drink-temperature_Hot_0_0': 32,
        'Hot-drink-temperature_Very-hot_0_0': 33,
        'Lambmutton-intake_Less-than-once-a-week_0_0': 34,
        'Lambmutton-intake_Never_0_0': 35,
        'Major-dietary-changes-in-the-last-5-years_No_0_0': 36,
        'Major-dietary-changes-in-the-last-5-years_Yes-because-of-other-reasons_0_0': 37,
        'Maternal-smoking-around-birth_No_0_0': 38,
        'Milk-type-used_Semiskimmed_0_0': 39,
        'Mineral-and-other-dietary-supplements_None-of-the-above_0_0': 40,
        'Mother-still-alive_No_0_0': 41,
        'Mother-still-alive_Yes_0_0': 42,
        'Mouthteeth-dental-problems_None-of-the-above_0_0': 43,
        'Nonoily-fish-intake_Once-a-week_0_0': 44,
        'Oily-fish-intake_Less-than-once-a-week_0_0': 45,
        'Oily-fish-intake_Once-a-week_0_0': 46,
        'Pain-types-experienced-in-last-month_Headache_0_0': 47,
        'Pain-types-experienced-in-last-month_None-of-the-above_0_0': 48,
        'Part-of-a-multiple-birth_No_0_0': 49,
        'Past-tobacco-smoking_I-have-never-smoked_0_0': 50,
        'Past-tobacco-smoking_Smoked-on-most-or-all-days_0_0': 51,
        'Pork-intake_Less-than-once-a-week_0_0': 52,
        'Pork-intake_Never_0_0': 53,
        'Poultry-intake_24-times-a-week_0_0': 54,
        'Poultry-intake_Once-a-week_0_0': 55,
        'Processed-meat-intake_Less-than-once-a-week_0_0': 56,
        'Processed-meat-intake_Once-a-week_0_0': 57,
        'Salt-added-to-food_Neverrarely_0_0': 58,
        'Salt-added-to-food_Sometimes_0_0': 59,
        'Shortness-of-breath-walking-on-level-ground_No_0_0': 60,
        'Skin-colour_Fair_0_0': 61,
        'Smoking-status_Never_0_0': 62,
        'Smoking-status_Previous_0_0': 63,
        'Smokingsmokers-in-household_No_0_0': 64,
        'Spread-type_Other-type-of-spreadmargarine_0_0': 65,
        'Types-of-physical-activity-in-last-4-weeks_Other-exercises-eg-swimming-cycling-keep-fit-bowling_0_1': 66,
        'Types-of-physical-activity-in-last-4-weeks_Walking-for-pleasure-not-as-a-means-of-transport_0_0': 67,
        'Types-of-transport-used-excluding-work_Carmotor-vehicle_0_0': 68,
        'Types-of-transport-used-excluding-work_Walk_0_1': 69,
        'Usual-walking-pace_Brisk-pace_0_0': 70,
        'Usual-walking-pace_Steady-average-pace_0_0': 71,
        'Variation-in-diet_Sometimes_0_0': 72,
        'Vitamin-and-mineral-supplements_None-of-the-above_0_0': 73,
        'Wears-glasses-or-contact-lenses_Yes_0_0': 74,
        'Weight-change-compared-with-1-year-ago_No-weigh-about-the-same_0_0': 75,
        'Weight-change-compared-with-1-year-ago_Yes-gained-weight_0_0': 76,
        'Wheeze-or-whistling-in-the-chest-in-last-year_No_0_0': 77,
    },
)

categorical_phenotypes_134 = TensorMap(
        'categorical-phenotypes-134', Interpretation.CATEGORICAL, path_prefix='categorical', annotation_units=64,
        channel_map={
            'Adopted-as-a-child_No_0_0': 0, 'Adopted-as-a-child_No_2_0': 1,
            'Alcohol-intake-frequency_Once-or-twice-a-week_0_0': 2,
            'Alcohol-intake-versus-10-years-previously_About-the-same_0_0': 3,
            'Alcohol-intake-versus-10-years-previously_Less-nowadays_2_0': 4,
            'Alcohol-usually-taken-with-meals_Yes_2_0': 5, 'Beef-intake_Less-than-once-a-week_0_0': 6,
            'Beef-intake_Less-than-once-a-week_2_0': 7,
            'Bread-type_Wholemeal-or-wholegrain_0_0': 10,
            'Bread-type_Wholemeal-or-wholegrain_2_0': 11,
            'Breastfed-as-a-baby_Yes_0_0': 12, 'Breathing-problems-during-period-of-job_No_0_0': 13,
            'Breathing-problems-during-period-of-job_No_0_1': 14,
            'Cereal-type_Oat-cereal-eg-Ready-Brek-porridge_0_0': 17,
            'Cheese-intake_24-times-a-week_0_0': 18, 'Cheese-intake_24-times-a-week_2_0': 19,
            'Coffee-type_Instant-coffee_0_0': 20, 'Comparative-body-size-at-age-10_Thinner_0_0': 21,
            'Comparative-height-size-at-age-10_About-average_0_0': 22,
            'Country-of-birth-UKelsewhere_England_0_0': 23, 'Country-of-birth-UKelsewhere_England_2_0': 24,
            'Current-tobacco-smoking_No_0_0': 25,
            'Current-tobacco-smoking_No_2_0': 26,
            'Drive-faster-than-motorway-speed-limit_Neverrarely_2_0': 29,
            'Drive-faster-than-motorway-speed-limit_Sometimes_0_0': 30,
            'Ease-of-skin-tanning_Get-moderately-tanned_0_0': 31,
            'Ever-had-breast-cancer-screening-mammogram_Yes_0_0': 32,
            'Ever-had-breast-cancer-screening-mammogram_Yes_2_0': 33,
            'Ever-had-cervical-smear-test_Yes_0_0': 34, 'Ever-had-cervical-smear-test_Yes_2_0': 35,
            'Ever-taken-oral-contraceptive-pill_Yes_0_0': 36,
            'Ever-taken-oral-contraceptive-pill_Yes_2_0': 37,
            'Ever-used-hormonereplacement-therapy-HRT_No_0_0': 38,
            'Eye-problemsdisorders_None-of-the-above_0_0': 39,
            'Eye-problemsdisorders_None-of-the-above_2_0': 40, 'FI1-numeric-addition-test_15_2_0': 41,
            'FI3-word-interpolation_Adult_2_0': 42, 'FI4-positional-arithmetic_6_2_0': 43,
            'FI6-conditional-arithmetic_69_2_0': 44, 'FI7-synonym_Cease_2_0': 45,
            'Facial-ageing_Younger-than-you-are_0_0': 46, 'Facial-ageing_Younger-than-you-are_2_0': 47,
            'Father-still-alive_No_0_0': 48, 'Father-still-alive_No_2_0': 49,
            'Father-still-alive_Yes_0_0': 50, 'Fracturedbroken-bones-in-last-5-years_No_0_0': 51,
            'Fracturedbroken-bones-in-last-5-years_No_2_0': 52, 'Genetic-sex_Female_0_0': 53,
            'Had-menopause_Yes_2_0': 54, 'Hair-colour-natural-before-greying_Dark-brown_0_0': 55,
            'Handedness-chiralitylaterality_Righthanded_0_0': 56,
            'Handedness-chiralitylaterality_Righthanded_2_0': 57, 'Hearing-difficultyproblems_No_0_0': 58,
            'Hearing-difficultyproblems_No_2_0': 59, 'Hot-drink-temperature_Hot_0_0': 60,
            'Hot-drink-temperature_Hot_2_0': 61, 'Lambmutton-intake_Less-than-once-a-week_0_0': 62,
            'Lambmutton-intake_Less-than-once-a-week_2_0': 63,
            'Longstanding-illness-disability-or-infirmity_No_0_0': 64,
            'Longstanding-illness-disability-or-infirmity_No_2_0': 65,
            'Major-dietary-changes-in-the-last-5-years_No_0_0': 66,
            'Major-dietary-changes-in-the-last-5-years_No_2_0': 67,
            'Maternal-smoking-around-birth_No_0_0': 68,
            'Milk-type-used_Semiskimmed_0_0': 73,
            'Milk-type-used_Semiskimmed_2_0': 74,
            'Milk-type-used_Skimmed_0_0': 75,
            'Mineral-and-other-dietary-supplements_None-of-the-above_0_0': 76,
            'Mineral-and-other-dietary-supplements_None-of-the-above_2_0': 77,
            'Mother-still-alive_No_2_0': 78,
            'Mother-still-alive_Yes_0_0': 79,
            'Mouthteeth-dental-problems_None-of-the-above_0_0': 80,
            'Mouthteeth-dental-problems_None-of-the-above_2_0': 81, 'Noisy-workplace_No_2_0': 82,
            'Nonoily-fish-intake_Once-a-week_0_0': 83, 'Nonoily-fish-intake_Once-a-week_2_0': 84,
            'Oily-fish-intake_Less-than-once-a-week_0_0': 85,
            'Overall-health-rating_Good_0_0': 88,
            'Overall-health-rating_Good_2_0': 89,
            'Pain-types-experienced-in-last-month_None-of-the-above_0_0': 92,
            'Part-of-a-multiple-birth_No_0_0': 93, 'Part-of-a-multiple-birth_No_2_0': 94,
            'Past-tobacco-smoking_I-have-never-smoked_0_0': 95,
            'Past-tobacco-smoking_I-have-never-smoked_2_0': 96, 'Pork-intake_Less-than-once-a-week_0_0': 97,
            'Pork-intake_Less-than-once-a-week_2_0': 98, 'Poultry-intake_24-times-a-week_0_0': 99,
            'Poultry-intake_24-times-a-week_2_0': 100, 'Processed-meat-intake_Less-than-once-a-week_2_0': 101,
            'Salt-added-to-food_Neverrarely_0_0': 102, 'Salt-added-to-food_Neverrarely_2_0': 103,
            'Shortness-of-breath-walking-on-level-ground_No_2_0': 104, 'Skin-colour_Fair_0_0': 105,
            'Smoking-status_Never_0_0': 106, 'Smoking-status_Never_2_0': 107,
            'Smokingsmokers-in-household_No_0_0': 108, 'Smokingsmokers-in-household_No_2_0': 109,
            'Spread-type_Butterspreadable-butter_2_0': 110,
            'Spread-type_Other-type-of-spreadmargarine_0_0': 111,
            'Types-of-physical-activity-in-last-4-weeks_Other-exercises-eg-swimming-cycling-keep-fit-bowling_0_1': 112,
            'Types-of-physical-activity-in-last-4-weeks_Other-exercises-eg-swimming-cycling-keep-fit-bowling_2_1': 113,
            'Types-of-physical-activity-in-last-4-weeks_Walking-for-pleasure-not-as-a-means-of-transport_0_0': 114,
            'Types-of-physical-activity-in-last-4-weeks_Walking-for-pleasure-not-as-a-means-of-transport_2_0': 115,
            'Types-of-transport-used-excluding-work_Carmotor-vehicle_0_0': 116,
            'Types-of-transport-used-excluding-work_Carmotor-vehicle_2_0': 117,
            'Types-of-transport-used-excluding-work_Walk_0_1': 118,
            'Types-of-transport-used-excluding-work_Walk_2_1': 119,
            'UK-Biobank-assessment-centre_Cheadle-imaging_2_0': 120, 'Usual-walking-pace_Brisk-pace_0_0': 121,
            'Usual-walking-pace_Brisk-pace_2_0': 122, 'Usual-walking-pace_Steady-average-pace_0_0': 123,
            'Usual-walking-pace_Steady-average-pace_2_0': 124, 'Variation-in-diet_Sometimes_0_0': 125,
            'Variation-in-diet_Sometimes_2_0': 126,
            'Vitamin-and-mineral-supplements_None-of-the-above_0_0': 129,
            'Vitamin-and-mineral-supplements_None-of-the-above_2_0': 130,
            'Wears-glasses-or-contact-lenses_Yes_0_0': 131, 'Wears-glasses-or-contact-lenses_Yes_2_0': 132,
            'Weight-change-compared-with-1-year-ago_No-weigh-about-the-same_0_0': 133,
            'Weight-change-compared-with-1-year-ago_No-weigh-about-the-same_2_0': 8,
            'Wheeze-or-whistling-in-the-chest-in-last-year_No_0_0': 9,
            'Wheeze-or-whistling-in-the-chest-in-last-year_No_2_0': 15,
            'Worked-with-materials-containing-asbestos_Rarelynever_0_0': 16,
            'Worked-with-materials-containing-asbestos_Rarelynever_0_1': 27,
            'Worked-with-paints-thinners-or-glues_Rarelynever_0_0': 28,
            'Worked-with-paints-thinners-or-glues_Rarelynever_0_1': 69,
            'Worked-with-pesticides_Rarelynever_0_0': 70, 'Worked-with-pesticides_Rarelynever_0_1': 71,
            'Workplace-full-of-chemical-or-other-fumes_Rarelynever_0_0': 72,
            'Workplace-full-of-chemical-or-other-fumes_Rarelynever_0_1': 86,
            'Workplace-had-a-lot-of-cigarette-smoke-from-other-people-smoking_Rarelynever_0_0': 87,
            'Workplace-had-a-lot-of-diesel-exhaust_Rarelynever_0_0': 90,
            'Workplace-had-a-lot-of-diesel-exhaust_Rarelynever_0_1': 91,
            'Workplace-very-dusty_Rarelynever_0_0': 127, 'Workplace-very-dusty_Rarelynever_0_1': 128,
        },
)
