from ml4cvd.tensor_from_file import normalized_first_date, TMAPS
from ml4cvd.TensorMap import TensorMap, make_range_validator, Interpretation
from ml4cvd.defines import MRI_SEGMENTED_CHANNEL_MAP, ECG_CHAR_2_IDX, StorageType
from ml4cvd.metrics import weighted_crossentropy, ignore_zeros_logcosh, y_true_times_mse, y_true_squared_times_mse, y_true_cubed_times_mse, y_true_squared_times_logcosh


diploid_cm = {'homozygous_reference': 0, 'heterozygous': 1, 'homozygous_variant': 2}
TMAPS['rs3829740'] = TensorMap('rs3829740', Interpretation.CATEGORICAL, channel_map=diploid_cm)
TMAPS['rs2234962'] = TensorMap('rs2234962', Interpretation.CATEGORICAL, channel_map=diploid_cm)
TMAPS['rs2042995'] = TensorMap('rs2042995', Interpretation.CATEGORICAL, channel_map=diploid_cm)

TMAPS['rs3829740_weighted'] = TensorMap('rs3829740', Interpretation.CATEGORICAL, channel_map=diploid_cm, loss=weighted_crossentropy([1, 1, 1.5], 'rs3829740'))
TMAPS['rs2234962_weighted'] = TensorMap('rs2234962', Interpretation.CATEGORICAL, channel_map=diploid_cm, loss=weighted_crossentropy([.8, 1, 1.5], 'rs2234962'))
TMAPS['rs2042995_weighted'] = TensorMap('rs2042995', Interpretation.CATEGORICAL, channel_map=diploid_cm, loss=weighted_crossentropy([.6, 1.5, 2], 'rs2042995'))


TMAPS['akap9_lof'] = TensorMap('AKAP9', Interpretation.CATEGORICAL, channel_map={'no_akap9_lof': 0, 'akap9_lof': 1})
TMAPS['dsc2_lof'] = TensorMap('DSC2', Interpretation.CATEGORICAL, channel_map={'no_dsc2_lof': 0, 'dsc2_lof': 1})
TMAPS['ryr2_lof'] = TensorMap('RYR2', Interpretation.CATEGORICAL, channel_map={'no_ryr2_lof': 0, 'ryr2_lof': 1})
TMAPS['ttn_lof'] = TensorMap('TTN', Interpretation.CATEGORICAL, channel_map={'no_ttn_lof': 0, 'ttn_lof': 1})


TMAPS['ecg_semi_coarse'] = TensorMap(
    'ecg_semi_coarse', Interpretation.CATEGORICAL, loss=weighted_crossentropy([1.0, 1.0, 2.0, 4.0, 16.0, 20.0], 'ecg_semi_coarse'),
    channel_map={'Normal_sinus_rhythm': 0, 'Sinus_bradycardia': 1, 'Marked_sinus_bradycardia': 2, 'Other_sinus_rhythm': 3, 'Atrial_fibrillation': 4, 'Other_rhythm': 5},
)


TMAPS['ecg_semi_coarse_with_poor'] = TensorMap(
    'ecg_semi_coarse_with_poor', Interpretation.CATEGORICAL, loss=weighted_crossentropy([1.0, 2.0, 3.0, 3.0, 20.0, 20.0], 'ecg_semi_coarse_with_poor'),
    channel_map={'Normal_sinus_rhythm': 0, 'Sinus_bradycardia': 1, 'Marked_sinus_bradycardia': 2, 'Other_sinus_rhythm': 3, 'Atrial_fibrillation': 4, 'Other_rhythm': 5},
)

TMAPS['ecg_normal'] = TensorMap(
    'ecg_normal', Interpretation.CATEGORICAL, loss=weighted_crossentropy([2.0, 3.0, 3.0, 3.0], 'ecg_normal'),
    channel_map={'Normal_ECG': 0, 'Abnormal_ECG': 1, 'Borderline_ECG': 2, 'Otherwise_normal_ECG': 3},
)
TMAPS['ecg_infarct'] = TensorMap(
    'ecg_infarct', Interpretation.CATEGORICAL, channel_map={'no_infarct': 0, 'infarct': 1},
    loss=weighted_crossentropy([1.0, 8.0], 'ecg_infarct'),
)
TMAPS['ecg_poor_data'] = TensorMap(
    'ecg_poor_data', Interpretation.CATEGORICAL, channel_map={'no_poor_data_quality': 0, 'poor_data_quality': 1},
    loss=weighted_crossentropy([1.0, 8.0], 'ecg_poor_data'),
)
TMAPS['ecg_block'] = TensorMap(
    'ecg_block', Interpretation.CATEGORICAL, channel_map={'no_block': 0, 'block': 1},
    loss=weighted_crossentropy([1.0, 8.0], 'ecg_block'),
)

TMAPS['ecg_rest_next_char'] = TensorMap('ecg_rest_next_char', Interpretation.LANGUAGE, shape=(len(ECG_CHAR_2_IDX),), channel_map=ECG_CHAR_2_IDX, activation='softmax', loss='categorical_crossentropy', loss_weight=2.0)
TMAPS['ecg_rest_text'] = TensorMap('ecg_rest_text', Interpretation.LANGUAGE, shape=(100, len(ECG_CHAR_2_IDX)), path_prefix='ukb_ecg_rest', channel_map={'context': 0, 'alphabet': 1}, dependent_map=TMAPS['ecg_rest_next_char'])

TMAPS['p-axis'] = TensorMap(
    'PAxis', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'PAxis': 0}, loss='logcosh', validator=make_range_validator(-50, 130),
    normalization={'mean': 48.7, 'std': 23.1},
)
TMAPS['p-duration'] = TensorMap(
    'PDuration', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'PDuration': 0}, loss='logcosh', validator=make_range_validator(30, 140),
    normalization={'mean': 96.1, 'std': 18.85},
)
TMAPS['p-offset'] = TensorMap(
    'POffset', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'POffset': 0}, loss='logcosh', validator=make_range_validator(200, 500),
    normalization={'mean': 369.1, 'std': 28.42},
)
TMAPS['p-onset'] = TensorMap(
    'POnset', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'POnset': 0}, loss='logcosh', validator=make_range_validator(120, 400),
    normalization={'mean': 275.1, 'std': 26.420},
)
TMAPS['pp-interval'] = TensorMap(
    'PPInterval', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'PPInterval': 0}, loss='logcosh', validator=make_range_validator(300, 1800),
    normalization={'mean': 1036.1, 'std': 185.0},
)
TMAPS['pq-interval'] = TensorMap(
    'PQInterval', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'PQInterval': 0}, loss='logcosh', validator=make_range_validator(70, 400),
    normalization={'mean': 165.9, 'std': 26.3},
)
TMAPS['q-offset'] = TensorMap(
    'QOffset', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'QOffset': 0}, loss='logcosh', validator=make_range_validator(300, 600),
    normalization={'mean': 525.1, 'std': 13.52},
)
TMAPS['q-onset'] = TensorMap(
    'QOnset', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'QOnset': 0}, loss='logcosh', validator=make_range_validator(370, 600),
    normalization={'mean': 435.1, 'std': 11.420},
)
TMAPS['qrs-complexes'] = TensorMap(
    'QRSComplexes', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'QRSComplexes': 0}, loss='logcosh', validator=make_range_validator(0, 60),
    normalization={'mean': 8.0, 'std': 20.0},
)
TMAPS['qrs-duration'] = TensorMap(
    'QRSDuration', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'QRSDuration': 0}, loss='logcosh', validator=make_range_validator(45, 175),
    normalization={'mean': 89.53, 'std': 12.21},
)
TMAPS['qrs-num'] = TensorMap(
    'QRSNum', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'QRSNum': 0}, loss='logcosh', validator=make_range_validator(2, 30),
    normalization={'mean': 9.61, 'std': 1.64},
)
TMAPS['qt-interval'] = TensorMap(
    'QTInterval', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'QTInterval': 0}, loss='logcosh', validator=make_range_validator(300, 600),
    normalization={'mean': 426.1, 'std': 32.24},
)
TMAPS['qt-interval-quintiles'] = TensorMap(
    'QTInterval', Interpretation.DISCRETIZED, path_prefix='ukb_ecg_rest',
    channel_map={'QTInterval': 0}, normalization={'mean': 426.1, 'std': 32.24},
    discretization_bounds=[-0.842, -0.253, 0.253, 0.842],
)
TMAPS['qtc-interval'] = TensorMap(
    'QTCInterval', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'QTCInterval': 0}, loss='logcosh', validator=make_range_validator(300, 600),
    normalization={'mean': 419.1, 'std': 20.7},
)
TMAPS['r-axis'] = TensorMap(
    'RAxis', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'RAxis': 0}, loss='logcosh', validator=make_range_validator(-100, 200),
    normalization={'mean': 25.7, 'std': 36.6},
)
TMAPS['rr-interval'] = TensorMap(
    'RRInterval', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'RRInterval': 0}, loss='logcosh', validator=make_range_validator(400, 2000),
    normalization={'mean': 1040.61, 'std': 175.5},
)
TMAPS['ventricular-rate'] = TensorMap(
    'VentricularRate', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'VentricularRate': 0}, validator=make_range_validator(30, 150),
    loss='logcosh', normalization={'mean': 59.3, 'std': 10.6},
)
TMAPS['t-offset'] = TensorMap(
    'TOffset', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'TOffset': 0}, loss='logcosh', validator=make_range_validator(700, 1000),
    normalization={'mean': 860.7, 'std': 32.52},
)
TMAPS['t-axis'] = TensorMap(
    'TAxis', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'TAxis': 0}, loss='logcosh', validator=make_range_validator(-100, 200),
    normalization={'mean': 40.8, 'std': 32.6},
)

TMAPS['af_prs'] = TensorMap('AF_PRS_LDscore', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'AF_PRS_LDscore': 0}, normalization={'mean': -1.0, 'std': 0.4})
TMAPS['charge'] = TensorMap(
    'charge', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'charge': 0}, normalization={'mean': 12.0, 'std': 2.0},
    validator=make_range_validator(0, 20),
)

TMAPS['qtc-intervalp'] = TensorMap(
    'QTCInterval', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'QTCInterval': 0}, loss='logcosh', validator=make_range_validator(100, 900),
    parents=[TMAPS['qt-interval'], TMAPS['rr-interval']], normalization={'mean': 419.1, 'std': 20.7},
)
TMAPS['qrs-durationpp'] = TensorMap(
    'QRSDuration', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'QRSDuration': 0}, loss='logcosh', validator=make_range_validator(45, 175),
    normalization={'mean': 89.53, 'std': 12.21},
    parents=[TMAPS['qtc-intervalp']],
)

TMAPS['p-axis-sentinel'] = TensorMap(
    'PAxis', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'PAxis': 0}, sentinel=0, metrics=['logcosh'],
    normalization={'mean': 48.7, 'std': 23.1},
)
TMAPS['p-duration-sentinel'] = TensorMap(
    'PDuration', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'PDuration': 0}, sentinel=0, metrics=['logcosh'],
    normalization={'mean': 96.1, 'std': 18.85},
)
TMAPS['p-offset-sentinel'] = TensorMap(
    'POffset', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'POffset': 0}, sentinel=0, metrics=['logcosh'],
    normalization={'mean': 369.1, 'std': 28.42},
)
TMAPS['p-onset-sentinel'] = TensorMap(
    'POnset', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'POnset': 0}, sentinel=0, metrics=['logcosh'],
    normalization={'mean': 275.1, 'std': 26.420},
)
TMAPS['pp-interval-sentinel'] = TensorMap(
    'PPInterval', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'PPInterval': 0}, sentinel=0, metrics=['logcosh'],
    normalization={'mean': 1036.1, 'std': 185.0},
)
TMAPS['pq-interval-sentinel'] = TensorMap(
    'PQInterval', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'PQInterval': 0}, sentinel=0, metrics=['logcosh'],
    normalization={'mean': 165.9, 'std': 26.3},
)
TMAPS['qrs-duration-sentinel'] = TensorMap(
    'QRSDuration', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'QRSDuration': 0}, sentinel=0,
    normalization={'mean': 89.53, 'std': 12.21},
)
TMAPS['qt-interval-sentinel'] = TensorMap(
    'QTInterval', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'QTInterval': 0}, sentinel=0,
    normalization={'mean': 426.1, 'std': 32.24},
)
TMAPS['qtc-interval-sentinel'] = TensorMap(
    'QTCInterval', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'QTCInterval': 0}, sentinel=0,
    normalization={'mean': 419.1, 'std': 20.7},
)
TMAPS['qtc-intervalp-sentinel'] = TensorMap(
    'QTCInterval', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'QTCInterval': 0}, sentinel=0,
    normalization={'mean': 419.1, 'std': 20.7},
    parents=[TMAPS['qt-interval'], TMAPS['rr-interval']],
)
TMAPS['qtc-intervalp-sentinel'] = TensorMap(
    'QTCInterval', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'QTCInterval': 0}, sentinel=0,
    normalization={'mean': 419.1, 'std': 20.7},
    parents=[TMAPS['qt-interval'], TMAPS['rr-interval']],
)
TMAPS['r-axis-sentinel'] = TensorMap('RAxis', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'RAxis': 0}, sentinel=0, normalization={'mean': 25.7, 'std': 36.6})
TMAPS['rr-interval-sentinel'] = TensorMap(
    'RRInterval', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'RRInterval': 0}, sentinel=0,
    normalization={'mean': 1040.61, 'std': 175.5},
)
TMAPS['t-axis-sentinel'] = TensorMap('TAxis', Interpretation.CONTINUOUS, path_prefix='ukb_ecg_rest', channel_map={'TAxis': 0}, sentinel=0, normalization={'mean': 40.8, 'std': 32.6})


TMAPS['bb_baseline'] = TensorMap(
    'bb_baseline', Interpretation.CATEGORICAL, channel_map={'no_bb_baseline': 0, 'bb_baseline': 1},
    loss=weighted_crossentropy([0.0453, 0.9547], 'bb_baseline'),
)
TMAPS['ccb_baseline'] = TensorMap(
    'ccb_baseline', Interpretation.CATEGORICAL, channel_map={'no_ccb_baseline': 0, 'ccb_baseline': 1},
    loss=weighted_crossentropy([0.0044, 0.9956], 'ccb_baseline'),
)
TMAPS['class1_baseline'] = TensorMap(
    'class1_baseline', Interpretation.CATEGORICAL, channel_map={'no_class1_baseline': 0, 'class1_baseline': 1},
    loss=weighted_crossentropy([0.0023, 0.9977], 'class1_baseline'),
)
TMAPS['class3_baseline'] = TensorMap(
    'class3_baseline', Interpretation.CATEGORICAL, channel_map={'no_class3_baseline': 0, 'class3_baseline': 1},
    loss=weighted_crossentropy([0.0011, 0.9989], 'class3_baseline'),
)
TMAPS['qtc_drug_def_baseline'] = TensorMap(
    'qtc_drug_def_baseline', Interpretation.CATEGORICAL,
    channel_map={'no_qtc_drug_def_baseline': 0, 'qtc_drug_def_baseline': 1},
    loss=weighted_crossentropy([0.0210, 0.9790], 'qtc_drug_def_baseline'),
)
TMAPS['qtc_drug_poss_baseline'] = TensorMap(
    'qtc_drug_poss_baseline', Interpretation.CATEGORICAL,
    channel_map={'no_qtc_drug_poss_baseline': 0, 'qtc_drug_poss_baseline': 1},
    loss=weighted_crossentropy([0.0189, 0.9811], 'qtc_drug_poss_baseline'),
)
TMAPS['combined_qtc_drug_baseline'] = TensorMap(
    'combined_qtc_drug_baseline', Interpretation.CATEGORICAL,
    channel_map={'no_combined_qtc_drug_baseline': 0, 'combined_qtc_drug_baseline': 1},
    loss=weighted_crossentropy([0.0389, 0.9611], 'combined_qtc_drug_baseline'),
)

TMAPS['class1_baseline'] = TensorMap('class1_baseline', Interpretation.CATEGORICAL, channel_map={'no_class1_baseline': 0, 'class1_baseline': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0023, 0.9977], 'class1_baseline'))
TMAPS['bb_baseline'] = TensorMap('bb_baseline', Interpretation.CATEGORICAL, channel_map={'no_bb_baseline': 0, 'bb_baseline': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0453, 0.9547], 'bb_baseline'))
TMAPS['class3_baseline'] = TensorMap('class3_baseline', Interpretation.CATEGORICAL, channel_map={'no_class3_baseline': 0, 'class3_baseline': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0011, 0.9989], 'class3_baseline'))
TMAPS['ccb_baseline'] = TensorMap('ccb_baseline', Interpretation.CATEGORICAL, channel_map={'no_ccb_baseline': 0, 'ccb_baseline': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0044, 0.9956], 'ccb_baseline'))
TMAPS['qtc_drug_def_baseline'] = TensorMap('qtc_drug_def_baseline', Interpretation.CATEGORICAL, channel_map={'no_qtc_drug_def_baseline': 0, 'qtc_drug_def_baseline': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0210, 0.9790], 'qtc_drug_def_baseline'))
TMAPS['qtc_drug_poss_baseline'] = TensorMap('qtc_drug_poss_baseline', Interpretation.CATEGORICAL, channel_map={'no_qtc_drug_poss_baseline': 0, 'qtc_drug_poss_baseline': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0189, 0.9811], 'qtc_drug_poss_baseline'))
TMAPS['class1_fu'] = TensorMap('class1_fu', Interpretation.CATEGORICAL, channel_map={'no_class1_fu': 0, 'class1_fu': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0018, 0.9982], 'class1_fu'))
TMAPS['bb_fu'] = TensorMap('bb_fu', Interpretation.CATEGORICAL, channel_map={'no_bb_fu': 0, 'bb_fu': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0306, 0.9694], 'bb_fu'))
TMAPS['class3_fu'] = TensorMap('class3_fu', Interpretation.CATEGORICAL, channel_map={'no_class3_fu': 0, 'class3_fu': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0006, 0.9994], 'class3_fu'))
TMAPS['ccb_fu'] = TensorMap('ccb_fu', Interpretation.CATEGORICAL, channel_map={'no_ccb_fu': 0, 'ccb_fu': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0035, 0.9965], 'ccb_fu'))
TMAPS['qtc_drug_def_fu'] = TensorMap('qtc_drug_def_fu', Interpretation.CATEGORICAL, channel_map={'no_qtc_drug_def_fu': 0, 'qtc_drug_def_fu': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0140, 0.9860], 'qtc_drug_def_fu'))
TMAPS['qtc_drug_poss_fu'] = TensorMap('qtc_drug_poss_fu', Interpretation.CATEGORICAL, channel_map={'no_qtc_drug_poss_fu': 0, 'qtc_drug_poss_fu': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0127, 0.9873], 'qtc_drug_poss_fu'))
TMAPS['qtc_drug_def_any'] = TensorMap('qtc_drug_def_any', Interpretation.CATEGORICAL, channel_map={'no_qtc_drug_def_any': 0, 'qtc_drug_def_any': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0302, 0.9698], 'qtc_drug_def_any'))
TMAPS['qtc_drug_poss_any'] = TensorMap('qtc_drug_poss_any', Interpretation.CATEGORICAL, channel_map={'no_qtc_drug_poss_any': 0, 'qtc_drug_poss_any': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0267, 0.9733], 'qtc_drug_poss_any'))
TMAPS['any_class1'] = TensorMap('any_class1', Interpretation.CATEGORICAL, channel_map={'no_any_class1': 0, 'any_class1': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0031, 0.9969], 'any_class1'))
TMAPS['any_bb'] = TensorMap('any_bb', Interpretation.CATEGORICAL, channel_map={'no_any_bb': 0, 'any_bb': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0602, 0.9398], 'any_bb'))
TMAPS['any_class3'] = TensorMap('any_class3', Interpretation.CATEGORICAL, channel_map={'no_any_class3': 0, 'any_class3': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0013, 0.9987], 'any_class3'))
TMAPS['any_ccb'] = TensorMap('any_ccb', Interpretation.CATEGORICAL, channel_map={'no_any_ccb': 0, 'any_ccb': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0062, 0.9938], 'any_ccb'))
TMAPS['combined_qtc_drug_baseline'] = TensorMap('combined_qtc_drug_baseline', Interpretation.CATEGORICAL, channel_map={'no_combined_qtc_drug_baseline': 0, 'combined_qtc_drug_baseline': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0389, 0.9611], 'combined_qtc_drug_baseline'))
TMAPS['combined_qtc_drug_fu'] = TensorMap('combined_qtc_drug_fu', Interpretation.CATEGORICAL, channel_map={'no_combined_qtc_drug_fu': 0, 'combined_qtc_drug_fu': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0260, 0.9740], 'combined_qtc_drug_fu'))
TMAPS['combined_qtc_drug_any'] = TensorMap('combined_qtc_drug_any', Interpretation.CATEGORICAL, channel_map={'no_combined_qtc_drug_any': 0, 'combined_qtc_drug_any': 1}, loss_weight=100.0, loss=weighted_crossentropy([0.0546, 0.9454], 'combined_qtc_drug_any'))

TMAPS['ecg-bike-max-hr-no0'] = TensorMap(
    'bike_max_hr', Interpretation.CONTINUOUS, channel_map={'bike_max_hr': 0},
    loss=ignore_zeros_logcosh, metrics=['logcosh'], normalization={'mean': 110.03, 'std': 20.04},
)
TMAPS['ecg-bike-resting-hr-no0'] = TensorMap(
    'bike_resting_hr', Interpretation.CONTINUOUS, channel_map={'bike_resting_hr': 0},
    loss=ignore_zeros_logcosh, metrics=['logcosh'], normalization={'mean': 71.2, 'std': 12.57},
)
TMAPS['ecg-bike-max-pred-hr-no0'] = TensorMap(
    'bike_max_pred_hr', Interpretation.CONTINUOUS, channel_map={'bike_max_pred_hr': 0},
    loss=ignore_zeros_logcosh, metrics=['logcosh'], normalization={'mean': 167.5, 'std': 5.78},
)

TMAPS['weight_kg'] = TensorMap('weight_kg',  Interpretation.CONTINUOUS, normalization={'mean': 76.54286701805927, 'std': 15.467605416933122}, loss='logcosh', channel_map={'weight_kg': 0})
TMAPS['height_cm'] = TensorMap('height_cm',  Interpretation.CONTINUOUS, normalization={'mean': 169.18064748408653, 'std': 9.265265197273026}, loss='logcosh', channel_map={'height_cm': 0})
TMAPS['bmi_bsa'] = TensorMap('bmi',  Interpretation.CONTINUOUS, normalization={'mean': 26.65499238706321, 'std': 4.512077188749083}, loss='logcosh', channel_map={'bmi': 0})
TMAPS['bsa_mosteller'] = TensorMap('bsa_mosteller',  Interpretation.CONTINUOUS, normalization={'mean': 1.8894831981880114, 'std': 0.22169301057810176}, loss='logcosh', channel_map={'bsa_mosteller': 0})
TMAPS['bsa_dubois'] = TensorMap('bsa_dubois',  Interpretation.CONTINUOUS, normalization={'mean': 1.8671809970639703, 'std': 0.20913930961120797}, loss='logcosh', channel_map={'bsa_dubois': 0})

TMAPS['lv_mass'] = TensorMap(
    'lv_mass', Interpretation.CONTINUOUS, activation='linear', loss='logcosh', validator=make_range_validator(0, 500),
    channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8},
)
TMAPS['lv_mass_no0'] = TensorMap(
    'lv_mass', Interpretation.CONTINUOUS, activation='linear', loss=ignore_zeros_logcosh,
    channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8},
)

TMAPS['lv_mass_sentinel'] = TensorMap(
    'lv_mass', Interpretation.CONTINUOUS, activation='linear', sentinel=0,
    channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8},
)
TMAPS['LVM_sentinel'] = TensorMap(
    'LVM',  Interpretation.CONTINUOUS, normalization={'mean': 89.70372484725051, 'std': 24.803669503436304}, sentinel=0,
    validator=make_range_validator(-1, 300), channel_map={'LVM': 0},
)
TMAPS['lv_mass_prediction'] = TensorMap(
    'lv_mass_sentinel_prediction', Interpretation.CONTINUOUS, activation='linear', loss='logcosh', loss_weight=10.0,
    validator=make_range_validator(0, 300), channel_map={'lv_mass_sentinel_prediction': 0},
    normalization={'mean': 89.7, 'std': 24.8},
)
TMAPS['lv_mass_dubois_index_prediction'] = TensorMap(
    'lv_mass_dubois_index_sentinel_prediction', Interpretation.CONTINUOUS, activation='linear', loss='logcosh',
    validator=make_range_validator(0, 300), loss_weight=10.0,
    channel_map={'lv_mass_dubois_index_sentinel_prediction': 0}, normalization={'mean': 89.7, 'std': 24.8},
)
TMAPS['lv_mass_mosteller_index_prediction'] = TensorMap(
    'lv_mass_mosteller_index_sentinel_prediction', Interpretation.CONTINUOUS, activation='linear', loss='logcosh',
    validator=make_range_validator(0, 300), loss_weight=10.0,
    channel_map={'lv_mass_mosteller_index_sentinel_prediction': 0},
    normalization={'mean': 89.7, 'std': 24.8},
)

TMAPS['LVM_prediction'] = TensorMap(
    'LVM_sentinel_prediction',  Interpretation.CONTINUOUS, normalization={'mean': 89.70372484725051, 'std': 24.803669503436304},
    validator=make_range_validator(0, 300), channel_map={'LVM_sentinel_prediction': 0},
)

TMAPS['lvm_dubois_index_prediction'] = TensorMap(
    'lvm_dubois_index_sentinel_prediction', Interpretation.CONTINUOUS, activation='linear', loss='logcosh',
    validator=make_range_validator(0, 300), channel_map={'lvm_dubois_index_sentinel_prediction': 0},
    normalization={'mean': 42.0, 'std': 8.0},
)
TMAPS['lvm_mosteller_index_prediction'] = TensorMap(
    'lvm_mosteller_index_sentinel_prediction', Interpretation.CONTINUOUS, activation='linear', loss='logcosh',
    validator=make_range_validator(0, 300), channel_map={'lvm_mosteller_index_sentinel_prediction': 0},
    normalization={'mean': 42.0, 'std': 8.0},
)

TMAPS['LVM_prediction_sentinel'] = TensorMap(
    'LVM_sentinel_prediction',  Interpretation.CONTINUOUS, sentinel=0, channel_map={'LVM_sentinel_prediction': 0},
    normalization={'mean': 89.70372484725051, 'std': 24.803669503436304},
)
TMAPS['lvm_dubois_index_prediction_sentinel'] = TensorMap(
    'lvm_dubois_index_sentinel_prediction', Interpretation.CONTINUOUS, activation='linear', loss='logcosh',
    sentinel=0, channel_map={'lvm_dubois_index_sentinel_prediction': 0},
    normalization={'mean': 89.7, 'std': 24.8},
)
TMAPS['lvm_mosteller_index_prediction_sentinel'] = TensorMap(
    'lvm_mosteller_index_sentinel_prediction', Interpretation.CONTINUOUS, activation='linear', loss='logcosh',
    sentinel=0, channel_map={'lvm_mosteller_index_sentinel_prediction': 0},
    normalization={'mean': 89.7, 'std': 24.8},
)



TMAPS['end_systole_volume'] = TensorMap(
    'end_systole_volume', Interpretation.CONTINUOUS, activation='linear', validator=make_range_validator(0, 300),
    loss='logcosh', channel_map={'end_systole_volume': 0},
    normalization={'mean': 47.0, 'std': 10.0},
)
TMAPS['end_diastole_volume'] = TensorMap(
    'end_diastole_volume', Interpretation.CONTINUOUS, activation='linear', validator=make_range_validator(0, 400),
    loss='logcosh', channel_map={'end_diastole_volume': 0},
    normalization={'mean': 142.0, 'std': 21.0},
)
TMAPS['ejection_fraction'] = TensorMap(
    'ejection_fraction', Interpretation.CONTINUOUS, activation='linear', validator=make_range_validator(0.2, 0.9),
    normalization={'mean': 0.50, 'std': 0.046},
    loss='logcosh', loss_weight=1.0, channel_map={'ejection_fraction': 0},
)


# Apply correction from Sanghvi et al.Journal of Cardiovascular Magnetic Resonance 2016
TMAPS['corrected_extracted_lvedv'] = TensorMap(
    'corrected_extracted_lvedv', Interpretation.CONTINUOUS, activation='linear', validator=make_range_validator(0, 400),
    loss='logcosh', channel_map={'corrected_extracted_lvedv': 0},
    normalization={'mean': 142.0, 'std': 21.0},
)
TMAPS['corrected_extracted_lvef'] = TensorMap(
    'corrected_extracted_lvef', Interpretation.CONTINUOUS, activation='linear', validator=make_range_validator(0.2, 0.9),
    normalization={'mean': 0.50, 'std': 0.046},
    loss='logcosh', channel_map={'corrected_extracted_lvef': 0},
)
TMAPS['corrected_extracted_lvesv'] = TensorMap(
    'corrected_extracted_lvesv', Interpretation.CONTINUOUS, activation='linear', validator=make_range_validator(0, 300),
    loss='logcosh', channel_map={'corrected_extracted_lvesv': 0},
    normalization={'mean': 47.0, 'std': 10.0},
)

TMAPS['corrected_extracted_lvesv_sentinel'] = TensorMap(
    'corrected_extracted_lvesv', Interpretation.CONTINUOUS, activation='linear', sentinel=0.0,
    channel_map={'corrected_extracted_lvesv': 0}, normalization={'mean': 47.0, 'std': 10.0},
)
TMAPS['corrected_extracted_lvedv_sentinel'] = TensorMap(
    'corrected_extracted_lvedv', Interpretation.CONTINUOUS, activation='linear', sentinel=0.0,
    channel_map={'corrected_extracted_lvedv': 0}, normalization={'mean': 142.0, 'std': 21.0},
)
TMAPS['corrected_extracted_lvef_sentinel'] = TensorMap(
    'corrected_extracted_lvef', Interpretation.CONTINUOUS, activation='linear', sentinel=0.0,
    normalization={'mean': 0.50, 'std': 0.046}, channel_map={'corrected_extracted_lvef': 0},
)
TMAPS['corrected_extracted_lvef_sentinel'] = TensorMap(
    'corrected_extracted_lvef', Interpretation.CONTINUOUS, activation='linear', sentinel=0.0,
    normalization={'mean': 0.50, 'std': 0.046}, channel_map={'corrected_extracted_lvef': 0},
)

TMAPS['LA_2Ch_vol_max'] = TensorMap(
    'LA_2Ch_vol_max',  Interpretation.CONTINUOUS, normalization={'mean': 63.45582391534391, 'std': 22.548034481265972},
    validator=make_range_validator(0, 400), loss='logcosh', channel_map={'LA_2Ch_vol_max': 0},
)
TMAPS['LA_2Ch_vol_min'] = TensorMap(
    'LA_2Ch_vol_min',  Interpretation.CONTINUOUS, normalization={'mean': 28.308681904761904, 'std': 15.842444310837582},
    validator=make_range_validator(0, 200), loss='logcosh', channel_map={'LA_2Ch_vol_min': 0},
)
TMAPS['LA_4Ch_vol_max'] = TensorMap(
    'LA_4Ch_vol_max',  Interpretation.CONTINUOUS, normalization={'mean': 74.53903305263158, 'std': 25.448756860639776},
    validator=make_range_validator(0, 400), loss='logcosh', channel_map={'LA_4Ch_vol_max': 0},
)
TMAPS['LA_4Ch_vol_min'] = TensorMap(
    'LA_4Ch_vol_min',  Interpretation.CONTINUOUS, normalization={'mean': 31.014961894736846, 'std': 17.146722819760804},
    validator=make_range_validator(0, 200), loss='logcosh', channel_map={'LA_4Ch_vol_min': 0},
)
TMAPS['LA_Biplan_vol_max'] = TensorMap(
    'LA_Biplan_vol_max',  Interpretation.CONTINUOUS, normalization={'mean': 67.86355108225109, 'std': 21.793845470012105},
    validator=make_range_validator(0, 400), loss='logcosh', channel_map={'LA_Biplan_vol_max': 0},
)
TMAPS['LA_Biplan_vol_min'] = TensorMap(
    'LA_Biplan_vol_min',  Interpretation.CONTINUOUS, normalization={'mean': 28.79685670995671, 'std': 15.43219634139272},
    validator=make_range_validator(0, 300), loss='logcosh', channel_map={'LA_Biplan_vol_min': 0},
)
TMAPS['LVEDV'] = TensorMap(
    'LVEDV',  Interpretation.CONTINUOUS, normalization={'mean': 144.1479505192425, 'std': 34.39409859908663}, loss='logcosh',
    validator=make_range_validator(0, 500), channel_map={'LVEDV': 0},
)
TMAPS['LVESV'] = TensorMap(
    'LVESV',  Interpretation.CONTINUOUS, normalization={'mean': 59.58324862553452, 'std': 21.186976544044025}, loss='logcosh',
    validator=make_range_validator(0, 400), channel_map={'LVESV': 0},
)
TMAPS['LVM'] = TensorMap(
    'LVM',  Interpretation.CONTINUOUS, normalization={'mean': 89.70372484725051, 'std': 24.803669503436304}, loss='logcosh',
    validator=make_range_validator(0, 400), channel_map={'LVM': 0},
)
TMAPS['LVSV'] = TensorMap(
    'LVSV',  Interpretation.CONTINUOUS, normalization={'mean': 84.85198120147119, 'std': 19.2700091046526}, loss='logcosh',
    validator=make_range_validator(0, 400), channel_map={'LVSV': 0},
)
TMAPS['RA_4Ch_vol_max'] = TensorMap(
    'RA_4Ch_vol_max',  Interpretation.CONTINUOUS, normalization={'mean': 79.22289586811351, 'std': 26.504015552539048},
    validator=make_range_validator(0, 500), loss='logcosh', channel_map={'RA_4Ch_vol_max': 0},
)
TMAPS['RA_4Ch_vol_min'] = TensorMap(
    'RA_4Ch_vol_min',  Interpretation.CONTINUOUS, normalization={'mean': 46.25831176961603, 'std': 20.002160080524803},
    validator=make_range_validator(0, 400), loss='logcosh', channel_map={'RA_4Ch_vol_min': 0},
)
TMAPS['RVEDV'] = TensorMap(
    'RVEDV',  Interpretation.CONTINUOUS, normalization={'mean': 152.41239853151131, 'std': 37.15198900632509}, loss='logcosh',
    validator=make_range_validator(0, 500), channel_map={'RVEDV': 0},
)
TMAPS['RVEF'] = TensorMap(
    'RVEF',  Interpretation.CONTINUOUS, normalization={'mean': 56.404863078182565, 'std': 6.526231365539632}, loss='logcosh',
    validator=make_range_validator(10, 200), channel_map={'RVEF': 0},
)
TMAPS['RVESV'] = TensorMap(
    'RVESV',  Interpretation.CONTINUOUS, normalization={'mean': 67.61379869467673, 'std': 22.853189258914284}, loss='logcosh',
    validator=make_range_validator(0, 300), channel_map={'RVESV': 0},
)
TMAPS['RVSV'] = TensorMap(
    'RVSV',  Interpretation.CONTINUOUS, normalization={'mean': 85.0908258288989, 'std': 19.30893645374548}, loss='logcosh',
    validator=make_range_validator(0, 200), channel_map={'RVSV': 0},
)
TMAPS['LAQC'] = TensorMap(
    'LAQC',  Interpretation.CONTINUOUS, normalization={'mean': 1.2657977883096367, 'std': 0.5561369836438385}, loss='logcosh',
    validator=make_range_validator(0, 200), channel_map={'LAQC': 0},
)
TMAPS['LVQC'] = TensorMap(
    'LVQC',  Interpretation.CONTINUOUS, normalization={'mean': 1.1737756714060033, 'std': 0.4620420984104567}, loss='logcosh',
    validator=make_range_validator(0, 200), channel_map={'LVQC': 0},
)
TMAPS['RAQC'] = TensorMap(
    'RAQC',  Interpretation.CONTINUOUS, normalization={'mean': 1.1860189573459716, 'std': 0.4791815490882246}, loss='logcosh',
    validator=make_range_validator(0, 200), channel_map={'RAQC': 0},
)
TMAPS['RVQC'] = TensorMap(
    'RVQC',  Interpretation.CONTINUOUS, normalization={'mean': 1.179699842022117, 'std': 0.4648958893626213}, loss='logcosh',
    validator=make_range_validator(0, 200), channel_map={'RVQC': 0},
)

TMAPS['myocardial_mass'] = TensorMap(
    'myocardium_mass',  Interpretation.CONTINUOUS, validator=make_range_validator(0, 400), loss='logcosh', path_prefix='continuous',
    channel_map={'myocardium_mass': 0}, normalization={'mean': 89.70, 'std': 24.80},
)
TMAPS['myocardial_mass_noheritable'] = TensorMap(
    'inferred_myocardial_mass_noheritable',  Interpretation.CONTINUOUS, path_prefix='continuous',
    loss='logcosh', validator=make_range_validator(0, 400), normalization={'mean': 89.70, 'std': 24.80},
    channel_map={'inferred_myocardial_mass_noheritable': 0},
)
TMAPS['myocardial_mass_noheritable_sentinel'] = TensorMap(
    'inferred_myocardial_mass_noheritable',  Interpretation.CONTINUOUS, sentinel=0, loss='logcosh',
    normalization={'mean': 89.70, 'std': 24.80}, path_prefix='continuous',
    channel_map={'inferred_myocardial_mass_noheritable': 0},
)

TMAPS['myocardial_mass'] = TensorMap(
    'myocardium_mass',  Interpretation.CONTINUOUS, validator=make_range_validator(0, 400), loss='logcosh', path_prefix='continuous',
    channel_map={'myocardium_mass': 0}, normalization={'mean': 89.70, 'std': 24.80},
)


TMAPS['adjusted_myocardium_mass_sentinel'] = TensorMap(
    'adjusted_myocardium_mass', Interpretation.CONTINUOUS, validator=make_range_validator(0, 400), path_prefix='continuous',
    loss='logcosh', channel_map={'adjusted_myocardium_mass': 0}, normalization={'mean': 89.70, 'std': 24.80},
    sentinel=0.0,
)

TMAPS['adjusted_myocardium_mass_mse'] = TensorMap(
    'adjusted_myocardium_mass', Interpretation.CONTINUOUS, validator=make_range_validator(0, 400), path_prefix='continuous',
    loss='mse', channel_map={'adjusted_myocardium_mass': 0}, normalization={'mean': 89.70, 'std': 24.80},
)
TMAPS['adjusted_myocardium_mass_y_true_mse'] = TensorMap(
    'adjusted_myocardium_mass', Interpretation.CONTINUOUS, validator=make_range_validator(0, 400), path_prefix='continuous',
    loss=y_true_times_mse, channel_map={'adjusted_myocardium_mass': 0}, normalization={'mean': 89.70, 'std': 24.80},
)

TMAPS['adjusted_myocardium_mass_y_true_sqr_mse'] = TensorMap(
    'adjusted_myocardium_mass', Interpretation.CONTINUOUS, validator=make_range_validator(0, 400), path_prefix='continuous',
    loss=y_true_squared_times_mse, channel_map={'adjusted_myocardium_mass': 0}, normalization={'mean': 89.70, 'std': 24.80},
)

TMAPS['adjusted_myocardium_mass_y_true_cube_mse'] = TensorMap(
    'adjusted_myocardium_mass', Interpretation.CONTINUOUS, validator=make_range_validator(0, 400), path_prefix='continuous',
    loss=y_true_cubed_times_mse, channel_map={'adjusted_myocardium_mass': 0}, normalization={'mean': 89.70, 'std': 24.80},
)


TMAPS['adjusted_myocardium_mass_y_true_sqr_logcosh'] = TensorMap(
    'adjusted_myocardium_mass', Interpretation.CONTINUOUS, validator=make_range_validator(0, 400),
    loss=y_true_squared_times_logcosh, channel_map={'adjusted_myocardium_mass': 0}, path_prefix='continuous',
    normalization={'mean': 89.70, 'std': 24.80},
)


TMAPS['proton_fat'] = TensorMap(
    '22402_Proton-density-fat-fraction-PDFF_2_0', Interpretation.CONTINUOUS, channel_map={'22402_Proton-density-fat-fraction-PDFF_2_0': 0},
    activation='linear', loss='logcosh',  annotation_units=1, path_prefix='continuous',
    validator=make_range_validator(0, 100), normalization={'mean': 3.91012, 'std': 4.64437},
)
TMAPS['liver_fat'] = TensorMap(
    '22402_Liver-fat-percentage_2_0', Interpretation.CONTINUOUS, channel_map={'22402_Liver-fat-percentage_2_0': 0},
    activation='linear', loss='logcosh',  annotation_units=1, path_prefix='continuous',
    validator=make_range_validator(0, 100), normalization={'mean': 3.91012, 'std': 4.64437},
)
TMAPS['liver_fat_sentinel'] = TensorMap(
    '22402_Liver-fat-percentage_2_0', Interpretation.CONTINUOUS, channel_map={'22402_Liver-fat-percentage_2_0': 0},
    normalization={'mean': 3.91012, 'std': 4.64437}, activation='linear', sentinel=0.0, path_prefix='continuous',
)
TMAPS['liver_fat_echo_predicted'] = TensorMap(
    'liver_fat_sentinel_prediction', Interpretation.CONTINUOUS, channel_map={'liver_fat_sentinel_prediction': 0},
    validator=make_range_validator(0, 100), normalization={'mean': 3.91012, 'std': 4.64437}, path_prefix='continuous', activation='linear', loss='logcosh',
)
TMAPS['liver_fat_echo_predicted_sentinel'] = TensorMap(
    'liver_fat_sentinel_prediction', Interpretation.CONTINUOUS, channel_map={'liver_fat_sentinel_prediction': 0},
    normalization={'mean': 3.91012, 'std': 4.64437}, activation='linear', path_prefix='continuous', sentinel=0.0,
)

TMAPS['gre_mullti_echo_10_te_liver'] = TensorMap('gre_mullti_echo_10_te_liver', shape=(160, 160, 10), loss='logcosh', normalization={'zero_mean_std1': 1.0})
TMAPS['gre_mullti_echo_10_te_liver_12bit'] = TensorMap('gre_mullti_echo_10_te_liver_12bit', shape=(160, 160, 10), loss='logcosh', normalization={'zero_mean_std1': 1.0})
TMAPS['lms_ideal_optimised_low_flip_6dyn'] = TensorMap('lms_ideal_optimised_low_flip_6dyn', shape=(232, 256, 36), loss='logcosh', normalization={'zero_mean_std1': 1.0})
TMAPS['lms_ideal_optimised_low_flip_6dyn_12bit'] = TensorMap('lms_ideal_optimised_low_flip_6dyn_12bit', shape=(232, 256, 36), loss='logcosh', normalization={'zero_mean_std1': 1.0})
TMAPS['lms_ideal_optimised_low_flip_6dyn_4slice'] = TensorMap('lms_ideal_optimised_low_flip_6dyn_4slice', shape=(232, 256, 4), loss='logcosh', normalization={'zero_mean_std1': 1.0})

TMAPS['shmolli_192i'] = TensorMap('shmolli_192i', shape=(288, 384, 7), normalization={'zero_mean_std1': 1.0})
TMAPS['shmolli_192i_liver'] = TensorMap('shmolli_192i_liver', shape=(288, 384, 7), normalization={'zero_mean_std1': 1.0})
TMAPS['shmolli_192i_12bit'] = TensorMap('shmolli_192i_12bit', shape=(288, 384, 7), normalization={'zero_mean_std1': 1.0})
TMAPS['shmolli_192i_fitparams'] = TensorMap('shmolli_192i_fitparams', shape=(288, 384, 7), normalization={'zero_mean_std1': 1.0})
TMAPS['shmolli_192i_t1map'] = TensorMap('shmolli_192i_t1map', shape=(288, 384, 2), normalization={'zero_mean_std1': 1.0})

TMAPS['sax_pixel_width'] = TensorMap(
    'mri_pixel_width_cine_segmented_sax_inlinevf', Interpretation.CONTINUOUS, annotation_units=2, channel_map={'mri_pixel_width_cine_segmented_sax_inlinevf': 0},
    validator=make_range_validator(0, 4), normalization={'mean': 1.83, 'std': 0.1},
)
TMAPS['sax_pixel_height'] = TensorMap(
    'mri_pixel_height_segmented_sax_inlinevf', Interpretation.CONTINUOUS, annotation_units=2, channel_map={'mri_pixel_height_cine_segmented_sax_inlinevf': 0},
    validator=make_range_validator(0, 4), normalization={'mean': 1.83, 'std': 0.1},
)

TMAPS['ejection_fractionp'] = TensorMap(
    'ejection_fraction', Interpretation.CONTINUOUS, activation='linear',
    normalization={'mean': 0.50, 'std': 0.046},
    loss='logcosh', loss_weight=1.0, channel_map={'ejection_fraction': 0},
    parents=[TMAPS['end_systole_volume'], TMAPS['end_diastole_volume']],
)

TMAPS['cine_segmented_sax_b1'] = TensorMap('cine_segmented_sax_b1', shape=(256, 256, 50), loss='mse')
TMAPS['cine_segmented_sax_b2'] = TensorMap('cine_segmented_sax_b2', shape=(256, 256, 50), loss='mse')
TMAPS['cine_segmented_sax_b4'] = TensorMap('cine_segmented_sax_b4', shape=(256, 256, 50), loss='mse')
TMAPS['cine_segmented_sax_b6'] = TensorMap('cine_segmented_sax_b6', shape=(256, 256, 50), loss='mse')

TMAPS['cine_segmented_lax_2ch'] = TensorMap('cine_segmented_lax_2ch', shape=(256, 256, 50), normalization={'zero_mean_std1': True})
TMAPS['cine_segmented_lax_3ch'] = TensorMap('cine_segmented_lax_3ch', shape=(256, 256, 50), normalization={'zero_mean_std1': True})
TMAPS['cine_segmented_lax_4ch'] = TensorMap('cine_segmented_lax_4ch', shape=(256, 256, 50), normalization={'zero_mean_std1': True})

TMAPS['cine_segmented_lax_2ch_4d'] = TensorMap('cine_segmented_lax_2ch_4d', shape=(256, 256, 50, 1), normalization={'zero_mean_std1': True})
TMAPS['cine_segmented_lax_3ch_4d'] = TensorMap('cine_segmented_lax_3ch_4d', shape=(256, 256, 50, 1), normalization={'zero_mean_std1': True})
TMAPS['cine_segmented_lax_4ch_4d'] = TensorMap('cine_segmented_lax_4ch_4d', shape=(256, 256, 50, 1), normalization={'zero_mean_std1': True})

TMAPS['lax-view-detect'] = TensorMap(
    'lax-view-detect', Interpretation.CATEGORICAL,
    channel_map={
        'cine_segmented_lax_2ch': 0, 'cine_segmented_lax_3ch': 1,
        'cine_segmented_lax_4ch': 2,
    },
)

TMAPS['sax-view-detect'] = TensorMap(
    'sax-view-detect', Interpretation.CATEGORICAL,
    channel_map={
        'cine_segmented_sax_b1': 0, 'cine_segmented_sax_b2': 1,
        'cine_segmented_sax_b3': 2, 'cine_segmented_sax_b4': 3,
        'cine_segmented_sax_b5': 4, 'cine_segmented_sax_b6': 5,
        'cine_segmented_sax_b7': 6, 'cine_segmented_sax_b8': 7,
        'cine_segmented_sax_b9': 8, 'cine_segmented_sax_b10': 9,
        'cine_segmented_sax_b11': 10,
    },
)

TMAPS['slax-view-detect'] = TensorMap(
    'slax-view-detect', Interpretation.CATEGORICAL,
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


TMAPS['genetic_pca_1'] = TensorMap(
    '22009_Genetic-principal-components_0_1', Interpretation.CONTINUOUS, path_prefix='continuous', normalization={'mean': -0.014422761536727896, 'std': 10.57799283718005},
    loss='logcosh', channel_map={'22009_Genetic-principal-components_0_1': 0},
)
TMAPS['genetic_pca_2'] = TensorMap(
    '22009_Genetic-principal-components_0_2', Interpretation.CONTINUOUS, path_prefix='continuous', #normalization={'mean': -0.014422761536727896, 'std': 10.57799283718005},
    loss='logcosh', activation='linear', channel_map={'22009_Genetic-principal-components_0_2': 0},
)
TMAPS['genetic_pca_3'] = TensorMap(
    '22009_Genetic-principal-components_0_3', Interpretation.CONTINUOUS, path_prefix='continuous', #normalization={'mean': -0.014422761536727896, 'std': 10.57799283718005},
    loss='logcosh', activation='linear', channel_map={'22009_Genetic-principal-components_0_3': 0},
)
TMAPS['genetic_pca_4'] = TensorMap(
    '22009_Genetic-principal-components_0_4', Interpretation.CONTINUOUS, path_prefix='continuous', #normalization={'mean': -0.014422761536727896, 'std': 10.57799283718005},
    loss='logcosh', activation='linear', channel_map={'22009_Genetic-principal-components_0_4': 0},
)
TMAPS['genetic_pca_5'] = TensorMap(
    '22009_Genetic-principal-components_0_5', Interpretation.CONTINUOUS, path_prefix='continuous', #normalization={'mean': -0.014422761536727896, 'std': 10.57799283718005},
    loss='logcosh', activation='linear', channel_map={'22009_Genetic-principal-components_0_5': 0},
)
TMAPS['genetic_pca_all5'] = TensorMap(
    'genetic_pca_all5', Interpretation.CONTINUOUS, path_prefix='continuous', normalization={'mean': -0.014422761536727896, 'std': 10.57799283718005},
    loss='logcosh', annotation_units=5, shape=(5,), activation='linear',
    channel_map={
        '22009_Genetic-principal-components_0_0': 0, '22009_Genetic-principal-components_0_1': 1,
        '22009_Genetic-principal-components_0_2': 2, '22009_Genetic-principal-components_0_3': 3,
        '22009_Genetic-principal-components_0_4': 4,
    },
)

TMAPS['genetic_caucasian'] = TensorMap(
    'Genetic-ethnic-grouping_Caucasian_0_0', Interpretation.CATEGORICAL,
    storage_type=StorageType.CATEGORICAL_FLAG, path_prefix='categorical',
    channel_map={'no_caucasian': 0, 'Genetic-ethnic-grouping_Caucasian_0_0': 1},
)
TMAPS['genetic_caucasian_weighted'] = TensorMap(
    'Genetic-ethnic-grouping_Caucasian_0_0', Interpretation.CATEGORICAL, storage_type=StorageType.CATEGORICAL_FLAG,
    path_prefix='categorical', channel_map={'no_caucasian': 0, 'Genetic-ethnic-grouping_Caucasian_0_0': 1},
    loss=weighted_crossentropy([10.0, 1.0], 'caucasian_loss'),
)

TMAPS['mothers_age'] = TensorMap(
    'mothers_age_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'mother_age': 0, 'mother_alive': 2, 'mother_dead': 3, 'not-missing': 1},
    normalization={'mean': 75.555, 'std': 11.977}, annotation_units = 4,
)

TMAPS['fathers_age'] = TensorMap(
    'fathers_age_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'father_age': 0, 'father_alive': 2, 'father_dead': 3, 'not-missing': 1},
    normalization={'mean':70.928, 'std': 12.746}, annotation_units = 4,
)

TMAPS['genetic_sex'] = TensorMap(
    'Genetic-sex_Male_0_0', Interpretation.CATEGORICAL, storage_type=StorageType.CATEGORICAL_FLAG, path_prefix='categorical', annotation_units=2,
    channel_map={'Genetic-sex_Female_0_0': 0, 'Genetic-sex_Male_0_0': 1}, loss='categorical_crossentropy',
)
TMAPS['sex'] = TensorMap(
    'Sex_Male_0_0', Interpretation.CATEGORICAL, storage_type=StorageType.CATEGORICAL_FLAG, path_prefix='categorical', annotation_units=2,
    channel_map={'Sex_Female_0_0': 0, 'Sex_Male_0_0': 1}, loss='categorical_crossentropy',
)
TMAPS['bmi'] = TensorMap(
    '23104_Body-mass-index-BMI_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'23104_Body-mass-index-BMI_0_0': 0}, annotation_units=1,
    validator=make_range_validator(0, 300), normalization={'mean': 27.432061533712652, 'std': 4.785244772462738}, loss='logcosh',
)
TMAPS['bmi_ukb'] = TensorMap(
    'bmi', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'23104_Body-mass-index-BMI_0_0': 0}, annotation_units=1,
    validator=make_range_validator(0, 300), normalization={'mean': 27.432061533712652, 'std': 4.785244772462738}, loss='logcosh',
)
TMAPS['birth_year'] = TensorMap(
    '22200_Year-of-birth_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'22200_Year-of-birth_0_0': 0}, annotation_units=1, loss='logcosh',
    validator=make_range_validator(1901, 2025), normalization={'mean': 1952.0639129359386, 'std': 7.656326148519739},
)
TMAPS['birth_year_34'] = TensorMap(
    '34_Year-of-birth_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'34_Year-of-birth_0_0': 0}, annotation_units=1, loss='logcosh',
    validator=make_range_validator(1901, 2025), normalization = {'mean': 1952.0639129359386, 'std': 7.656326148519739},
)
TMAPS['age_0'] = TensorMap(
    '21003_Age-when-attended-assessment-centre_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh', validator=make_range_validator(1, 120),
    normalization={'mean': 56.52847159208494, 'std': 8.095287610193827}, channel_map={'21003_Age-when-attended-assessment-centre_0_0': 0},
)
TMAPS['age_1'] = TensorMap(
    '21003_Age-when-attended-assessment-centre_1_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh', validator=make_range_validator(1, 120),
    normalization={'mean': 61.4476555588322, 'std': 7.3992113757847005}, channel_map={'21003_Age-when-attended-assessment-centre_1_0': 0},
)
TMAPS['age_2'] = TensorMap(
    '21003_Age-when-attended-assessment-centre_2_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh', validator=make_range_validator(1, 120),
    normalization={'mean': 63.35798891483556, 'std': 7.554638350423902}, channel_map={'21003_Age-when-attended-assessment-centre_2_0': 0},
)

TMAPS['brain_volume'] = TensorMap(
    '25010_Volume-of-brain-greywhite-matter_2_0', Interpretation.CONTINUOUS, path_prefix='continuous', normalization={'mean': 1165940.0, 'std': 111511.0},
    channel_map={'25010_Volume-of-brain-greywhite-matter_2_0': 0}, loss='logcosh', loss_weight=0.1,
)

TMAPS['sodium'] = TensorMap(
    '30530_Sodium-in-urine_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'30530_Sodium-in-urine_0_0': 0},
    normalization={'mean': 77.45323967267045, 'std': 44.441236848463774}, annotation_units=1, loss='logcosh',
)
TMAPS['potassium'] = TensorMap(
    '30520_Potassium-in-urine_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'30520_Potassium-in-urine_0_0': 0},
    normalization={'mean': 63.06182700345117, 'std': 33.84208704773539}, annotation_units=1, loss='logcosh',
)
TMAPS['cholesterol_hdl'] = TensorMap(
    '30760_HDL-cholesterol_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'30760_HDL-cholesterol_0_0': 0},
    normalization={'mean': 1.4480129055069355, 'std': 0.3823115953478376}, annotation_units=1, loss='logcosh',
)
TMAPS['cholesterol'] = TensorMap(
    '30690_Cholesterol_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'30690_Cholesterol_0_0': 0},
    normalization={'mean': 5.692381214399044, 'std': 1.1449409331668705}, annotation_units=1, loss='logcosh',
)

TMAPS['cigarettes'] = TensorMap('2887_Number-of-cigarettes-previously-smoked-daily_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'2887_Number-of-cigarettes-previously-smoked-daily_0_0': 0}, normalization = {'mean': 18.92662147068755, 'std':10.590930376362259 }, annotation_units=1)
TMAPS['alcohol'] = TensorMap('5364_Average-weekly-intake-of-other-alcoholic-drinks_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'5364_Average-weekly-intake-of-other-alcoholic-drinks_0_0': 0}, normalization = {'mean': 0.03852570253005904, 'std':0.512608370266108 }, annotation_units=1)


def alcohol_channel_map(instance=0, array_idx=0):
    return {
        f'Alcohol-intake-frequency_Never_{instance}_{array_idx}': 0,
        f'Alcohol-intake-frequency_Special-occasions-only_{instance}_{array_idx}': 1,
        f'Alcohol-intake-frequency_One-to-three-times-a-month_{instance}_{array_idx}': 2,
        f'Alcohol-intake-frequency_Once-or-twice-a-week_{instance}_{array_idx}': 3,
        f'Alcohol-intake-frequency_Three-or-four-times-a-week_{instance}_{array_idx}': 4,
        f'Alcohol-intake-frequency_Daily-or-almost-daily_{instance}_{array_idx}': 5,
    }


TMAPS['alcohol_0'] = TensorMap('alcohol_0', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_channel_map(instance=0))
TMAPS['alcohol_1'] = TensorMap('alcohol_1', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_channel_map(instance=1))
TMAPS['alcohol_2'] = TensorMap('alcohol_2', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_channel_map(instance=2))


def alcohol_status_map(instance=0, array_idx=0):
    return {
        f'Alcohol-drinker-status_Never_{instance}_{array_idx}': 0,
        f'Alcohol-drinker-status_Previous_{instance}_{array_idx}': 1,
        f'Alcohol-drinker-status_Current_{instance}_{array_idx}': 2,
    }


TMAPS['alcohol_status_0'] = TensorMap('alcohol_status_0', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_status_map(instance=0))
TMAPS['alcohol_status_1'] = TensorMap('alcohol_status_1', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_status_map(instance=1))
TMAPS['alcohol_status_2'] = TensorMap('alcohol_status_2', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_status_map(instance=2))


def alcohol_meals_map(instance=0, array_idx=0):
    return {
        f'Alcohol-usually-taken-with-meals_No_{instance}_{array_idx}': 0,
        f'Alcohol-usually-taken-with-meals_It-varies_{instance}_{array_idx}': 1,
        f'Alcohol-usually-taken-with-meals_Yes_{instance}_{array_idx}': 2,
    }


TMAPS['alcohol_meals_0'] = TensorMap('alcohol_meals_0', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_meals_map(instance=0))
TMAPS['alcohol_meals_1'] = TensorMap('alcohol_meals_1', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_meals_map(instance=1))
TMAPS['alcohol_meals_2'] = TensorMap('alcohol_meals_2', Interpretation.CATEGORICAL, path_prefix='categorical', channel_map=alcohol_meals_map(instance=2))

TMAPS['coffee'] = TensorMap(
    '1498_Coffee-intake_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'1498_Coffee-intake_0_0': 0},
    normalization={'mean': 2.015086529948216, 'std': 2.0914960998390497}, annotation_units=1,
)
TMAPS['water'] = TensorMap(
    '1528_Water-intake_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'1528_Water-intake_0_0': 0},
    normalization={'mean': 2.7322977785723324, 'std': 2.261996814128837}, annotation_units=1,
)
TMAPS['meat'] = TensorMap(
    '3680_Age-when-last-ate-meat_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'3680_Age-when-last-ate-meat_0_0': 0},
    normalization={'mean': 29.74062983480561, 'std': 14.417292213873964}, annotation_units=1,
)
TMAPS['walks'] = TensorMap(
    '864_Number-of-daysweek-walked-10-minutes_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'864_Number-of-daysweek-walked-10-minutes_0_0': 0},
    normalization={'mean': 5.369732285440756, 'std': 1.9564911925721618}, annotation_units=1,
)
TMAPS['walk_duration'] = TensorMap(
    '874_Duration-of-walks_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'874_Duration-of-walks_0_0': 0},
    normalization={'mean': 61.64092215093373, 'std': 78.79522990818906}, annotation_units=1,
)
TMAPS['physical_activities'] = TensorMap(
    '884_Number-of-daysweek-of-moderate-physical-activity-10-minutes_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'884_Number-of-daysweek-of-moderate-physical-activity-10-minutes_0_0': 0 },
    normalization={'mean': 3.6258833281089258, 'std': 2.3343738999823676}, annotation_units=1,
)
TMAPS['physical_activity'] = TensorMap(
    '894_Duration-of-moderate-activity_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'894_Duration-of-moderate-activity_0_0': 0 },
    normalization={'mean': 66.2862593866103, 'std': 77.28681218835422}, annotation_units=1,
)
TMAPS['physical_activity_vigorous'] = TensorMap(
    '904_Number-of-daysweek-of-vigorous-physical-activity-10-minutes_0_0', Interpretation.CONTINUOUS,
    channel_map={'904_Number-of-daysweek-of-vigorous-physical-activity-10-minutes_0_0': 0}, path_prefix='continuous',
    normalization={'mean': 1.838718301735063, 'std': 1.9593505421480895}, annotation_units=1,
)
TMAPS['physical_activity_vigorous_duration'] = TensorMap(
    '914_Duration-of-vigorous-activity_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'914_Duration-of-vigorous-activity_0_0': 0},
    normalization={'mean': 44.854488382965144, 'std': 48.159967071781466}, annotation_units=1,
)
TMAPS['tv'] = TensorMap(
    '1070_Time-spent-watching-television-TV_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'1070_Time-spent-watching-television-TV_0_0': 0},
    normalization={'mean': 2.7753595642790914, 'std': 1.7135478462887321}, annotation_units=1,
)
TMAPS['computer'] = TensorMap(
    '1080_Time-spent-using-computer_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'1080_Time-spent-using-computer_0_0': 0},
    normalization={'mean': 0.9781465855433753, 'std': 1.4444414103121512}, annotation_units=1,
)
TMAPS['car'] = TensorMap(
    '1090_Time-spent-driving_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', channel_map={'1090_Time-spent-driving_0_0': 0},
    normalization={'mean': 0.8219851505445748, 'std': 1.304094814200189}, annotation_units=1,
)
TMAPS['summer'] = TensorMap(
    '1050_Time-spend-outdoors-in-summer_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'1050_Time-spend-outdoors-in-summer_0_0': 0},
    normalization={'mean': 3.774492304870845, 'std': 2.430483731404539}, annotation_units=1,
)
TMAPS['winter'] = TensorMap(
    '1060_Time-spent-outdoors-in-winter_0_0', Interpretation.CONTINUOUS, path_prefix='continuous',
    channel_map={'1060_Time-spent-outdoors-in-winter_0_0': 0},
    normalization={'mean': 1.8629686916635555, 'std': 1.88916218603397}, annotation_units=1,
)

TMAPS['systolic_blood_pressure_0'] = TensorMap(
    '4080_Systolic-blood-pressure-automated-reading_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh',
    channel_map={'4080_Systolic-blood-pressure-automated-reading_0_0': 0}, validator=make_range_validator(40, 400),
    normalization={'mean': 137.79964191990328, 'std': 19.292863700283757},
)
TMAPS['diastolic_blood_pressure_0'] = TensorMap(
    '4079_Diastolic-blood-pressure-automated-reading_0_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh',
    channel_map={'4079_Diastolic-blood-pressure-automated-reading_0_0': 0}, validator=make_range_validator(20, 300),
    normalization={'mean': 82.20657551284782, 'std': 10.496040770224475},
)

TMAPS['systolic_blood_pressure_1'] = TensorMap(
    '4080_Systolic-blood-pressure-automated-reading_1_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh',
    channel_map={'4080_Systolic-blood-pressure-automated-reading_1_0': 0}, validator=make_range_validator(40, 400),
    normalization={'mean': 137.79964191990328, 'std': 19.292863700283757},
)
TMAPS['diastolic_blood_pressure_1'] = TensorMap(
    '4079_Diastolic-blood-pressure-automated-reading_1_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh',
    channel_map={'4079_Diastolic-blood-pressure-automated-reading_1_0': 0}, validator=make_range_validator(20, 300),
    normalization={'mean': 82.20657551284782, 'std': 10.496040770224475},
)

TMAPS['systolic_blood_pressure_2'] = TensorMap(
    '4080_Systolic-blood-pressure-automated-reading_2_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh',
    channel_map={'4080_Systolic-blood-pressure-automated-reading_2_0': 0}, validator=make_range_validator(40, 400),
    normalization={'mean': 137.79964191990328, 'std': 19.292863700283757},
)
TMAPS['diastolic_blood_pressure_2'] = TensorMap(
    '4079_Diastolic-blood-pressure-automated-reading_2_0', Interpretation.CONTINUOUS, path_prefix='continuous', loss='logcosh',
    channel_map={'4079_Diastolic-blood-pressure-automated-reading_2_0': 0}, validator=make_range_validator(20, 300),
    normalization={'mean': 82.20657551284782, 'std': 10.496040770224475},
)


TMAPS['categorical-phenotypes-25'] = TensorMap(
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

TMAPS['categorical-phenotypes-36'] = TensorMap(
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

TMAPS['categorical-phenotypes-78'] = TensorMap(
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

TMAPS['categorical-phenotypes-134'] = TensorMap(
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


TMAPS['ecg-bike-max-hr'] = TensorMap(
    'max_hr', path_prefix='ecg_bike', loss='logcosh', metrics=['mape'],
    normalization={'mean': 110.03, 'std': 20.04}, shape=(1,),
    tensor_from_file=normalized_first_date,
)
TMAPS['ecg-bike-resting-hr'] = TensorMap(
    'resting_hr', Interpretation.CONTINUOUS, path_prefix='ecg_bike', loss='logcosh', shape=(1,),
    metrics=['mape'], normalization={'mean': 71.2, 'std': 12.57},
    tensor_from_file=normalized_first_date,
)
TMAPS['ecg-bike-age'] = TensorMap(
    'age', Interpretation.CONTINUOUS, path_prefix='ecg_bike', loss='logcosh', metrics=['mape'], shape=(1,),
    normalization={'mean': 60, 'std': 7.65},
    tensor_from_file=normalized_first_date,
)
TMAPS['ecg-bike-max-pred-hr'] = TensorMap(
    'max_pred_hr', Interpretation.CONTINUOUS, path_prefix='ecg_bike', loss='logcosh', metrics=['mape'], shape=(1,),
    normalization={'mean': 167.5, 'std': 5.81},
    tensor_from_file=normalized_first_date,
)
TMAPS['ecg-bike-trend-hr'] = TensorMap(
    'trend_heartrate', Interpretation.CONTINUOUS, shape=(106, 1), path_prefix='ecg_bike',
    tensor_from_file=normalized_first_date,
)
TMAPS['ecg-bike-trend-load'] = TensorMap(
    'trend_load', Interpretation.CONTINUOUS, shape=(106, 1), path_prefix='ecg_bike',
    tensor_from_file=normalized_first_date,
)
TMAPS['ecg-bike-trend-grade'] = TensorMap(
    'trend_grade', Interpretation.CONTINUOUS, shape=(106, 1), path_prefix='ecg_bike',
    tensor_from_file=normalized_first_date,
)

TMAPS['ecg-bike-raw-trend-hr'] = TensorMap(
    'trend_heartrate', Interpretation.CONTINUOUS, shape=(87,), path_prefix='ukb_ecg_bike',
    tensor_from_file=normalized_first_date,
)
TMAPS['ecg-bike-raw-trend-load'] = TensorMap(
    'trend_load', Interpretation.CONTINUOUS, shape=(87,), path_prefix='ukb_ecg_bike',
    tensor_from_file=normalized_first_date,
)
TMAPS['ecg-bike-raw-trend-grade'] = TensorMap(
    'trend_grade', Interpretation.CONTINUOUS, shape=(87,), path_prefix='ukb_ecg_bike',
    tensor_from_file=normalized_first_date,
)
TMAPS['ecg-bike-raw-trend-artifact'] = TensorMap(
    'trend_artifact', Interpretation.CONTINUOUS, shape=(87,), path_prefix='ukb_ecg_bike',
    tensor_from_file=normalized_first_date,
)
TMAPS['ecg-bike-raw-trend-mets'] = TensorMap(
    'trend_mets', Interpretation.CONTINUOUS, shape=(87,), path_prefix='ukb_ecg_bike',
    tensor_from_file=normalized_first_date,
)
TMAPS['ecg-bike-raw-trend-pacecount'] = TensorMap(
    'trend_pacecount', Interpretation.CONTINUOUS, shape=(87,), path_prefix='ukb_ecg_bike',
    tensor_from_file=normalized_first_date,
)
TMAPS['ecg-bike-raw-trend-phasename'] = TensorMap(
    'trend_phasename', Interpretation.CONTINUOUS, shape=(87,), path_prefix='ukb_ecg_bike',
    tensor_from_file=normalized_first_date,
)
TMAPS['ecg-bike-raw-trend-phasetime'] = TensorMap(
    'trend_phasetime', Interpretation.CONTINUOUS, shape=(87,), path_prefix='ukb_ecg_bike',
    tensor_from_file=normalized_first_date,
)
TMAPS['ecg-bike-raw-trend-time'] = TensorMap(
    'trend_time', Interpretation.CONTINUOUS, shape=(87,), path_prefix='ukb_ecg_bike',
    tensor_from_file=normalized_first_date,
)
TMAPS['ecg-bike-raw-trend-vecount'] = TensorMap(
    'trend_vecount', Interpretation.CONTINUOUS, shape=(87,), path_prefix='ukb_ecg_bike',
    tensor_from_file=normalized_first_date,
)
TMAPS['ecg-bike-raw-full'] = TensorMap(
    'full', Interpretation.CONTINUOUS, shape=(216500, 3), path_prefix='ukb_ecg_bike',
    tensor_from_file=normalized_first_date,
)
