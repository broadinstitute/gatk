import numpy as np

from ml4cvd.TensorMap import TensorMap

from ml4cvd.metrics import weighted_crossentropy, ignore_zeros_l2, ignore_zeros_logcosh
from ml4cvd.defines import MRI_SEGMENTED, MRI_ZOOM_MASK, ECG_BIKE_FULL_SIZE, ECG_BIKE_MEDIAN_SIZE, ECG_BIKE_STRIP_SIZE, ECG_CHAR_2_IDX, IMPUTATION_RANDOM, ECG_BIKE_RECOVERY_SIZE


def _get_lead_cm(length):
    lead_cm = {}
    lead_weights = []
    for i in range(length):
        wave_val = i - (length//2)
        lead_cm['w'+str(wave_val).replace('-', '_')] = i
        lead_weights.append((np.abs(wave_val+1)/(length/2)) + 1.0)
    return lead_cm, lead_weights


diploid_cm = {'homozygous_reference': 0, 'heterozygous': 1, 'homozygous_variant': 2}

TMAPS = dict()

TMAPS['rs3829740'] = TensorMap('rs3829740', group='categorical_index', channel_map=diploid_cm)
TMAPS['rs2234962'] = TensorMap('rs2234962', group='categorical_index', channel_map=diploid_cm)
TMAPS['rs2042995'] = TensorMap('rs2042995', group='categorical_index', channel_map=diploid_cm)

TMAPS['rs3829740_weighted'] = TensorMap('rs3829740', group='categorical_index', channel_map=diploid_cm, loss=weighted_crossentropy([1, 1, 1.5], 'rs3829740'))
TMAPS['rs2234962_weighted'] = TensorMap('rs2234962', group='categorical_index', channel_map=diploid_cm, loss=weighted_crossentropy([.8, 1, 1.5], 'rs2234962'))
TMAPS['rs2042995_weighted'] = TensorMap('rs2042995', group='categorical_index', channel_map=diploid_cm, loss=weighted_crossentropy([.6, 1.5, 2], 'rs2042995'))


TMAPS['akap9_lof'] = TensorMap('AKAP9', group='categorical_flag', channel_map={'no_akap9_lof': 0, 'akap9_lof': 1})
TMAPS['dsc2_lof'] = TensorMap('DSC2', group='categorical_flag', channel_map={'no_dsc2_lof': 0, 'dsc2_lof': 1})
TMAPS['ryr2_lof'] = TensorMap('RYR2', group='categorical_flag', channel_map={'no_ryr2_lof': 0, 'ryr2_lof': 1})
TMAPS['ttn_lof'] = TensorMap('TTN', group='categorical_flag', channel_map={'no_ttn_lof': 0, 'ttn_lof': 1})
TMAPS['ttntv'] = TensorMap('has_ttntv', group='categorical_flag', channel_map={'no_TTN_tv': 0, 'TTN_tv': 1})

TMAPS['ecg_rest'] = TensorMap('strip', shape=(5000, 12), group='ecg_rest',
        channel_map={'strip_I': 0, 'strip_II': 1, 'strip_III': 2, 'strip_V1': 3, 'strip_V2': 4, 'strip_V3': 5,
                     'strip_V4': 6, 'strip_V5': 7, 'strip_V6': 8, 'strip_aVF': 9, 'strip_aVL': 10, 'strip_aVR': 11})

TMAPS['ecg_rest_fft'] = TensorMap('ecg_rest_fft', shape=(5000, 12), group='ecg_rest',
        channel_map={'strip_I': 0, 'strip_II': 1, 'strip_III': 2, 'strip_V1': 3, 'strip_V2': 4, 'strip_V3': 5,
                     'strip_V4': 6, 'strip_V5': 7, 'strip_V6': 8, 'strip_aVF': 9, 'strip_aVL': 10, 'strip_aVR': 11})

TMAPS['ecg_rest_stack'] = TensorMap('strip', shape=(600, 12, 8), group='ecg_rest',
        channel_map={'strip_I': 0, 'strip_II': 1, 'strip_III': 2, 'strip_V1': 3, 'strip_V2': 4, 'strip_V3': 5,
                     'strip_V4': 6, 'strip_V5': 7, 'strip_V6': 8, 'strip_aVF': 9, 'strip_aVL': 10, 'strip_aVR': 11})
TMAPS['ecg_rest_median'] = TensorMap('median', group='ecg_rest', shape=(600, 12), loss='logcosh', activation='linear',
                  metrics=['mse', 'mae', 'logcosh'],
                  channel_map={'median_I': 0, 'median_II': 1, 'median_III': 2, 'median_V1': 3, 'median_V2': 4,
                               'median_V3': 5, 'median_V4': 6, 'median_V5': 7, 'median_V6': 8, 'median_aVF': 9,
                               'median_aVL': 10, 'median_aVR': 11})

TMAPS['ecg_rest_median_stack'] = TensorMap('median', group='ecg_rest', shape=(600, 12, 1), activation='linear',
                                           metrics=['mse', 'mae', 'logcosh'], loss='logcosh', loss_weight=1.0,
                  channel_map={'median_I': 0, 'median_II': 1, 'median_III': 2, 'median_V1': 3, 'median_V2': 4,
                               'median_V3': 5, 'median_V4': 6, 'median_V5': 7, 'median_V6': 8, 'median_aVF': 9,
                               'median_aVL': 10, 'median_aVR': 11})

TMAPS['ecg_median_1lead'] = TensorMap('median', group='ecg_rest', shape=(600, 1), loss='logcosh', loss_weight=10.0,
                                      activation='linear', metrics=['mse', 'mae', 'logcosh'], channel_map={'lead': 0})
TMAPS['ecg_rest_1lead'] = TensorMap('strip', shape=(600, 8), group='ecg_rest', channel_map={'lead': 0}, dependent_map=TMAPS['ecg_median_1lead'])
TMAPS['ecg_median_1lead_categorical'] = TensorMap('median', group='categorical', shape=(600, 32), activation='softmax', 
                                                  channel_map=_get_lead_cm(32)[0],
                                                  loss=weighted_crossentropy(_get_lead_cm(32)[1], 'ecg_median_categorical'))
TMAPS['ecg_rest_1lead_categorical'] = TensorMap('strip', shape=(600, 8), group='ecg_rest',
                                                channel_map={'window0': 0, 'window1': 1, 'window2': 2, 'window3': 3, 'window4': 4, 'window5': 5, 'window6': 6, 'window7': 7},
                                                dependent_map=TMAPS['ecg_median_1lead_categorical'])

TMAPS['ecg_rhythm'] = TensorMap('ecg_rhythm', group='categorical', loss=weighted_crossentropy([2.0, 3.0, 3.0, 6.0], 'ecg_rhythm'),
                  channel_map={'Normal_sinus_rhythm': 0, 'Sinus_bradycardia': 1, 'Marked_sinus_bradycardia': 2, 'Atrial_fibrillation': 3})
TMAPS['ecg_coarse'] = TensorMap('ecg_coarse', group='categorical', loss=weighted_crossentropy([1.0, 15.0, 5.0], 'ecg_coarse'),
                                channel_map={'Sinus_rhythm': 0, 'Atrial_fibrillation': 1, 'Other_rhythm': 2})
TMAPS['ecg_semi_coarse'] = TensorMap('ecg_semi_coarse', group='categorical', loss=weighted_crossentropy([1.0, 1.0, 2.0, 4.0, 16.0, 20.0], 'ecg_semi_coarse'),
                                     channel_map={'Normal_sinus_rhythm': 0, 'Sinus_bradycardia': 1, 'Marked_sinus_bradycardia': 2, 'Other_sinus_rhythm': 3, 'Atrial_fibrillation': 4, 'Other_rhythm': 5})
TMAPS['ecg_semi_coarse_with_poor'] = TensorMap('ecg_semi_coarse_with_poor', group='categorical', loss=weighted_crossentropy([1.0, 1.0, 2.0, 4.0, 16.0, 20.0], 'ecg_semi_coarse_with_poor'),
                                     channel_map={'Normal_sinus_rhythm': 0, 'Sinus_bradycardia': 1, 'Marked_sinus_bradycardia': 2, 'Other_sinus_rhythm': 3, 'Atrial_fibrillation': 4, 'Other_rhythm': 5})
TMAPS['ecg_normal'] = TensorMap('ecg_normal', group='categorical', loss=weighted_crossentropy([2.0, 3.0, 3.0, 3.0], 'ecg_normal'),
                  channel_map={'Normal_ECG': 0, 'Abnormal_ECG': 1, 'Borderline_ECG': 2, 'Otherwise_normal_ECG': 3})
TMAPS['ecg_infarct'] = TensorMap('ecg_infarct', group='categorical', channel_map={'no_infarct': 0, 'infarct': 1},
                             loss=weighted_crossentropy([1.0, 8.0], 'ecg_infarct'))
TMAPS['ecg_poor_data'] = TensorMap('ecg_poor_data', group='categorical', channel_map={'no_poor_data_quality': 0, 'poor_data_quality': 1},
                             loss=weighted_crossentropy([1.0, 8.0], 'ecg_poor_data'))
TMAPS['ecg_block'] = TensorMap('ecg_block', group='categorical', channel_map={'no_block': 0, 'block': 1},
                             loss=weighted_crossentropy([1.0, 8.0], 'ecg_block'))


TMAPS['acute_mi'] = TensorMap('acute_mi', group='ecg_categorical_interpretation', channel_map={'no_acute_mi': 0, 'ACUTE MI': 1},
                              loss=weighted_crossentropy([0.1, 10.0], 'acute_mi'))

TMAPS['anterior_blocks'] = TensorMap('anterior_blocks', group='ecg_categorical_interpretation',
                                     channel_map={'no_anterior_blocks': 0, 'Left anterior fascicular block': 1, 'Left posterior fascicular block': 2},
                                     loss=weighted_crossentropy([0.1, 10.0, 10.0], 'anterior_blocks'))

TMAPS['av_block'] = TensorMap('av_block', group='ecg_categorical_interpretation', channel_map={'no_av_block': 0, 'st degree AV block': 1},
                              loss=weighted_crossentropy([0.1, 10.0], 'av_block'))

TMAPS['fine_rhythms'] = TensorMap('fine_rhythms', group='ecg_categorical_interpretation',
                                  loss=weighted_crossentropy([0.5, 2.0, 0.1, 10.0, 10.0, 10.0, 15.0, 2.0, 10.0, 0.5, 0.2, 5.0], 'fine_rhythms'),
                                   channel_map={'no_fine_rhythms': 0, 'Normal sinus rhythm with sinus arrhythmia': 1, 'Normal sinus rhythm': 2,
                                                'Sinus rhythm with fusion complexes': 3, 'Sinus rhythm with marked sinus arrhythmia': 4,
                                                'Sinus rhythm with short PR': 5, 'Sinus rhythm with sinus arrhythmia': 6,
                                                'Sinus rhythm with 1st degree AV block': 7, 'Sinus tachycardia': 8,
                                                'Marked sinus bradycardia': 9, 'Sinus bradycardia': 10, 'Atrial fibrillation': 11})

TMAPS['incomplete_right_bundle_branch_block'] = TensorMap('incomplete_right_bundle_branch_block', group='ecg_categorical_interpretation',
                              channel_map={'no_incomplete_right_bundle_branch_block': 0, 'Incomplete right bundle branch block': 1},
                              loss=weighted_crossentropy([0.1, 10.0], 'incomplete_right_bundle_branch_block'))

TMAPS['infarcts'] = TensorMap('infarcts', group='ecg_categorical_interpretation',
                              channel_map={'no_infarcts': 0, 'Anterior infarct': 1, 'Anteroseptal infarct': 2, 'Inferior infarct': 3, 'Lateral infarct': 4, 'Septal infarct': 5},
                              loss=weighted_crossentropy([0.1, 4.0, 6.0, 7.0, 6.0, 4.0], 'infarcts'))

TMAPS['left_atrial_enlargement'] = TensorMap('left_atrial_enlargement', group='ecg_categorical_interpretation',
                              channel_map={'no_left_atrial_enlargement': 0, 'Left atrial enlargement': 1},
                              loss=weighted_crossentropy([0.1, 10.0], 'left_atrial_enlargement'))

TMAPS['left_ventricular_hypertrophy'] = TensorMap('left_ventricular_hypertrophy', group='ecg_categorical_interpretation',
                              channel_map={'no_left_ventricular_hypertrophy': 0, 'Left ventricular hypertrophy': 1},
                              loss=weighted_crossentropy([0.1, 10.0], 'left_ventricular_hypertrophy'))

TMAPS['lvh_fine'] = TensorMap('lvh_fine', group='ecg_categorical_interpretation', loss=weighted_crossentropy([0.5, 2.0, 8.0, 8.0, 8.0], 'lvh_fine'),
                              channel_map={'no_lvh_fine': 0, 'Minimal voltage criteria for LVH may be normal variant': 1,
                                           'Moderate voltage criteria for LVH may be normal variant': 2, 'Voltage criteria for LVH may be normal variant': 3,
                                           'Left ventricular hypertrophy': 4})

TMAPS['poor_data_quality'] = TensorMap('poor_data_quality', group='ecg_categorical_interpretation', channel_map={'no_poor_data_quality': 0, 'Poor data quality': 1},
                                       loss=weighted_crossentropy([0.1, 3.0], 'poor_data_quality'))

TMAPS['premature_atrial_complexes'] = TensorMap('premature_atrial_complexes', group='ecg_categorical_interpretation',
                                                channel_map={'no_premature_atrial_complexes': 0, 'premature atrial complexes': 1},
                                                loss=weighted_crossentropy([0.1, 10.0], 'premature_atrial_complexes'))

TMAPS['premature_supraventricular_complexes'] = TensorMap('premature_supraventricular_complexes', group='ecg_categorical_interpretation',
                                                channel_map={'no_premature_supraventricular_complexes': 0, 'premature supraventricular complexes': 1},
                                                loss=weighted_crossentropy([0.1, 10.0], 'premature_supraventricular_complexes'))

TMAPS['premature_ventricular_complexes'] = TensorMap('premature_ventricular_complexes', group='ecg_categorical_interpretation',
                                                channel_map={'no_premature_ventricular_complexes': 0, 'premature ventricular complexes': 1},
                                                loss=weighted_crossentropy([0.1, 10.0], 'premature_ventricular_complexes'))

TMAPS['prolonged_qt'] = TensorMap('prolonged_qt', group='ecg_categorical_interpretation', channel_map={'no_prolonged_qt': 0, 'Prolonged QT': 1},
                                  loss=weighted_crossentropy([0.1, 10.0], 'prolonged_qt'))


TMAPS['ecg_rhythmp'] = TensorMap('ecg_rhythm', group='categorical', loss=weighted_crossentropy([2.0, 3.0, 3.0, 6.0], 'ecg_rhythmp'), activation='softmax',
                  channel_map={'Normal_sinus_rhythm': 0, 'Sinus_bradycardia': 1, 'Marked_sinus_bradycardia': 2, 'Atrial_fibrillation': 3}, parents=['output_median_ecg_rest'])
TMAPS['ecg_normalp'] = TensorMap('ecg_normal', group='categorical', loss=weighted_crossentropy([2.0, 3.0, 3.0, 3.0], 'ecg_normalp'), activation='softmax',
                  channel_map={'Normal_ECG': 0, 'Abnormal_ECG': 1, 'Borderline_ECG': 2, 'Otherwise_normal_ECG': 3}, parents=['output_median_ecg_rest'])
TMAPS['ecg_infarctp'] = TensorMap('ecg_infarct', group='categorical', channel_map={'no_infarct': 0, 'infarct': 1}, activation='softmax',
                             loss=weighted_crossentropy([1.0, 6.0], 'ecg_infarctp'), parents=['output_median_ecg_rest'])


TMAPS['ecg_rest_next_char'] = TensorMap('ecg_rest_next_char', shape=(len(ECG_CHAR_2_IDX),), channel_map=ECG_CHAR_2_IDX, activation='softmax', loss='categorical_crossentropy', loss_weight=2.0)
TMAPS['ecg_rest_text'] = TensorMap('ecg_rest_text', shape=(100, len(ECG_CHAR_2_IDX)), group='ecg_text', channel_map={'context': 0, 'alphabet': 1}, dependent_map=TMAPS['ecg_rest_next_char'])


TMAPS['p-axis'] = TensorMap('PAxis', group='continuous', channel_map={'PAxis': 0}, loss='logcosh',
                        normalization={'mean': 48.7, 'std': 23.1})
TMAPS['p-duration'] = TensorMap('PDuration', group='continuous', channel_map={'PDuration': 0}, loss='logcosh',
                            normalization={'mean': 96.1, 'std': 18.85})
TMAPS['p-offset'] = TensorMap('POffset', group='continuous', channel_map={'POffset': 0}, loss='logcosh',
                          normalization={'mean': 369.1, 'std': 28.42})
TMAPS['p-onset'] = TensorMap('POnset', group='continuous', channel_map={'POnset': 0}, loss='logcosh',
                         normalization={'mean': 275.1, 'std': 26.420})
TMAPS['pp-interval'] = TensorMap('PPInterval', group='continuous', channel_map={'PPInterval': 0}, loss='logcosh',
                             normalization={'mean': 1036.1, 'std': 185.0})
TMAPS['pq-interval'] = TensorMap('PQInterval', group='continuous', channel_map={'PQInterval': 0}, loss='logcosh',
                             normalization={'mean': 165.9, 'std': 26.3})
TMAPS['q-offset'] = TensorMap('QOffset', group='continuous', channel_map={'QOffset': 0}, loss='logcosh',
                          normalization={'mean': 525.1, 'std': 13.52})
TMAPS['q-onset'] = TensorMap('QOnset', group='continuous', channel_map={'QOnset': 0}, loss='logcosh',
                         normalization={'mean': 435.1, 'std': 11.420})
TMAPS['qrs-complexes'] = TensorMap('QRSComplexes', group='continuous', channel_map={'QRSDuration': 0}, loss='logcosh',
                               normalization={'mean': 8.0, 'std': 20.0})
TMAPS['qrs-duration'] = TensorMap('QRSDuration', group='continuous', channel_map={'QRSDuration': 0}, loss='logcosh',
                              normalization={'mean': 89.53, 'std': 12.21})
TMAPS['qrs-num'] = TensorMap('QRSNum', group='continuous', channel_map={'QRSNum': 0}, loss='logcosh',
                         normalization={'mean': 9.61, 'std': 1.64})
TMAPS['qt-interval'] = TensorMap('QTInterval', group='continuous', channel_map={'QTInterval': 0}, loss='logcosh',
                             normalization={'mean': 426.1, 'std': 32.24})
TMAPS['qtc-interval'] = TensorMap('QTCInterval', group='continuous', channel_map={'QTCInterval': 0}, loss='logcosh',
                              normalization={'mean': 419.1, 'std': 20.7})
TMAPS['r-axis'] = TensorMap('RAxis', group='continuous', channel_map={'RAxis': 0}, loss='logcosh',
                        normalization={'mean': 25.7, 'std': 36.6})
TMAPS['rr-interval'] = TensorMap('RRInterval', group='continuous', channel_map={'RRInterval': 0}, loss='logcosh',
                             normalization={'mean': 1040.61, 'std': 175.5})
TMAPS['ventricular-rate'] = TensorMap('VentricularRate', group='continuous', channel_map={'VentricularRate': 0},
                                  loss='logcosh', normalization={'mean': 59.3, 'std': 10.6})
TMAPS['t-offset'] = TensorMap('TOffset', group='continuous', channel_map={'TOffset': 0}, loss='logcosh',
                          normalization={'mean': 860.7, 'std': 32.52})
TMAPS['t-axis'] = TensorMap('TAxis', group='continuous', channel_map={'TAxis': 0}, loss='logcosh',
                        normalization={'mean': 40.8, 'std': 32.6})

TMAPS['p-axisp'] = TensorMap('PAxis', group='continuous', channel_map={'PAxis': 0}, loss='logcosh',
                        normalization={'mean': 48.7, 'std': 23.1}, parents=['output_median_ecg_rest'])
TMAPS['p-durationp'] = TensorMap('PDuration', group='continuous', channel_map={'PDuration': 0}, loss='logcosh',
                            normalization={'mean': 96.1, 'std': 18.85}, parents=['output_median_ecg_rest'])
TMAPS['p-offsetp'] = TensorMap('POffset', group='continuous', channel_map={'POffset': 0}, loss='logcosh',
                          normalization={'mean': 369.1, 'std': 28.42}, parents=['output_median_ecg_rest'])
TMAPS['p-onsetp'] = TensorMap('POnset', group='continuous', channel_map={'POnset': 0}, loss='logcosh',
                         normalization={'mean': 275.1, 'std': 26.420}, parents=['output_median_ecg_rest'])
TMAPS['pp-intervalp'] = TensorMap('PPInterval', group='continuous', channel_map={'PPInterval': 0}, loss='logcosh',
                             normalization={'mean': 1036.1, 'std': 185.0}, parents=['output_median_ecg_rest'])
TMAPS['pq-intervalp'] = TensorMap('PQInterval', group='continuous', channel_map={'PQInterval': 0}, loss='logcosh',
                             normalization={'mean': 165.9, 'std': 26.3}, parents=['output_median_ecg_rest'])
TMAPS['q-offsetp'] = TensorMap('QOffset', group='continuous', channel_map={'QOffset': 0}, loss='logcosh',
                          normalization={'mean': 525.1, 'std': 13.52}, parents=['output_median_ecg_rest'])
TMAPS['q-onsetp'] = TensorMap('QOnset', group='continuous', channel_map={'QOnset': 0}, loss='logcosh',
                         normalization={'mean': 435.1, 'std': 11.420}, parents=['output_median_ecg_rest'])
TMAPS['qrs-complexesp'] = TensorMap('QRSComplexes', group='continuous', channel_map={'QRSDuration': 0}, loss='logcosh',
                               normalization={'mean': 8.0, 'std': 20.0}, parents=['output_median_ecg_rest'])
TMAPS['qrs-durationp'] = TensorMap('QRSDuration', group='continuous', channel_map={'QRSDuration': 0}, loss='logcosh',
                              normalization={'mean': 89.53, 'std': 12.21}, parents=['output_median_ecg_rest'])
TMAPS['qrs-nump'] = TensorMap('QRSNum', group='continuous', channel_map={'QRSNum': 0}, loss='logcosh',
                         normalization={'mean': 9.61, 'std': 1.64}, parents=['output_median_ecg_rest'])
TMAPS['qt-intervalp'] = TensorMap('QTInterval', group='continuous', channel_map={'QTInterval': 0}, loss='logcosh',
                             normalization={'mean': 426.1, 'std': 32.24}, parents=['output_median_ecg_rest'])
TMAPS['qtc-intervalp'] = TensorMap('QTCInterval', group='continuous', channel_map={'QTCInterval': 0}, loss='logcosh',
                              normalization={'mean': 419.1, 'std': 20.7}, parents=['output_QTInterval_continuous', 'output_RRInterval_continuous'])
TMAPS['r-axisp'] = TensorMap('RAxis', group='continuous', channel_map={'RAxis': 0}, loss='logcosh',
                        normalization={'mean': 25.7, 'std': 36.6}, parents=['output_median_ecg_rest'])
TMAPS['rr-intervalp'] = TensorMap('RRInterval', group='continuous', channel_map={'RRInterval': 0}, loss='logcosh',
                             normalization={'mean': 1040.61, 'std': 175.5}, parents=['output_median_ecg_rest'])
TMAPS['ventricular-ratep'] = TensorMap('VentricularRate', group='continuous', channel_map={'VentricularRate': 0},
                                  loss='logcosh', normalization={'mean': 59.3, 'std': 10.6}, parents=['output_median_ecg_rest'])
TMAPS['t-offsetp'] = TensorMap('TOffset', group='continuous', channel_map={'TOffset': 0}, loss='logcosh',
                          normalization={'mean': 860.7, 'std': 32.52}, parents=['output_median_ecg_rest'])
TMAPS['t-axisp'] = TensorMap('TAxis', group='continuous', channel_map={'TAxis': 0}, loss='logcosh',
                        normalization={'mean': 40.8, 'std': 32.6}, parents=['output_median_ecg_rest'])

TMAPS['charge'] = TensorMap('charge', group='continuous', channel_map={'charge': 0}, normalization={'mean': 12.0, 'std': 2.0})
TMAPS['af_prs'] = TensorMap('AF_PRS_LDscore', group='continuous', channel_map={'AF_PRS_LDscore': 0}, normalization={'mean': -1.0, 'std': 0.4})


TMAPS['p-axis-no0'] = TensorMap('PAxis', group='continuous', channel_map={'PAxis': 0}, loss=ignore_zeros_logcosh, metrics=['logcosh'],
                        normalization={'mean': 48.7, 'std': 23.1})
TMAPS['p-duration-no0'] = TensorMap('PDuration', group='continuous', channel_map={'PDuration': 0}, loss=ignore_zeros_logcosh, metrics=['logcosh'],
                            normalization={'mean': 96.1, 'std': 18.85})
TMAPS['p-offset-no0'] = TensorMap('POffset', group='continuous', channel_map={'POffset': 0}, loss=ignore_zeros_logcosh, metrics=['logcosh'],
                          normalization={'mean': 369.1, 'std': 28.42})
TMAPS['p-onset-no0'] = TensorMap('POnset', group='continuous', channel_map={'POnset': 0}, loss=ignore_zeros_logcosh, metrics=['logcosh'],
                         normalization={'mean': 275.1, 'std': 26.420})
TMAPS['pp-interval-no0'] = TensorMap('PPInterval', group='continuous', channel_map={'PPInterval': 0}, loss=ignore_zeros_logcosh, metrics=['logcosh'],
                             normalization={'mean': 1036.1, 'std': 185.0})
TMAPS['pq-interval-no0'] = TensorMap('PQInterval', group='continuous', channel_map={'PQInterval': 0}, loss=ignore_zeros_logcosh, metrics=['logcosh'],
                             normalization={'mean': 165.9, 'std': 26.3})


TMAPS['p-axis-sentinel'] = TensorMap('PAxis', group='continuous', channel_map={'PAxis': 0}, sentinel=0, metrics=['logcosh'],
                        normalization={'mean': 48.7, 'std': 23.1})
TMAPS['p-duration-sentinel'] = TensorMap('PDuration', group='continuous', channel_map={'PDuration': 0}, sentinel=0, metrics=['logcosh'],
                            normalization={'mean': 96.1, 'std': 18.85})
TMAPS['p-offset-sentinel'] = TensorMap('POffset', group='continuous', channel_map={'POffset': 0}, sentinel=0, metrics=['logcosh'],
                          normalization={'mean': 369.1, 'std': 28.42})
TMAPS['p-onset-sentinel'] = TensorMap('POnset', group='continuous', channel_map={'POnset': 0}, sentinel=0, metrics=['logcosh'],
                         normalization={'mean': 275.1, 'std': 26.420})
TMAPS['pp-interval-sentinel'] = TensorMap('PPInterval', group='continuous', channel_map={'PPInterval': 0},  sentinel=0, metrics=['logcosh'],
                             normalization={'mean': 1036.1, 'std': 185.0})
TMAPS['pq-interval-sentinel'] = TensorMap('PQInterval', group='continuous', channel_map={'PQInterval': 0}, sentinel=0, metrics=['logcosh'],
                                     normalization={'mean': 165.9, 'std': 26.3})

TMAPS['ecg-bike-max-hr-no0'] = TensorMap('bike_max_hr', group='continuous', channel_map={'bike_max_hr': 0},
                                         loss=ignore_zeros_logcosh, metrics=['logcosh'], normalization={'mean': 110.03, 'std': 20.04})
TMAPS['ecg-bike-resting-hr-no0'] = TensorMap('bike_resting_hr', group='continuous', channel_map={'bike_resting_hr': 0},
                                             loss=ignore_zeros_logcosh, metrics=['logcosh'], normalization={'mean': 71.2, 'std': 12.57})
TMAPS['ecg-bike-max-pred-hr-no0'] = TensorMap('bike_max_pred_hr', group='continuous', channel_map={'bike_max_pred_hr': 0},
                                             loss=ignore_zeros_logcosh, metrics=['logcosh'], normalization={'mean': 167.5, 'std': 5.78})

TMAPS['ecg_bike_0'] = TensorMap('full_0', shape=ECG_BIKE_FULL_SIZE, group='ecg_bike', channel_map={'I': 0, '2': 1, '3': 2})
TMAPS['ecg_bike_1'] = TensorMap('full_1', shape=ECG_BIKE_FULL_SIZE, group='ecg_bike', channel_map={'I': 0, '2': 1, '3': 2})
TMAPS['ecg_bike_recovery'] = TensorMap('full', shape=ECG_BIKE_RECOVERY_SIZE, group='ecg_bike_recovery', channel_map={'lead_I': 0, 'lead_2': 1, 'lead_3': 2})
TMAPS['ecg_bike_m0'] = TensorMap('median_0', shape=ECG_BIKE_MEDIAN_SIZE, group='ecg_bike',
                             channel_map={'I': 0, '2': 1, '3': 2})
TMAPS['ecg_bike_m1'] = TensorMap('median_1', shape=ECG_BIKE_MEDIAN_SIZE, group='ecg_bike',
                             channel_map={'I': 0, '2': 1, '3': 2})
TMAPS['ecg_bike_s0'] = TensorMap('strip_0', shape=ECG_BIKE_STRIP_SIZE, group='ecg_bike',
                             channel_map={'I': 0, '2': 1, '3': 2})
TMAPS['ecg_bike_s1'] = TensorMap('strip_1', shape=ECG_BIKE_STRIP_SIZE, group='ecg_bike',
                             channel_map={'I': 0, '2': 1, '3': 2})

TMAPS['aligned_distance'] = TensorMap('aligned_distance', group='continuous', channel_map={'zeros': 0})

TMAPS['lv_mass'] = TensorMap('lv_mass', group='continuous', activation='linear', loss='logcosh',
                             channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8})
TMAPS['lv_mass_no0'] = TensorMap('lv_mass', group='continuous', activation='linear', loss=ignore_zeros_logcosh,
                             channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8})

TMAPS['lv_mass_sentinel'] = TensorMap('lv_mass', group='continuous', activation='linear', sentinel=0,
                                      channel_map={'lv_mass': 0}, normalization={'mean': 89.7, 'std': 24.8})

TMAPS['end_systole_volume'] = TensorMap('end_systole_volume', group='continuous', activation='linear',
                                    loss='logcosh', channel_map={'end_systole_volume': 0},
                                    normalization={'mean': 47.0, 'std': 10.0})
TMAPS['end_diastole_volume'] = TensorMap('end_diastole_volume', group='continuous', activation='linear',
                                     loss='logcosh', channel_map={'end_diastole_volume': 0},
                                     normalization={'mean': 142.0, 'std': 21.0})
TMAPS['ejection_fraction'] = TensorMap('ejection_fraction', group='continuous', activation='linear',
                                   normalization={'mean': 0.50, 'std': 0.046},
                                   loss='logcosh', loss_weight=1.0, channel_map={'ejection_fraction': 0})


# Apply correction from Sanghvi et al.Journal of Cardiovascular Magnetic Resonance 2016
TMAPS['corrected_extracted_lvedv'] = TensorMap('corrected_extracted_lvedv', group='continuous', activation='linear',
                                     loss='logcosh', channel_map={'corrected_extracted_lvedv': 0},
                                     normalization={'mean': 142.0, 'std': 21.0})
TMAPS['corrected_extracted_lvef'] = TensorMap('corrected_extracted_lvef', group='continuous', activation='linear',
                                   normalization={'mean': 0.50, 'std': 0.046},
                                   loss='logcosh', channel_map={'corrected_extracted_lvef': 0})
TMAPS['corrected_extracted_lvesv'] = TensorMap('corrected_extracted_lvesv', group='continuous', activation='linear',
                                    loss='logcosh', channel_map={'corrected_extracted_lvesv': 0},
                                    normalization={'mean': 47.0, 'std': 10.0})

TMAPS['corrected_extracted_lvesv_sentinel'] = TensorMap('corrected_extracted_lvesv', group='continuous', activation='linear', sentinel=0.0,
                                                        channel_map={'corrected_extracted_lvesv': 0}, normalization={'mean': 47.0, 'std': 10.0})
TMAPS['corrected_extracted_lvedv_sentinel'] = TensorMap('corrected_extracted_lvedv', group='continuous', activation='linear', sentinel=0.0,
                                                        channel_map={'corrected_extracted_lvedv': 0}, normalization={'mean': 142.0, 'std': 21.0})
TMAPS['corrected_extracted_lvef_sentinel'] = TensorMap('corrected_extracted_lvef', group='continuous', activation='linear', sentinel=0.0,
                                                       normalization={'mean': 0.50, 'std': 0.046}, channel_map={'corrected_extracted_lvef': 0})


TMAPS['liver_fat'] = TensorMap('22402_Liver-fat-percentage_2_0', group='continuous', channel_map={'22402_Liver-fat-percentage_2_0': 0},
                               normalization={'mean': 3.91012, 'std': 4.64437}, activation='linear', loss='logcosh', loss_weight=1.0)
TMAPS['liver_fat_sentinel'] = TensorMap('22402_Liver-fat-percentage_2_0', group='continuous', channel_map={'22402_Liver-fat-percentage_2_0': 0},
                               normalization={'mean': 3.91012, 'std': 4.64437}, activation='linear', sentinel=0.0)
TMAPS['liver_fat_echo_predicted'] = TensorMap('liver_fat_sentinel_prediction', group='continuous', channel_map={'liver_fat_sentinel_prediction': 0},
                               normalization={'mean': 3.91012, 'std': 4.64437}, activation='linear', loss='logcosh')
TMAPS['liver_fat_echo_predicted_sentinel'] = TensorMap('liver_fat_sentinel_prediction', group='continuous', channel_map={'liver_fat_sentinel_prediction': 0},
                               normalization={'mean': 3.91012, 'std': 4.64437}, activation='linear', sentinel=0.0)

TMAPS['gre_mullti_echo_10_te_liver'] = TensorMap('gre_mullti_echo_10_te_liver', (160, 160, 10), group='root_array', loss='logcosh')
TMAPS['gre_mullti_echo_10_te_liver_12bit'] = TensorMap('gre_mullti_echo_10_te_liver_12bit', (160, 160, 10), group='root_array', loss='logcosh')
TMAPS['lms_ideal_optimised_low_flip_6dyn'] = TensorMap('lms_ideal_optimised_low_flip_6dyn', (232, 256, 36), group='root_array', loss='logcosh')
TMAPS['lms_ideal_optimised_low_flip_6dyn_12bit'] = TensorMap('lms_ideal_optimised_low_flip_6dyn_12bit', (232, 256, 36), group='root_array', loss='logcosh')
TMAPS['lms_ideal_optimised_low_flip_6dyn_4slice'] = TensorMap('lms_ideal_optimised_low_flip_6dyn_4slice', (232, 256, 4), loss='logcosh')

TMAPS['t1_p2_1mm_fov256_sag_ti_880'] = TensorMap('t1_p2_1mm_fov256_sag_ti_880', (256, 256, 416), group='root_array')
TMAPS['t1_brain_208z'] = TensorMap('t1_brain_208z', (256, 256, 208))
TMAPS['t1_brain_208z_half'] = TensorMap('t1_brain_208z_half', (256 // 2, 256 // 2, 208 // 2))
TMAPS['t1_brain_208z_quarter'] = TensorMap('t1_brain_208z_quarter', (256 // 4, 256 // 4, 208 // 4))
TMAPS['t1_brain_208z_3d'] = TensorMap('t1_brain_208z_3d', (256, 256, 208, 1))
TMAPS['t1_brain_208z_half_3d'] = TensorMap('t1_brain_208z_half_3d', (256 // 2, 256 // 2, 208 // 2, 1))
TMAPS['t1_brain_208z_quarter_3d'] = TensorMap('t1_brain_208z_quarter_3d', (256 // 4, 256 // 4, 208 // 4, 1))

TMAPS['shmolli_192i'] = TensorMap('shmolli_192i', (288, 384, 7), group='root_array')
TMAPS['shmolli_192i_12bit'] = TensorMap('shmolli_192i_12bit', (288, 384, 7), group='root_array')
TMAPS['shmolli_192i_fitparams'] = TensorMap('shmolli_192i_fitparams', (288, 384, 7), group='root_array')
TMAPS['shmolli_192i_t1map'] = TensorMap('shmolli_192i_t1map', (288, 384, 2), group='root_array')

TMAPS['mri_pixel_width'] = TensorMap('mri_pixel_width', group='continuous', annotation_units=1, channel_map={'mri_pixel_width': 0}, normalization={'mean': 1.83, 'std': 0.1})
TMAPS['mri_pixel_height'] = TensorMap('mri_pixel_height', group='continuous', annotation_units=1, channel_map={'mri_pixel_height': 0}, normalization={'mean': 1.83, 'std': 0.1})


TMAPS['end_systole_volumep'] = TensorMap('end_systole_volume', group='continuous', activation='linear',
                                     loss='logcosh', channel_map={'end_systole_volume': 0},
                                     normalization={'mean': 47.0, 'std': 10.0},
                                     parents=['output_' + MRI_SEGMENTED + '_categorical'])
TMAPS['end_diastole_volumep'] = TensorMap('end_diastole_volume', group='continuous', activation='linear',
                                      loss='logcosh', channel_map={'end_diastole_volume': 0},
                                      normalization={'mean': 142.0, 'std': 21.0},
                                      parents=['output_' + MRI_SEGMENTED + '_categorical'])
TMAPS['end_systole_volumepz'] = TensorMap('end_systole_volume', group='continuous', activation='linear',
                                      loss='logcosh', channel_map={'end_systole_volume': 0},
                                      normalization={'mean': 47.0, 'std': 10.0},
                                      parents=['output_' + MRI_ZOOM_MASK + '_categorical'])
TMAPS['end_diastole_volumepz'] = TensorMap('end_diastole_volume', group='continuous', activation='linear',
                                       loss='logcosh', channel_map={'end_diastole_volume': 0},
                                       normalization={'mean': 142.0, 'std': 21.0},
                                       parents=['output_' + MRI_ZOOM_MASK + '_categorical'])
TMAPS['ejection_fractionp'] = TensorMap('ejection_fraction', group='continuous', activation='linear',
                                    normalization={'mean': 0.50, 'std': 0.046},
                                    loss='logcosh', loss_weight=1.0, channel_map={'ejection_fraction': 0},
                                    parents=['output_end_systole_volume_continuous',
                                             'output_end_diastole_volume_continuous'])

TMAPS['cine_segmented_sax_b1'] = TensorMap('cine_segmented_sax_b1', (256, 256, 50), group='root_array', loss='mse')
TMAPS['cine_segmented_sax_b2'] = TensorMap('cine_segmented_sax_b2', (256, 256, 50), group='root_array', loss='mse')
TMAPS['cine_segmented_sax_b4'] = TensorMap('cine_segmented_sax_b4', (256, 256, 50), group='root_array', loss='mse')
TMAPS['cine_segmented_sax_b6'] = TensorMap('cine_segmented_sax_b6', (256, 256, 50), group='root_array', loss='mse')

TMAPS['cine_segmented_lax_2ch'] = TensorMap('cine_segmented_lax_2ch', (256, 256, 50), group='root_array', loss='logcosh')
TMAPS['cine_segmented_lax_3ch'] = TensorMap('cine_segmented_lax_3ch', (256, 256, 50), group='root_array', loss='logcosh')
TMAPS['cine_segmented_lax_4ch'] = TensorMap('cine_segmented_lax_4ch', (256, 256, 50), group='root_array', loss='logcosh')
TMAPS['lax-view-detect'] = TensorMap('lax-view-detect', group='categorical',
                                 channel_map={'cine_segmented_lax_2ch': 0, 'cine_segmented_lax_3ch': 1,
                                              'cine_segmented_lax_4ch': 2})

TMAPS['sax-view-detect'] = TensorMap('sax-view-detect', group='categorical',
                                 channel_map={'cine_segmented_sax_b1': 0, 'cine_segmented_sax_b2': 1,
                                              'cine_segmented_sax_b3': 2, 'cine_segmented_sax_b4': 3,
                                              'cine_segmented_sax_b5': 4, 'cine_segmented_sax_b6': 5,
                                              'cine_segmented_sax_b7': 6, 'cine_segmented_sax_b8': 7,
                                              'cine_segmented_sax_b9': 8, 'cine_segmented_sax_b10': 9,
                                              'cine_segmented_sax_b11': 10})

TMAPS['slax-view-detect'] = TensorMap('slax-view-detect', group='categorical',
                                  channel_map={'cine_segmented_lax_2ch': 11, 'cine_segmented_lax_3ch': 12,
                                               'cine_segmented_lax_4ch': 13, 'cine_segmented_sax_b1': 0,
                                               'cine_segmented_sax_b2': 1, 'cine_segmented_sax_b3': 2,
                                               'cine_segmented_sax_b4': 3, 'cine_segmented_sax_b5': 4,
                                               'cine_segmented_sax_b6': 5, 'cine_segmented_sax_b7': 6,
                                               'cine_segmented_sax_b8': 7, 'cine_segmented_sax_b9': 8,
                                               'cine_segmented_sax_b10': 9, 'cine_segmented_sax_b11': 10})

TMAPS['mothers_age'] = TensorMap('mothers_age_0', group='continuous',
                                 channel_map={'mother_age': 0, 'mother_alive': 2, 'mother_dead': 3, 'not-missing': 1},
                                 normalization={'mean': 75.555, 'std': 11.977}, annotation_units = 4)

TMAPS['fathers_age'] = TensorMap('fathers_age_0', group='continuous',
                                 channel_map={'father_age': 0, 'father_alive': 2, 'father_dead': 3, 'not-missing': 1},
                                 normalization={'mean':70.928, 'std': 12.746}, annotation_units = 4)

TMAPS['genetic_sex'] = TensorMap('genetic_sex', group='categorical', annotation_units=2, channel_map={'Genetic-sex_Female_0_0': 0, 'Genetic-sex_Male_0_0': 1}, loss='categorical_crossentropy')
TMAPS['sex'] = TensorMap('sex', group='categorical', annotation_units=2, channel_map={'Sex_Female_0_0': 0, 'Sex_Male_0_0': 1}, loss='categorical_crossentropy')
TMAPS['bmi'] = TensorMap('23104_Body-mass-index-BMI_0_0', group='continuous', channel_map={'23104_Body-mass-index-BMI_0_0':0}, normalization = {'mean': 27.432061533712652, 'std': 4.785244772462738}, annotation_units=1, loss='logcosh')
TMAPS['birth_year'] = TensorMap('22200_Year-of-birth_0_0', group='continuous', channel_map={'22200_Year-of-birth_0_0': 0}, normalization = {'mean': 1952.0639129359386, 'std': 7.656326148519739}, annotation_units=1, loss='logcosh', loss_weight=1.0)
TMAPS['birth_year_34'] = TensorMap('34_Year-of-birth_0_0', group='continuous', channel_map={'34_Year-of-birth_0_0': 0}, normalization = {'mean': 1952.0639129359386, 'std': 7.656326148519739}, annotation_units=1, loss='logcosh', loss_weight=1.0)
TMAPS['brain_volume'] = TensorMap('25010_Volume-of-brain-greywhite-matter_2_0', group='continuous', normalization={'mean': 1165940.0, 'std': 111511.0}, channel_map={'25010_Volume-of-brain-greywhite-matter_2_0': 0}, loss='logcosh', loss_weight=0.1)

TMAPS['sodium'] = TensorMap('30530_Sodium-in-urine', group='continuous', channel_map={'30530_Sodium-in-urine_0_0': 0},
                            normalization={'mean': 77.45323967267045, 'std': 44.441236848463774}, annotation_units=1, loss='logcosh')
TMAPS['potassium'] = TensorMap('30520_Potassium-in-urine', group='continuous', channel_map={'30520_Potassium-in-urine_0_0': 0},
                               normalization={'mean': 63.06182700345117, 'std': 33.84208704773539}, annotation_units=1, loss='logcosh')
TMAPS['cholesterol_hdl'] = TensorMap('30760_HDL-cholesterol', group='continuous', channel_map={'30760_HDL-cholesterol_0_0': 0},
                                     normalization={'mean': 1.4480129055069355, 'std': 0.3823115953478376}, annotation_units=1, loss='logcosh')
TMAPS['cholesterol'] = TensorMap('30690_Cholesterol', group='continuous', channel_map={'30690_Cholesterol_0_0': 0},
                                 normalization={'mean': 5.692381214399044, 'std': 1.1449409331668705}, annotation_units=1, loss='logcosh')

TMAPS['cigarettes'] = TensorMap('2887_Number-of-cigarettes-previously-smoked-daily_0_0', group='continuous', channel_map={'2887_Number-of-cigarettes-previously-smoked-daily_0_0': 0}, normalization = {'mean': 18.92662147068755, 'std':10.590930376362259 }, annotation_units=1)
TMAPS['alcohol'] = TensorMap('5364_Average-weekly-intake-of-other-alcoholic-drinks_0_0', group='continuous', channel_map={'5364_Average-weekly-intake-of-other-alcoholic-drinks_0_0': 0}, normalization = {'mean': 0.03852570253005904, 'std':0.512608370266108 }, annotation_units=1)
TMAPS['coffee'] = TensorMap('1498_Coffee-intake_0_0', group='continuous', channel_map={'1498_Coffee-intake_0_0': 0},
                            normalization={'mean': 2.015086529948216, 'std': 2.0914960998390497}, annotation_units=1)
TMAPS['water'] = TensorMap('1528_Water-intake_0_0', group='continuous', channel_map={'1528_Water-intake_0_0': 0},
                            normalization={'mean': 2.7322977785723324, 'std': 2.261996814128837}, annotation_units=1)
TMAPS['meat'] = TensorMap('3680_Age-when-last-ate-meat_0_0', group='continuous',
                            channel_map={'3680_Age-when-last-ate-meat_0_0': 0},
                            normalization={'mean': 29.74062983480561, 'std': 14.417292213873964}, annotation_units=1)
TMAPS['walks'] = TensorMap('864_Number-of-daysweek-walked-10-minutes_0_0', group='continuous',
                           channel_map={'864_Number-of-daysweek-walked-10-minutes_0_0': 0},
                           normalization={'mean': 5.369732285440756, 'std': 1.9564911925721618}, annotation_units=1)
TMAPS['walk_duration'] = TensorMap('874_Duration-of-walks_0_0', group='continuous', channel_map={'874_Duration-of-walks_0_0': 0},
                           normalization={'mean': 61.64092215093373, 'std': 78.79522990818906}, annotation_units=1)
TMAPS['physical_activities'] = TensorMap('884_Number-of-daysweek-of-moderate-physical-activity-10-minutes_0_0', group='continuous',
                           channel_map={'884_Number-of-daysweek-of-moderate-physical-activity-10-minutes_0_0': 0 },
                           normalization={'mean': 3.6258833281089258, 'std': 2.3343738999823676}, annotation_units=1)
TMAPS['physical_activity'] = TensorMap('894_Duration-of-moderate-activity_0_0', group='continuous',
                           channel_map={'894_Duration-of-moderate-activity_0_0': 0 },
                           normalization={'mean': 66.2862593866103, 'std': 77.28681218835422}, annotation_units=1)
TMAPS['physical_activity_vigorous'] = TensorMap('904_Number-of-daysweek-of-vigorous-physical-activity-10-minutes_0_0', group='continuous',
                           channel_map={'904_Number-of-daysweek-of-vigorous-physical-activity-10-minutes_0_0': 0},
                           normalization={'mean': 1.838718301735063, 'std': 1.9593505421480895}, annotation_units=1)
TMAPS['physical_activity_vigorous_duration'] = TensorMap('914_Duration-of-vigorous-activity_0_0', group='continuous',
                           channel_map={'914_Duration-of-vigorous-activity_0_0': 0},
                           normalization={'mean': 44.854488382965144, 'std': 48.159967071781466}, annotation_units=1)
TMAPS['tv'] = TensorMap('1070_Time-spent-watching-television-TV_0_0', group='continuous',
                            channel_map={'1070_Time-spent-watching-television-TV_0_0': 0},
                            normalization={'mean': 2.7753595642790914, 'std': 1.7135478462887321}, annotation_units=1)
TMAPS['computer'] = TensorMap('1080_Time-spent-using-computer_0_0', group='continuous',
                            channel_map={'1080_Time-spent-using-computer_0_0': 0},
                            normalization={'mean': 0.9781465855433753, 'std': 1.4444414103121512}, annotation_units=1)
TMAPS['car'] = TensorMap('1090_Time-spent-driving_0_0', group='continuous', channel_map={'1090_Time-spent-driving_0_0': 0},
                            normalization={'mean': 0.8219851505445748, 'std': 1.304094814200189}, annotation_units=1)
TMAPS['summer'] = TensorMap('1050_Time-spend-outdoors-in-summer_0_0', group='continuous',
                            channel_map={'1050_Time-spend-outdoors-in-summer_0_0': 0},
                            normalization={'mean': 3.774492304870845, 'std': 2.430483731404539}, annotation_units=1)
TMAPS['winter'] = TensorMap('1060_Time-spent-outdoors-in-winter_0_0', group='continuous',
                            channel_map={'1060_Time-spent-outdoors-in-winter_0_0': 0},
                            normalization={'mean': 1.8629686916635555, 'std': 1.88916218603397}, annotation_units=1)

TMAPS['systolic_blood_pressure_0'] = TensorMap('4080_Systolic-blood-pressure-automated-reading_0_0', group='continuous', loss='logcosh',
                                               channel_map={'4080_Systolic-blood-pressure-automated-reading_0_0': 0},
                                               normalization={'mean': 137.79964191990328, 'std': 19.292863700283757})

TMAPS['diastolic_blood_pressure_0'] = TensorMap('4079_Diastolic-blood-pressure-automated-reading', group='continuous', loss='logcosh',
                                                channel_map={'4079_Diastolic-blood-pressure-automated-reading': 0},
                                                normalization={'mean': 82.20657551284782, 'std': 10.496040770224475})

# example of multi-field-continuous tensor map (note shape will be 1x8 to accommodate a not-missing channel for each value
# normalization must be dictionary of [mean, stdev] for each value. Requries an imputation method.
TMAPS['blood-pressure'] = TensorMap('blood-pressure', group='multi_field_continuous',
                          channel_map={'4080_Systolic-blood-pressure-automated-reading_0_0': 0,
                                       '4080_Systolic-blood-pressure-automated-reading_0_1': 1,
                                       '4079_Diastolic-blood-pressure-automated-reading_0_0': 2,
                                       '4079_Diastolic-blood-pressure-automated-reading_0_1': 3},
                          annotation_units=4,
                          normalization={'4080_Systolic-blood-pressure-automated-reading_0_0': [137.79964191990328, 19.292863700283757],
                                         '4080_Systolic-blood-pressure-automated-reading_0_1': [137.79964191990328, 19.292863700283757],
                                         '4079_Diastolic-blood-pressure-automated-reading_0_0': [82.20657551284782, 10.496040770224475],
                                         '4079_Diastolic-blood-pressure-automated-reading_0_1': [82.20657551284782, 10.496040770224475]},
                          imputation=IMPUTATION_RANDOM)

TMAPS['random-forest-fields'] = TensorMap('random-forest-fields', group='categorical',
                                  channel_map={'Medication-for-pain-relief-constipation-heartburn_Aspirin': 0,
                                               'Medication-for-pain-relief-constipation-heartburn_Do-not-know': 1,
                                               'Medication-for-pain-relief-constipation-heartburn_Ibuprofen-eg-Nurofen': 2,
                                               'Medication-for-pain-relief-constipation-heartburn_Laxatives-eg-Dulcolax-Senokot': 3,
                                               'Medication-for-pain-relief-constipation-heartburn_None-of-the-above': 4,
                                               'Medication-for-pain-relief-constipation-heartburn_Omeprazole-eg-Zanprol': 5,
                                               'Medication-for-pain-relief-constipation-heartburn_Paracetamol': 6,
                                               'Medication-for-pain-relief-constipation-heartburn_Ranitidine-eg-Zantac': 7,
                                               'Vascularheart-problems-diagnosed-by-doctor_None-of-the-above': 8,
                                               'Vascularheart-problems-diagnosed-by-doctor_Heart-attack': 9,
                                               'Vascularheart-problems-diagnosed-by-doctor_Angina': 10,
                                               'Vascularheart-problems-diagnosed-by-doctor_Stroke': 11,
                                               'Vascularheart-problems-diagnosed-by-doctor_High-blood-pressure': 12,
                                               'Vascularheart-problems-diagnosed-by-doctor_Prefer-not-to-answer': 13,
                                               'Had-other-major-operations_Yes--you-will-be-asked-about-this-later-by-an-interviewer':14,
                                               'Had-other-major-operations_No': 15,
                                               'Had-other-major-operations_Do-not-know': 16,
                                               'Had-other-major-operations_Prefer-not-to-answer': 17,
                                               'Sex_Female': 18,
                                               'Sex_Male': 19,
                                               'Mother-still-alive_Yes': 20,
                                               'Mother-still-alive_No': 21,
                                               'Mother-still-alive_Do-not-know': 22,
                                               'Mother-still-alive_Prefer-not-to-answer': 23,
                                               'Father-still-alive_Yes': 24,
                                               'Father-still-alive_No': 25,
                                               'Father-still-alive_Do-not-know': 26,
                                               'Father-still-alive_Prefer-not-to-answer': 27,
                                               'Adopted-as-a-child_Yes': 28,
                                               'Adopted-as-a-child_No': 29,
                                               'Adopted-as-a-child_Do-not-know': 30,
                                               'Adopted-as-a-child_Prefer-not-to-answer': 31
                                               })

TMAPS['categorical-phenotypes-25'] = TensorMap(
    'categorical-phenotypes-25', group='categorical',
    channel_map={'Adopted-as-a-child_No_0_0': 0,
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

                 })

TMAPS['categorical-phenotypes-36'] = TensorMap(
    'categorical-phenotypes-36', group='categorical',
    channel_map={'Adopted-as-a-child_No_0_0': 0,
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
                 'Wheeze-or-whistling-in-the-chest-in-last-year_No_0_0': 35
                 })

TMAPS['categorical-phenotypes-78'] = TensorMap(
    'categorical-phenotypes-78', group='categorical', annotation_units=64,
    channel_map={'Adopted-as-a-child_No_0_0': 0,
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
                 'Wheeze-or-whistling-in-the-chest-in-last-year_No_0_0': 77})

TMAPS['categorical-phenotypes-134'] = TensorMap(
        'categorical-phenotypes-134',
        group='categorical', annotation_units=64,
        channel_map={'Adopted-as-a-child_No_0_0': 0, 'Adopted-as-a-child_No_2_0': 1,
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
                     'Workplace-very-dusty_Rarelynever_0_0': 127, 'Workplace-very-dusty_Rarelynever_0_1': 128})
