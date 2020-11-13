from ml4h.TensorMap import TensorMap, Interpretation
from ml4h.defines import StorageType
from ml4h.metrics import weighted_crossentropy


diploid_cm = {'homozygous_reference': 0, 'heterozygous': 1, 'homozygous_variant': 2}
rs3829740 = TensorMap('rs3829740', Interpretation.CATEGORICAL, channel_map=diploid_cm)
rs2234962 = TensorMap('rs2234962', Interpretation.CATEGORICAL, channel_map=diploid_cm)
rs2042995 = TensorMap('rs2042995', Interpretation.CATEGORICAL, channel_map=diploid_cm)

rs3829740_weighted = TensorMap('rs3829740', Interpretation.CATEGORICAL, channel_map=diploid_cm, loss=weighted_crossentropy([1, 1, 1.5], 'rs3829740'))
rs2234962_weighted = TensorMap('rs2234962', Interpretation.CATEGORICAL, channel_map=diploid_cm, loss=weighted_crossentropy([.8, 1, 1.5], 'rs2234962'))
rs2042995_weighted = TensorMap('rs2042995', Interpretation.CATEGORICAL, channel_map=diploid_cm, loss=weighted_crossentropy([.6, 1.5, 2], 'rs2042995'))


akap9_lof = TensorMap('AKAP9', Interpretation.CATEGORICAL, channel_map={'no_akap9_lof': 0, 'akap9_lof': 1})
dsc2_lof = TensorMap('DSC2', Interpretation.CATEGORICAL, channel_map={'no_dsc2_lof': 0, 'dsc2_lof': 1})
ryr2_lof = TensorMap('RYR2', Interpretation.CATEGORICAL, channel_map={'no_ryr2_lof': 0, 'ryr2_lof': 1})
ttn_lof = TensorMap('TTN', Interpretation.CATEGORICAL, channel_map={'no_ttn_lof': 0, 'ttn_lof': 1})


def _ttn_tensor_from_file(tm, hd5, dependents={}):
    index = 0
    categorical_data = np.zeros(tm.shape, dtype=np.float32)
    if 'has_exome' not in hd5['categorical']:
        raise ValueError('Skipping people without exome sequencing.')
    if tm.name in hd5['categorical'] and int(hd5['categorical'][tm.name][0]) != 0:
        index = 1
    categorical_data[index] = 1.0
    return categorical_data


ttntv = TensorMap(
    'has_ttntv',  Interpretation.CATEGORICAL, channel_map={
    'no_TTN_tv': 0, 'TTN_tv': 1,
    }, tensor_from_file=_ttn_tensor_from_file,
)
ttntv_10x = TensorMap(
    'has_ttntv',  Interpretation.CATEGORICAL, channel_map={
    'no_TTN_tv': 0, 'TTN_tv': 1,
    }, loss_weight=10.0, tensor_from_file=_ttn_tensor_from_file,
)


bsa_mosteller = TensorMap('bsa_mosteller',  Interpretation.CONTINUOUS, normalization={'mean': 1.8894831981880114, 'std': 0.22169301057810176}, loss='logcosh', channel_map={'bsa_mosteller': 0})
bsa_dubois = TensorMap('bsa_dubois',  Interpretation.CONTINUOUS, normalization={'mean': 1.8671809970639703, 'std': 0.20913930961120797}, loss='logcosh', channel_map={'bsa_dubois': 0})




genetic_pca_1 = TensorMap(
    '22009_Genetic-principal-components_0_1', Interpretation.CONTINUOUS, path_prefix='continuous', normalization={'mean': -0.014422761536727896, 'std': 10.57799283718005},
    loss='logcosh', channel_map={'22009_Genetic-principal-components_0_1': 0},
)
genetic_pca_2 = TensorMap(
    '22009_Genetic-principal-components_0_2', Interpretation.CONTINUOUS, path_prefix='continuous', #normalization={'mean': -0.014422761536727896, 'std': 10.57799283718005},
    loss='logcosh', activation='linear', channel_map={'22009_Genetic-principal-components_0_2': 0},
)
genetic_pca_3 = TensorMap(
    '22009_Genetic-principal-components_0_3', Interpretation.CONTINUOUS, path_prefix='continuous', #normalization={'mean': -0.014422761536727896, 'std': 10.57799283718005},
    loss='logcosh', activation='linear', channel_map={'22009_Genetic-principal-components_0_3': 0},
)
genetic_pca_4 = TensorMap(
    '22009_Genetic-principal-components_0_4', Interpretation.CONTINUOUS, path_prefix='continuous', #normalization={'mean': -0.014422761536727896, 'std': 10.57799283718005},
    loss='logcosh', activation='linear', channel_map={'22009_Genetic-principal-components_0_4': 0},
)
genetic_pca_5 = TensorMap(
    '22009_Genetic-principal-components_0_5', Interpretation.CONTINUOUS, path_prefix='continuous', #normalization={'mean': -0.014422761536727896, 'std': 10.57799283718005},
    loss='logcosh', activation='linear', channel_map={'22009_Genetic-principal-components_0_5': 0},
)
genetic_pca_all5 = TensorMap(
    'genetic_pca_all5', Interpretation.CONTINUOUS, path_prefix='continuous', normalization={'mean': -0.014422761536727896, 'std': 10.57799283718005},
    loss='logcosh', annotation_units=5, shape=(5,), activation='linear',
    channel_map={
        '22009_Genetic-principal-components_0_0': 0, '22009_Genetic-principal-components_0_1': 1,
        '22009_Genetic-principal-components_0_2': 2, '22009_Genetic-principal-components_0_3': 3,
        '22009_Genetic-principal-components_0_4': 4,
    },
)

genetic_caucasian = TensorMap(
    'Genetic-ethnic-grouping_Caucasian_0_0', Interpretation.CATEGORICAL, path_prefix='categorical', storage_type=StorageType.CATEGORICAL_FLAG,
    channel_map={'no_caucasian': 0, 'Genetic-ethnic-grouping_Caucasian_0_0': 1})

genetic_caucasian_weighted = TensorMap(
    'Genetic-ethnic-grouping_Caucasian_0_0', Interpretation.CATEGORICAL, path_prefix='categorical', storage_type=StorageType.CATEGORICAL_FLAG,
    channel_map={'no_caucasian': 0, 'Genetic-ethnic-grouping_Caucasian_0_0': 1}, loss=weighted_crossentropy([10.0, 1.0], 'caucasian_loss'),
)
