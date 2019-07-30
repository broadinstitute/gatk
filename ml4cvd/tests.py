import unittest

from keras.losses import logcosh

from ml4cvd.TensorMap import TensorMap
from ml4cvd.arguments import parse_args
from ml4cvd.tensor_maps_by_script import TMAPS
from ml4cvd.DatabaseClient import BigQueryDatabaseClient, SqLiteDatabaseClient
from ml4cvd.recipes import test_multimodal_multitask, train_multimodal_multitask

ALL_TENSORS = '/mnt/disks/data/generated/tensors/test/2019-03-21/'
MODELS = '/mnt/ml4cvd/projects/jamesp/data/models/'


def _run_tests():
    suites = []
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestTensorMaps))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestTrainingModels))
    # TODO Add them to 'suites' when the pretrained model tests are updated to work with the tensors at ALL_TENSORS
    # suites.append(unittest.TestLoader().loadTestsFromTestCase(TestPretrainedModels))
    # suites.append(unittest.TestLoader().loadTestsFromTestCase(TestDatabaseClient))

    unittest.TextTestRunner(verbosity=3).run(unittest.TestSuite(suites))


class TestTensorMaps(unittest.TestCase):

    def test_tensor_map_equality(self):
        tensor_map_1a = TensorMap(name='tm', loss='logcosh', channel_map={'c1': 1, 'c2': 2}, metrics=[])
        tensor_map_1b = TensorMap(name='tm', loss='logcosh', channel_map={'c1': 1, 'c2': 2}, metrics=[])
        tensor_map_2a = TensorMap(name='tm', loss=logcosh, channel_map={'c1': 1, 'c2': 2}, metrics=[])
        tensor_map_2b = TensorMap(name='tm', loss=logcosh, channel_map={'c2': 2, 'c1': 1}, metrics=[])
        tensor_map_3 = TensorMap(name='tm', loss=logcosh, channel_map={'c1': 1, 'c2': 3}, metrics=[])
        tensor_map_4 = TensorMap(name='tm', loss=logcosh, channel_map={'c1': 1, 'c2': 3}, metrics=[all])
        tensor_map_5a = TensorMap(name='tm', loss=logcosh, channel_map={'c1': 1, 'c2': 3}, metrics=[all, any])
        tensor_map_5b = TensorMap(name='tm', loss=logcosh, channel_map={'c1': 1, 'c2': 3}, metrics=[any, all])
        tensor_map_6a = TensorMap(name='tm', loss=logcosh, channel_map={'c1': 1, 'c2': 3}, dependent_map=tensor_map_1a)
        tensor_map_6b = TensorMap(name='tm', loss=logcosh, channel_map={'c1': 1, 'c2': 3}, dependent_map=tensor_map_1b)

        self.assertEqual(tensor_map_1a, tensor_map_1b)
        self.assertEqual(tensor_map_2a, tensor_map_2b)
        self.assertEqual(tensor_map_1a, tensor_map_2a)
        self.assertNotEqual(tensor_map_2a, tensor_map_3)
        self.assertNotEqual(tensor_map_3, tensor_map_4)
        self.assertNotEqual(tensor_map_3, tensor_map_5a)
        self.assertNotEqual(tensor_map_4, tensor_map_5a)
        self.assertEqual(tensor_map_5a, tensor_map_5b)
        self.assertEqual(tensor_map_6a, tensor_map_6b)


class TestTrainingModels(unittest.TestCase):

    def test_train_categorical_mlp(self):
        delta = 1e-1
        args = parse_args()
        args.tensors = ALL_TENSORS
        args.input_tensors = ['categorical-phenotypes-78']
        args.output_tensors = ['coronary_artery_disease_soft', 'diabetes_type_2',
                               'hypertension', 'myocardial_infarction']
        args.epochs = 1
        args.batch_size = 32
        args.training_steps = 20
        args.validation_steps = 1
        args.test_steps = 32
        args.tensor_maps_in = [TMAPS[it] for it in args.input_tensors]
        args.tensor_maps_out = [TMAPS[ot] for ot in args.output_tensors]
        performances = train_multimodal_multitask(args)
        print('expected = ', performances)
        expected = {'no_coronary_artery_disease_soft': 0.528143258213825, 'coronary_artery_disease_soft': 0.528143258213825,
                    'no_diabetes_type_2': 0.6547365677800461, 'diabetes_type_2': 0.654736567780046, 'no_hypertension': 0.4729961761211761,
                    'hypertension': 0.4729961761211761, 'no_myocardial_infarction': 0.5480460307260938, 'myocardial_infarction': 0.5480460307260938}

        for k in performances:
            self.assertAlmostEqual(performances[k], expected[k], delta=delta)

    def test_train_mlp_cat36_pi(self):
        delta = 8e-1
        args = parse_args()
        args.tensors = ALL_TENSORS
        args.input_tensors = ['categorical-phenotypes-36']
        args.output_tensors = ['end_systole_volume', 'end_diastole_volume', 'ejection_fractionp',
                               'allergic_rhinitis_prevalent_incident', 'asthma_prevalent_incident',
                               'atrial_fibrillation_or_flutter_prevalent_incident', 'back_pain_prevalent_incident',
                               'breast_cancer_prevalent_incident', 'coronary_artery_disease_soft_prevalent_incident',
                               'diabetes_type_2_prevalent_incident',
                               'hypertension_prevalent_incident', 'myocardial_infarction_prevalent_incident']
        args.epochs = 1
        args.batch_size = 32
        args.training_steps = 20
        args.validation_steps = 1
        args.test_steps = 32
        args.tensor_maps_in = [TMAPS[it] for it in args.input_tensors]
        args.tensor_maps_out = [TMAPS[ot] for ot in args.output_tensors]
        performances = train_multimodal_multitask(args)
        print('expected = ', performances)
        expected = {'end_systole_volume_pearson': -0.18898451656690726,
                    'end_diastole_volume_pearson': 0.21810526919810447,
                    'ejection_fraction_pearson': 0.10870165257050499, 'no_allergic_rhinitis': 0.6144818514407859,
                    'prevalent_allergic_rhinitis': 0.7167318982387475, 'incident_allergic_rhinitis': 0.6180602006688962,
                    'no_asthma': 0.63458251953125, 'prevalent_asthma': 0.39872051795899494,
                    'incident_asthma': 0.6231338522393773, 'no_atrial_fibrillation_or_flutter': 0.5410833333333334,
                    'prevalent_atrial_fibrillation_or_flutter': 0.4638504611330698,
                    'incident_atrial_fibrillation_or_flutter': 0.633399209486166, 'no_back_pain': 0.5106227106227106,
                    'prevalent_back_pain': 0.5311572700296736, 'incident_back_pain': 0.6132478632478633,
                    'no_breast_cancer': 0.4686281065418478, 'prevalent_breast_cancer': 0.41631089217296113,
                    'incident_breast_cancer': 0.5118179810028718, 'no_coronary_artery_disease_soft': 0.527491408934708,
                    'prevalent_coronary_artery_disease_soft': 0.6513790600616949,
                    'incident_coronary_artery_disease_soft': 0.43841355846774194,
                    'no_diabetes_type_2': 0.5623635418471669, 'prevalent_diabetes_type_2': 0.48804846103470856,
                    'incident_diabetes_type_2': 0.5779872301611432, 'no_hypertension': 0.4664909906701802,
                    'prevalent_hypertension': 0.5800677998245739, 'incident_hypertension': 0.47476802434951704,
                    'no_myocardial_infarction': 0.6280876494023904,
                    'prevalent_myocardial_infarction': 0.7055550569864488,
                    'incident_myocardial_infarction': 0.5212917350848385}

        for k in performances:
            self.assertAlmostEqual(performances[k], expected[k], delta=delta)

    def test_train_mri_sax_zoom(self):
        delta = 7e-1
        args = parse_args()
        args.tensors = ALL_TENSORS
        args.input_tensors = ['sax_inlinevf_zoom_weighted']
        args.output_tensors = ['sax_inlinevf_zoom_mask_weighted', 'end_systole_volume', 'end_diastole_volume',
                               'ejection_fractionp', 'allergic_rhinitis', 'asthma', 'atrial_fibrillation_or_flutter',
                               'back_pain', 'breast_cancer', 'coronary_artery_disease_soft', 'diabetes_type_2',
                               'hypertension', 'myocardial_infarction']
        args.epochs = 1
        args.batch_size = 4
        args.training_steps = 24
        args.validation_steps = 1
        args.test_steps = 36
        args.t = 48
        args.pool_z = 2
        args.u_connect = True
        args.learning_rate = 0.0001
        args.tensor_maps_in = [TMAPS[it] for it in args.input_tensors]
        args.tensor_maps_out = [TMAPS[ot] for ot in args.output_tensors]
        performances = train_multimodal_multitask(args)
        print('expected = ', performances)
        expected = {'end_systole_volume_pearson': 0.010191148347492063,
                    'end_diastole_volume_pearson': 0.07746212713601273,
                    'ejection_fraction_pearson': 0.039482962469710864, 'no_allergic_rhinitis': 0.6583710407239819,
                    'allergic_rhinitis': 0.655697243932538, 'no_asthma': 0.4921536796536796,
                    'asthma': 0.49607683982683987, 'no_atrial_fibrillation_or_flutter': 0.3950320512820513,
                    'atrial_fibrillation_or_flutter': 0.3950320512820512, 'no_back_pain': 0.5713760117733627,
                    'back_pain': 0.5673289183222958, 'no_breast_cancer': 0.32195723684210525,
                    'breast_cancer': 0.3223684210526316, 'no_coronary_artery_disease_soft': 0.37866666666666665,
                    'coronary_artery_disease_soft': 0.37733333333333335, 'no_diabetes_type_2': 0.5410830999066294,
                    'diabetes_type_2': 0.5401493930905695, 'no_hypertension': 0.5034782608695653,
                    'hypertension': 0.5039613526570048, 'no_myocardial_infarction': 0.5632911392405063,
                    'myocardial_infarction': 0.564873417721519}

        for k in expected:
            self.assertAlmostEqual(performances[k], expected[k], delta=delta)

    def test_train_mri_systole_diastole(self):
        delta = 6e-1
        args = parse_args()
        args.tensors = ALL_TENSORS
        args.input_tensors = ['mri_systole_diastole_weighted']
        args.output_tensors = ['mri_systole_diastole_segmented_weighted', 'end_systole_volume', 'end_diastole_volume',
                               'ejection_fractionp', 'allergic_rhinitis', 'asthma', 'atrial_fibrillation_or_flutter',
                               'back_pain', 'breast_cancer', 'coronary_artery_disease_soft', 'diabetes_type_2',
                               'hypertension', 'myocardial_infarction']
        args.epochs = 1
        args.batch_size = 12
        args.training_steps = 96
        args.validation_steps = 1
        args.test_steps = 36
        args.u_connect = True
        args.learning_rate = 0.0005
        args.tensor_maps_in = [TMAPS[it] for it in args.input_tensors]
        args.tensor_maps_out = [TMAPS[ot] for ot in args.output_tensors]
        performances = train_multimodal_multitask(args)
        print('expected = ', performances)
        expected = {'background': 0.9997983811584907, 'ventricle': 0.3717467509311169, 'myocardium': 0.1753433358365212,
                    'end_systole_volume_pearson': 0.1376, 'end_diastole_volume_pearson': 0.13085844389045515,
                    'ejection_fraction_pearson': 0.09538719107146973, 'no_allergic_rhinitis': 0.46750762240688715,
                    'allergic_rhinitis': 0.4675076224068871, 'no_asthma': 0.5619121188712681,
                    'asthma': 0.5619121188712682, 'no_atrial_fibrillation_or_flutter': 0.6406823580220254,
                    'atrial_fibrillation_or_flutter': 0.6406823580220254, 'no_back_pain': 0.5318302387267905,
                    'back_pain': 0.5318302387267904, 'no_breast_cancer': 0.5274939903846154,
                    'breast_cancer': 0.5274939903846154, 'no_coronary_artery_disease_soft': 0.589689578713969,
                    'coronary_artery_disease_soft': 0.589689578713969, 'no_diabetes_type_2': 0.6338928856914469,
                    'diabetes_type_2': 0.6338928856914469, 'no_hypertension': 0.5819672131147541,
                    'hypertension': 0.5819672131147541, 'no_myocardial_infarction': 0.6285789335434726,
                    'myocardial_infarction': 0.6285789335434726}

        for k in expected:
            self.assertAlmostEqual(performances[k], expected[k], delta=delta)

    def test_train_mri_systole_diastole_pi(self):
        delta = 6e-1
        args = parse_args()
        args.tensors = ALL_TENSORS
        args.input_tensors = ['mri_systole_diastole_weighted']
        args.output_tensors = ['mri_systole_diastole_segmented_weighted', 'end_systole_volume', 'end_diastole_volume',
                               'ejection_fractionp', 'allergic_rhinitis_prevalent_incident',
                               'asthma_prevalent_incident', 'atrial_fibrillation_or_flutter_prevalent_incident',
                               'back_pain_prevalent_incident', 'breast_cancer_prevalent_incident',
                               'coronary_artery_disease_soft_prevalent_incident',
                               'diabetes_type_2_prevalent_incident', 'hypertension_prevalent_incident',
                               'myocardial_infarction_prevalent_incident']
        args.epochs = 1
        args.batch_size = 12
        args.training_steps = 96
        args.validation_steps = 1
        args.test_steps = 36
        args.u_connect = True
        args.learning_rate = 0.0005
        args.tensor_maps_in = [TMAPS[it] for it in args.input_tensors]
        args.tensor_maps_out = [TMAPS[ot] for ot in args.output_tensors]
        performances = train_multimodal_multitask(args)
        print('expected = ', performances)
        expected = {'background': 0.9994056020189077, 'ventricle': 0.3217267098772595, 'myocardium': 0.13751608768060483,
                    'end_systole_volume_pearson': -0.11803115474509973, 'end_diastole_volume_pearson': 0.07460391201856462,
                    'ejection_fraction_pearson': 0.014315858148787145, 'no_allergic_rhinitis': 0.4270341364261374,
                    'prevalent_allergic_rhinitis': 0.5918604651162791, 'incident_allergic_rhinitis': 0.437839186576009,
                    'no_asthma': 0.4421774889807788, 'prevalent_asthma': 0.4578408195429472, 'incident_asthma': 0.4612041884816754,
                    'no_atrial_fibrillation_or_flutter': 0.33988339451522354, 'prevalent_atrial_fibrillation_or_flutter': 0.3615023474178404,
                    'incident_atrial_fibrillation_or_flutter': 0.2651053864168618, 'no_back_pain': 0.5578817733990148,
                    'prevalent_back_pain': 0.307563025210084, 'incident_back_pain': 0.553459920988913, 'no_breast_cancer': 0.49939903846153844,
                    'prevalent_breast_cancer': 0.6088992974238876, 'incident_breast_cancer': 0.45130641330166266,
                    'no_coronary_artery_disease_soft': 0.4605321507760532, 'prevalent_coronary_artery_disease_soft': 0.3188862621486735,
                    'incident_coronary_artery_disease_soft': 0.49899026987332473, 'no_diabetes_type_2': 0.47961630695443647,
                    'prevalent_diabetes_type_2': 0.5485625485625485, 'incident_diabetes_type_2': 0.4376984126984127,
                    'no_hypertension': 0.6255949233209941, 'prevalent_hypertension': 0.4478044259066156, 'incident_hypertension': 0.6184678890849811,
                    'no_myocardial_infarction': 0.5187811925400578, 'prevalent_myocardial_infarction': 0.5925058548009368,
                    'incident_myocardial_infarction': 0.4287383177570093}

        for k in expected:
            self.assertAlmostEqual(performances[k], expected[k], delta=delta)


class TestPretrainedModels(unittest.TestCase):

    def test_ecg(self):
        delta = 9e-2
        args = parse_args()
        args.tensors = ALL_TENSORS
        args.model_file = MODELS + 'ecg_regress_diagnose.hd5'
        args.input_tensors = ['ecg_rest']
        args.output_tensors = ['ecg_normal', 'ecg_rhythm', 'p-axis', 'p-duration', 'p-offset', 'p-onset',
                               'pp-interval', 'pq-interval', 'q-offset', 'q-onset', 'qrs-duration', 'qrs-num',
                               'qt-interval', 'qtc-interval', 'ventricular-rate']
        args.test_steps = 128
        args.tensor_maps_in = [TMAPS[it] for it in args.input_tensors]
        args.tensor_maps_out = [TMAPS[ot] for ot in args.output_tensors]
        performances = test_multimodal_multitask(args)
        expected = {'Normal-ECG': 0.8752554358212905, 'Borderline-ECG': 0.5837728637728636,
                    'Abnormal-ECG': 0.723338215389989, 'Normal-sinus-rhythm': 0.9741548715787861,
                    'Otherwise-normal-ECG': 0.827044714759781, 'Sinus-bradycardia': 0.979314510040767,
                    'Marked-sinus-bradycardia': 0.9928738897936085, 'Atrial-Fibrillation': nan,
                    'PAxis_pearson': 0.3176830509728512, 'PDuration_pearson': 0.5991046519033024,
                    'POffset_pearson': 0.7529158679935809, 'POnset_pearson': 0.7180879470117637,
                    'PPInterval_pearson': 0.9042360806052834, 'PQInterval_pearson': 0.6737001705381603,
                    'QOffset_pearson': 0.2785414608613505, 'QOnset_pearson': -0.005761179589827403,
                    'QRSDuration_pearson': 0.47948940700236614, 'QRSNum_pearson': 0.9050081284476049,
                    'QTInterval_pearson': 0.8430877545949877, 'QTCInterval_pearson': 0.6979122540251536,
                    'VentricularRate_pearson': 0.9454089108881661}
        for k in expected:
            self.assertAlmostEqual(performances[k], expected[k], delta=delta)

    def test_segmenter(self):
        delta = 5e-2
        args = parse_args()
        args.tensors = ALL_TENSORS
        args.model_file = MODELS + '/seg_efp_from_zoom_unet.hd5'
        args.input_tensors = ['sax_inlinevf_zoom']
        args.output_tensors = ['sax_inlinevf_zoom_mask', 'end_systole_volume', 'end_diastole_volume',
                               'ejection_fractionp']
        args.batch_size = 3
        args.test_steps = 64
        args.tensor_maps_in = [TMAPS[it] for it in args.input_tensors]
        args.tensor_maps_out = [TMAPS[ot] for ot in args.output_tensors]
        performances = test_multimodal_multitask(args)
        expected = {'end_systole_volume_pearson': 0.6230096872484032, 'end_diastole_volume_pearson': 0.6701106883874123,
                    'ejection_fraction_pearson': 0.18809344508798623}
        for k in performances:
            self.assertAlmostEqual(performances[k], expected[k], delta=delta)

    def test_continuous_mlp(self):
        delta = 8e-2
        args = parse_args()
        args.tensors = ALL_TENSORS
        args.model_file = MODELS + 'mlp_cont.hd5'
        args.input_tensors = ['34_0', '1807_0', '1845_0', '1873_0', '1883_0', '2946_0', '3526_0', '3972_0', '3982_0',
                              '5057_0', '2355_0', '3809_0', '2217_0', '4689_0', '4700_0', '5430_0', '5901_0', '5923_0',
                              '5945_0', '2966_0', '2976_0', '3627_0', '3761_0', '3786_0', '3894_0', '3992_0', '4012_0',
                              '4022_0', '4056_0', '4269_0', '4272_0', '4276_0', '4279_0', '1568_0', '1578_0', '1588_0',
                              '1598_0', '1608_0', '4407_0', '4418_0', '4429_0', '4440_0', '4451_0', '4462_0', '5364_0',
                              '1289_0', '1299_0', '1309_0', '1319_0', '1438_0', '1458_0', '1488_0', '1498_0', '1528_0',
                              '3680_0', '864_0', '874_0', '884_0', '894_0', '904_0', '914_0', '1070_0', '1080_0',
                              '1090_0', '1050_0', '1060_0', '1737_0', '2277_0', '1269_0', '1279_0', '2867_0', '2887_0',
                              '2897_0', '2926_0', '3436_0', '3456_0', '3659_0', '2684_0', '2704_0', '2714_0', '2734_0',
                              '2744_0', '2754_0', '2764_0', '2794_0', '2804_0', '2824_0', '3536_0', '3546_0', '3581_0',
                              '3700_0', '3710_0', '3829_0', '3839_0', '3849_0', '3872_0', '3882_0', '2405_0', '129_0',
                              '130_0', '84_0', '87_0', '134_0', '135_0', '137_0', '92_0', '136_0']
        args.output_tensors = ['allergic_rhinitis', 'anxiety', 'asthma', 'atrial_fibrillation_or_flutter', 'back_pain',
                               'breast_cancer', 'cardiac_surgery', 'cervical_cancer', 'colorectal_cancer',
                               'coronary_artery_disease_hard', 'coronary_artery_disease_intermediate',
                               'coronary_artery_disease_soft', 'death', 'diabetes_all', 'diabetes_type_1',
                               'diabetes_type_2', 'enlarged_prostate', 'heart_failure', 'hypertension', 'lung_cancer',
                               'migraine', 'myocardial_infarction', 'osteoporosis', 'skin_cancer', 'stroke']
        args.batch_size = 64
        args.test_steps = 64
        args.tensor_maps_in = [TMAPS[it] for it in args.input_tensors]
        args.tensor_maps_out = [TMAPS[ot] for ot in args.output_tensors]
        performances = test_multimodal_multitask(args)
        expected = {'no_allergic_rhinitis': 0.5381671339010805, 'allergic_rhinitis': 0.5381671339010806,
                    'no_anxiety': 0.7003621730382293, 'anxiety': 0.7003621730382293, 'no_asthma': 0.7098805328824278,
                    'asthma': 0.7098805328824278, 'no_atrial_fibrillation_or_flutter': 0.762897004838907,
                    'atrial_fibrillation_or_flutter': 0.7628970048389071, 'no_back_pain': 0.6392667065243514,
                    'back_pain': 0.6392674392023774, 'no_breast_cancer': 0.891462981524304,
                    'breast_cancer': 0.8914629815243041, 'no_cardiac_surgery': 0.8447070865386153,
                    'cardiac_surgery': 0.8447070865386154, 'no_cervical_cancer': 0.972783099975436,
                    'cervical_cancer': 0.9727830999754361, 'no_colorectal_cancer': 0.8080598255089322,
                    'colorectal_cancer': 0.8080598255089323, 'no_coronary_artery_disease_hard': 0.8916216598793527,
                    'coronary_artery_disease_hard': 0.8916216598793529,
                    'no_coronary_artery_disease_intermediate': 0.8589982315216934,
                    'coronary_artery_disease_intermediate': 0.8589982315216935,
                    'no_coronary_artery_disease_soft': 0.8524390592171021,
                    'coronary_artery_disease_soft': 0.8524390592171023, 'no_death': 0.7505171302818825,
                    'death': 0.7505171302818824, 'no_diabetes_all': 0.9006059334335557,
                    'diabetes_all': 0.9006064125080604, 'no_diabetes_type_1': 0.9038895730706076,
                    'diabetes_type_1': 0.903892993979201, 'no_diabetes_type_2': 0.9007381407791559,
                    'diabetes_type_2': 0.9007381407791558, 'no_enlarged_prostate': 0.8586038157110603,
                    'enlarged_prostate': 0.8586038157110603, 'no_heart_failure': 0.8363361497680583,
                    'heart_failure': 0.8363361497680583, 'no_hypertension': 0.7928368461934274,
                    'hypertension': 0.7928369827530641, 'no_lung_cancer': 0.7398319527229255,
                    'lung_cancer': 0.7398319527229255, 'no_migraine': 0.6591549718124025,
                    'migraine': 0.6591539722522091, 'no_myocardial_infarction': 0.9108621456071485,
                    'myocardial_infarction': 0.9108621456071485, 'no_osteoporosis': 0.7302523030717248,
                    'osteoporosis': 0.730252303071725, 'no_skin_cancer': 0.8265695508855826,
                    'skin_cancer': 0.8265695508855826, 'no_stroke': 0.8313994023904383, 'stroke': 0.8313994023904382}

        for k in performances:
            self.assertAlmostEqual(performances[k], expected[k], delta=delta)

    def test_categorical_mlp(self):
        delta = 5e-2
        args = parse_args()
        args.tensors = ALL_TENSORS
        args.model_file = MODELS + 'mlp_cat965.hd5'
        args.input_tensors = ['categorical-phenotypes-965']
        args.output_tensors = ['allergic_rhinitis', 'asthma', 'atrial_fibrillation_or_flutter', 'back_pain',
                               'breast_cancer', 'coronary_artery_disease_soft', 'death', 'diabetes_type_2',
                               'enlarged_prostate', 'hypertension', 'lung_cancer', 'multiple_sclerosis',
                               'myocardial_infarction', 'osteoporosis', 'skin_cancer', 'stroke']
        args.batch_size = 64
        args.test_steps = 128
        args.tensor_maps_in = [TMAPS[it] for it in args.input_tensors]
        args.tensor_maps_out = [TMAPS[ot] for ot in args.output_tensors]
        performances = test_multimodal_multitask(args)
        expected = {'no_allergic_rhinitis': 0.8845998401962242, 'allergic_rhinitis': 0.8845998401962241,
                    'no_asthma': 0.9644744885626663, 'asthma': 0.9644745560061379,
                    'no_atrial_fibrillation_or_flutter': 0.7925845922900356,
                    'atrial_fibrillation_or_flutter': 0.7925848104801057, 'no_back_pain': 0.7540731705229888,
                    'back_pain': 0.7540731705229888, 'no_breast_cancer': 0.9141909324845224,
                    'breast_cancer': 0.9141909324845222, 'no_coronary_artery_disease_soft': 0.8970631307079905,
                    'coronary_artery_disease_soft': 0.8970631307079905, 'no_death': 0.8167174118043976,
                    'death': 0.8167174118043976, 'no_diabetes_type_2': 0.9277685489342402,
                    'diabetes_type_2': 0.9277683800796117, 'no_enlarged_prostate': 0.8618313955559356,
                    'enlarged_prostate': 0.8618313955559355, 'no_hypertension': 0.9678282448775123,
                    'hypertension': 0.9678282448775122, 'no_lung_cancer': 0.8199246986269229,
                    'lung_cancer': 0.8199246986269231, 'no_multiple_sclerosis': 0.8969885208312303,
                    'multiple_sclerosis': 0.8969867181979436, 'no_myocardial_infarction': 0.9110085617783363,
                    'myocardial_infarction': 0.9110085617783364, 'no_osteoporosis': 0.8676444220599344,
                    'osteoporosis': 0.8676444220599345, 'no_skin_cancer': 0.806383671647947,
                    'skin_cancer': 0.806383671647947, 'no_stroke': 0.8978971165452445, 'stroke': 0.8978971165452445}
        for k in performances:
            self.assertAlmostEqual(performances[k], expected[k], delta=delta)


class TestDatabaseClient(unittest.TestCase):
    def test_database_client(self):
        args = parse_args()

        bigquery_client = BigQueryDatabaseClient(credentials_file=args.bigquery_credentials_file)
        sqlite_client = SqLiteDatabaseClient(db_file=args.db)

        query = (
            'SELECT field, fieldid FROM {} '
            'WHERE fieldid BETWEEN 3120 AND 3190 '
            'ORDER BY field ASC '
            'LIMIT 7'
        )

        table = 'dictionary'
        dataset = args.bigquery_dataset

        bigquery_rows = bigquery_client.execute(query.format(f"`{dataset}.{table}`"))
        sqlite_rows = sqlite_client.execute(query.format(table))

        for bq_row, sql_row in zip(bigquery_rows, sqlite_rows):
            self.assertEqual((bq_row[0], bq_row[1]), sql_row)


# Back to the top!
if '__main__' == __name__:
    _run_tests()
