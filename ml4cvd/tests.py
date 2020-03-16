import unittest

from tensorflow.keras.losses import logcosh

from ml4cvd.TensorMap import TensorMap
from ml4cvd.arguments import parse_args
from ml4cvd.tensor_maps_by_script import TMAPS
from ml4cvd.recipes import test_multimodal_multitask, train_multimodal_multitask


ALL_TENSORS = '/mnt/ml4cvd/projects/tensors/sax-lax-ecg-rest-brain-1k/2019-11-06/'
ALL_TENSORS = '/mnt/disks/sax-lax-40k-lvm/2020-01-29/'
MODELS = '/mnt/ml4cvd/projects/models/for_testing/'


def _run_tests():
    suites = []
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestTensorMaps))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestTrainingModels))
    suites.append(unittest.TestLoader().loadTestsFromTestCase(TestPretrainedModels))
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
        delta = 2e-1
        args = parse_args()
        args.tensors = ALL_TENSORS
        args.input_tensors = ['categorical-phenotypes-134']
        args.output_tensors = ['coronary_artery_disease_soft', 'diabetes_type_2',
                               'hypertension', 'myocardial_infarction']
        args.epochs = 2
        args.optimizer = 'radam'
        args.batch_size = 64
        args.training_steps = 60
        args.validation_steps = 1
        args.test_steps = 64
        args.tensor_maps_in = [TMAPS[it] for it in args.input_tensors]
        args.tensor_maps_out = [TMAPS[ot] for ot in args.output_tensors]
        performances = train_multimodal_multitask(args)
        print('cat mlp expected = ', performances)
        expected = {'no_coronary_artery_disease_soft': 0.4945240833284755, 'coronary_artery_disease_soft': 0.4945240833284755,
                    'no_diabetes_type_2': 0.5470901959163001, 'diabetes_type_2': 0.5470910572382917, 'no_hypertension': 0.4933869347607721,
                    'hypertension': 0.49338678981369954, 'no_myocardial_infarction': 0.5612242431079916, 'myocardial_infarction': 0.5612242431079916}

        for k in performances:
            self.assertAlmostEqual(performances[k], expected[k], delta=delta)


class TestPretrainedModels(unittest.TestCase):
    def test_ecg_regress(self):
        delta = 15e-2
        args = parse_args()
        args.tensors = ALL_TENSORS
        args.optimizer = 'radam'
        args.model_file = MODELS + 'ecg_rest_regress.hd5'
        args.input_tensors = ['ecg_rest']
        args.output_tensors = ['p-axis', 'p-duration', 'p-offset', 'p-onset', 'pp-interval', 'pq-interval', 'q-offset', 'q-onset', 'qrs-complexes',
                               'qrs-duration', 'qt-interval', 'qtc-interval', 'r-axis', 'rr-interval', 't-offset', 't-axis']
        args.test_steps = 16
        args.batch_size = 24
        args.tensor_maps_in = [TMAPS[it] for it in args.input_tensors]
        args.tensor_maps_out = [TMAPS[ot] for ot in args.output_tensors]
        performances = test_multimodal_multitask(args)
        print('expected = ', performances)
        expected = {'PAxis_pearson': 0.6010829310005574, 'PDuration_pearson': 0.334397562979115, 'POffset_pearson': 0.8870224086564694,
                    'POnset_pearson': 0.915672533398803, 'PPInterval_pearson': 0.9843003466847959, 'PQInterval_pearson': 0.9024483040343809,
                    'QOffset_pearson': 0.7529739057566074, 'QOnset_pearson': 0.7417930219416248, 'QRSComplexes_pearson': 0.14711527713690073,
                    'QRSDuration_pearson': 0.9000439573470225, 'QTInterval_pearson': 0.9572200766158885, 'QTCInterval_pearson': 0.919422124977152,
                    'RAxis_pearson': 0.8858920261655899, 'RRInterval_pearson': 0.985806047245532, 'TOffset_pearson': 0.9446387492651126,
                    'TAxis_pearson': 0.6731483067923112}

        for k in expected:
            self.assertAlmostEqual(performances[k], expected[k], delta=delta)

    def test_ecg_rhythm(self):
        delta = 1e-1
        args = parse_args()
        args.optimizer = 'radam'
        args.tensors = ALL_TENSORS
        args.model_file = MODELS + 'ecg_rest_rhythm.hd5'
        args.input_tensors = ['ecg_rest']
        args.output_tensors = ['ecg_rhythm_poor']
        args.test_steps = 32
        args.batch_size = 32
        args.tensor_maps_in = [TMAPS[it] for it in args.input_tensors]
        args.tensor_maps_out = [TMAPS[ot] for ot in args.output_tensors]
        performances = test_multimodal_multitask(args)
        print('expected = ', performances)
        expected = {'Normal_sinus_rhythm': 0.9824029566392662, 'Sinus_bradycardia': 0.9918234160230809, 'Marked_sinus_bradycardia': 0.9975368393759994,
                    'Other_sinus_rhythm': 0.8898630136986301, 'Atrial_fibrillation': 0.9869892718557407, 'Other_rhythm': 0.7223864258347017}

        for k in expected:
            self.assertAlmostEqual(performances[k], expected[k], delta=delta)


# Back to the top!
if '__main__' == __name__:
    _run_tests()
