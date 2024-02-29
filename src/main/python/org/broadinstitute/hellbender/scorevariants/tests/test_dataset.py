from unittest import TestCase
from scorevariants.tests.utilities import (
    MockGATKTensorDirectory,
)
from scorevariants.dataset import (
    GATKTensorDataset,
)
from torch.utils.data import DataLoader
import numpy as np


class GATK2DTensorDatasetTests(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.check_keys = [
            "type",
            "label",
            "false_positive",
            "chrom",
            "pos",
            "best_practices",
            "read_tensor",
        ]
        cls.mock_directory = MockGATKTensorDirectory(tensors_per_directory=5)
        cls.dataset = GATKTensorDataset(
            cls.mock_directory.name, in_memory=False
        )
        cls.in_memory_dataset = GATKTensorDataset(
            cls.mock_directory.name, in_memory=True
        )
        cls.list_instantiated_dataset = GATKTensorDataset(
            cls.mock_directory.filenames[:16],
        )
        cls.datasets = [
            cls.dataset,
            cls.in_memory_dataset,
            cls.list_instantiated_dataset,
        ]
        cls.dataset_lengths = [
            20,
            20,
            16,
        ]

    def test_dataset_length(self):
        for dataset, exp_len in zip(self.datasets, self.dataset_lengths):
            with self.subTest():
                dataset_length = len(dataset)
                self.assertEqual(exp_len, dataset_length)

    def test_dataset_getitem_int_index(self):
        for dataset in self.datasets:
            with self.subTest():
                item = dataset[0]
                self.assertCountEqual(self.check_keys, item.keys())
                self.assertIsInstance(item["type"], float)
                self.assertIsInstance(item["label"], float)
                self.assertIsInstance(item["false_positive"], float)
                self.assertIsInstance(item["chrom"], str)
                self.assertIsInstance(item["pos"], np.int64)
                self.assertTupleEqual((7,), item["best_practices"].shape)
                self.assertTupleEqual(
                    (128, 128, 15), item["read_tensor"].shape
                )

    def test_dataset_getitem_slice_index(self):
        for dataset in self.datasets:
            with self.subTest():
                item = dataset[1:4]
                self.assertCountEqual(self.check_keys, item.keys())
                self.assertTupleEqual((3,), item["type"].shape)
                self.assertTupleEqual((3,), item["label"].shape)
                self.assertTupleEqual((3,), item["false_positive"].shape)
                self.assertEqual(3, len(item["chrom"]))
                self.assertTupleEqual((3,), item["pos"].shape)
                self.assertTupleEqual(
                    (3, 7),
                    item["best_practices"].shape,
                )
                self.assertTupleEqual(
                    (3, 128, 128, 15), item["read_tensor"].shape
                )

    def test_data_loader_integration(self):
        for dataset in self.datasets:
            with self.subTest():
                dataloader = DataLoader(dataset, batch_size=16)
                item = next(iter(dataloader))
                self.assertCountEqual(self.check_keys, item.keys())
                self.assertTupleEqual((16,), item["type"].shape)
                self.assertTupleEqual((16,), item["label"].shape)
                self.assertTupleEqual((16,), item["false_positive"].shape)
                self.assertEqual(16, len(item["chrom"]))
                self.assertTupleEqual((16,), item["pos"].shape)
                self.assertTupleEqual(
                    (
                        16,
                        7,
                    ),
                    item["best_practices"].shape,
                )
                self.assertTupleEqual(
                    (16, 128, 128, 15), item["read_tensor"].shape
                )


class GATK1DTensorDatasetTests(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mock_directory = MockGATKTensorDirectory(
            tensors_per_directory=5,
            mode="1d",
        )
        cls.dataset = GATKTensorDataset(
            cls.mock_directory.name,
            in_memory=False,
        )
        cls.in_memory_dataset = GATKTensorDataset(
            cls.mock_directory.name, in_memory=True
        )
        cls.list_instantiated_dataset = GATKTensorDataset(
            cls.mock_directory.filenames[:16],
        )
        cls.datasets = [
            cls.dataset,
            cls.in_memory_dataset,
            cls.list_instantiated_dataset,
        ]
        cls.dataset_lengths = [
            20,
            20,
            16,
        ]
        cls.check_keys = [
            "type",
            "label",
            "false_positive",
            "chrom",
            "pos",
            "best_practices",
            "reference",
        ]

    def test_dataset_length(self):
        for dataset, exp_len in zip(self.datasets, self.dataset_lengths):
            with self.subTest():
                dataset_length = len(dataset)
                self.assertEqual(exp_len, dataset_length)

    def test_dataset_getitem_int_index(self):
        for dataset in self.datasets:
            with self.subTest():
                item = dataset[0]
                self.assertCountEqual(self.check_keys, item.keys())
                self.assertIsInstance(item["type"], float)
                self.assertIsInstance(item["label"], float)
                self.assertIsInstance(item["false_positive"], float)
                self.assertIsInstance(item["chrom"], str)
                self.assertIsInstance(item["pos"], np.int64)
                self.assertTupleEqual((7,), item["best_practices"].shape)
                self.assertTupleEqual((128, 4), item["reference"].shape)

    def test_dataset_getitem_slice_index(self):
        for dataset in self.datasets:
            with self.subTest():
                item = dataset[1:4]
                self.assertCountEqual(self.check_keys, item.keys())
                self.assertTupleEqual((3,), item["type"].shape)
                self.assertTupleEqual((3,), item["label"].shape)
                self.assertTupleEqual((3,), item["false_positive"].shape)
                self.assertEqual(3, len(item["chrom"]))
                self.assertTupleEqual((3,), item["pos"].shape)
                self.assertTupleEqual(
                    (
                        3,
                        7,
                    ),
                    item["best_practices"].shape,
                )
                self.assertTupleEqual((3, 128, 4), item["reference"].shape)

    def test_data_loader_integration(self):
        for dataset in self.datasets:
            with self.subTest():
                dataloader = DataLoader(dataset, batch_size=16)
                item = next(iter(dataloader))
                self.assertCountEqual(self.check_keys, item.keys())
                self.assertTupleEqual((16,), item["type"].shape)
                self.assertTupleEqual((16,), item["label"].shape)
                self.assertTupleEqual((16,), item["false_positive"].shape)
                self.assertEqual(16, len(item["chrom"]))
                self.assertTupleEqual((16,), item["pos"].shape)
                self.assertTupleEqual(
                    (
                        16,
                        7,
                    ),
                    item["best_practices"].shape,
                )
                self.assertTupleEqual((16, 128, 4), item["reference"].shape)
