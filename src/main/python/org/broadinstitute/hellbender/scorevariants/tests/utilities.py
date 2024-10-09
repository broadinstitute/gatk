from io import BytesIO, StringIO
import h5py
import numpy as np
import os

from tempfile import TemporaryDirectory
from typing import Dict


class MockBAMTensor:
    """
    Mocks BAM tensors (2-D CNN input)
    """

    def __call__(self):
        return np.zeros((128, 128, 15))


class MockBestPractices:
    """
    Mocks Best practices (MLP input)
    """

    def __call__(self):
        return np.zeros(7)


class MockRefTensor:
    """
    Mocks Refereence Tensors (1-D CNN input)
    """

    def __call__(self):
        return np.zeros((128, 4))


class _MockGATKTensorHDF5:
    def _init_shared(self):
        self.mock_best_practices = MockBestPractices()()

    def _write_shared(self, f):
        f["best_practices"] = self.mock_best_practices


class Mock2DTensorHDF5(_MockGATKTensorHDF5):
    """
    Mocks the 2-D and best practices inputs and stored for a single variant
    """

    def __init__(self):
        self.mock_bam = MockBAMTensor()()
        self._init_shared()

    def __call__(self):
        bytes_io = BytesIO()
        with h5py.File(bytes_io, "w") as f:
            f["read_tensor"] = self.mock_bam
            self._write_shared(f)

        return bytes_io


class Mock1DTensorHDF5(_MockGATKTensorHDF5):
    """
    Mocks the 1-D and best practices inputs stored for a single variant
    """

    def __init__(self):
        self.mock_ref = MockRefTensor()()
        self._init_shared()

    def __call__(self):
        bytes_io = BytesIO()
        with h5py.File(bytes_io, "w") as f:
            f["reference"] = self.mock_ref
            self._write_shared(f)

        return bytes_io


class MockGATKTensorDirectory:
    """
    Mocks a GATKs tensor directory structure, namely:

        ├── INDEL
        ├── NOT_INDEL
        ├── NOT_SNP
        └── SNP

    Where each directory stored HDF5 files corresponding
    to a single variant.
    """

    class FileNameGenerator:
        def __init__(self, suffix=None):
            self.suffix = suffix if suffix else ""
            self.n = 0

        def __call__(self):
            output = self.n
            self.n += 1
            return f"{output}-chr{output}_{output}" + self.suffix

    directories = ["SNP", "NOT_SNP", "INDEL", "NOT_INDEL"]

    def __init__(
        self,
        mode="2d",
        tensors_per_directory=None,
        append_directory=None,
    ):
        self.mode = mode
        if append_directory is None:
            append_directory = ""
        self.dir = TemporaryDirectory()

        if tensors_per_directory is None:
            number_of_files = [10, 10, 10, 10]
        elif isinstance(tensors_per_directory, int):
            number_of_files = 4 * [tensors_per_directory]
        else:
            number_of_files = tensors_per_directory

        filename_gen = self.FileNameGenerator(".hd5")
        if mode == "2d":
            tensor_gen = Mock2DTensorHDF5()
        elif mode == "1d":
            tensor_gen = Mock1DTensorHDF5()
        else:
            raise ValueError(f"Invalid Choice for 'mode': '{mode}'")

        self.filenames = []
        for directory, number in zip(self.directories, number_of_files):
            full_dir = os.path.join(self.dir.name, append_directory, directory)
            os.makedirs(full_dir, exist_ok=True)
            for i in range(number):
                filename = os.path.join(full_dir, filename_gen())
                tensor_hdf5 = tensor_gen()
                with open(filename, "wb") as fh:
                    fh.write(tensor_hdf5.getvalue())
                self.filenames.append(filename)
        self.name = self.dir.name

    def __del__(self):
        self.dir.cleanup()


test_fasta = """>SEQUENCE_1
ACGTTCG
>SEQUENCE_2
TGCATGC
>SEQUENCE_3
NGCRTGC
"""


def _make_test_fasta(fasta=None):
    if fasta is None:
        fasta = test_fasta
    fasta_io = StringIO(fasta)
    return fasta_io


class MockVariantRecord:
    def __init__(
        self,
        pos: int,
        contig: str = None,
        info: Dict = None,
        ref=None,
        alts=None,
        chrom=None,
        *args,
    ):
        """
        Attributes:
            pos: position in the reference. 1-based inclusive.
            contig: name of the contig
        """
        self.pos = pos
        self.contig = contig
        if info is None:
            info = dict()
        self.info = info
        self.ref = ref
        self.alts = alts
        if chrom is None:
            chrom = contig
        self.chrom = chrom
