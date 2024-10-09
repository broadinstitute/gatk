from torch.utils.data import Dataset
from scorevariants.readers import ReferenceTensorReader
import numpy as np

class ReferenceDataset(Dataset):
    """
    Dataset that creates tensors on-the-fly from a reader

    """

    def __init__(self, reference_reader: ReferenceTensorReader):
        """
        Args:
            reference_reader: an intialized reader for the reference tensors
        """
        self.reference_reader = reference_reader
        self._index_variants()

    def _index_variants(self):
        """
        Index the variants in the dataset
        """
        self._index = np.array(self.reference_reader.get_variants())

    def __len__(self):
        return len(self._index)

    def __getitem__(self, idx):
        """
        Fetches an item(s) from the dataset

        Returns:
            A dict containing the following keys:
            * chrom : str for the chromosome of the variant
            * pos : int for the position of the variant(s) (1-based)
            * ref : str for the reference sequence of the variant
            * alt : str for the alt alleles
            * type : torch.long for the variant type (SNP/INDEL)
            * best_practices: torch.float32 for the annotation vector
            * reference: torch.float32 for the reference tensor

        """
        variant = self._index[idx]
        return self.reference_reader.variant_representation(variant)
