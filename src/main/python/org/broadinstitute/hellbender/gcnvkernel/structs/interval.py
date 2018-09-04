import numpy as np
from abc import abstractmethod
from typing import Dict
from copy import deepcopy


class Interval:
    """This class represents a genomic interval along with optional annotations.

    Note:
        Equality test and hashing is based on get_key() which excludes all annotations.
    """
    def __init__(self, contig: str, start: int, end: int):
        assert end >= start, "Interval end point must be greater or equal to its start point"
        self.contig = str(contig)
        self.start = int(start)
        self.end = int(end)
        self.annotations = dict()
        self._hash = hash(self.get_key())

    def get_key(self):
        return self.contig, self.start, self.end

    def add_annotation(self, key: str, annotation: 'IntervalAnnotation'):
        self.annotations[key] = annotation

    def get_annotation(self, key: str):
        return self.annotations[key].get_value()

    def get_padded(self, padding: int, keep_annotations=False) -> 'Interval':
        assert padding >= 0, "padding must be >= 0"
        padded_interval = Interval(self.contig, self.start - padding, self.end + padding)
        if keep_annotations:
            padded_interval.annotations = deepcopy(self.annotations)
        return padded_interval

    def overlaps_with(self, other):
        return (self.contig == other.contig) and (min(self.end, other.end) - max(self.start, other.start) > 0)

    def get_midpoint(self):
        return 0.5 * (self.start + self.end)

    def distance(self, other):
        if self.contig != other.contig:
            return np.inf
        else:
            return np.abs(self.get_midpoint() - other.get_midpoint())

    def __eq__(self, other: 'Interval'):
        return self.get_key() == other.get_key()

    def __ne__(self, other: 'Interval'):
        return self.get_key() != other.get_key()

    def __hash__(self):
        return self._hash

    def __str__(self):
        return str(self.get_key())

    def __repr__(self):
        return self.__str__()


class IntervalAnnotation:
    """Base class for all interval annotations."""
    def __init__(self, raw_value):
        self.parsed_value = self.parse(raw_value)

    def get_value(self):
        return self.parsed_value

    @staticmethod
    @abstractmethod
    def parse(raw_value):
        """Takes a raw value (e.g. a value that is directly read from a .tsv file) and casts it to
        the correct type."""
        pass

    @staticmethod
    @abstractmethod
    def get_key() -> str:
        """Returns a string identifier key for the annotation. The key is used for reading the annotation from
        and writing the annotation to a .tsv file."""
        pass

    def __repr__(self):
        return repr(self.get_value())

    def __str__(self):
        return str(self.get_value())


class GCContentAnnotation(IntervalAnnotation):
    """This class represents GC content annotation for an interval."""
    def __init__(self, gc_content):
        super().__init__(gc_content)

    @staticmethod
    def parse(raw_value):
        gc_content = float(raw_value)
        if not 0.0 <= gc_content <= 1.0:
            raise ValueError("GC content ({0}) must be a float in range [0, 1]".format(gc_content))
        return gc_content

    @staticmethod
    def get_key():
        return "GC_CONTENT"


interval_annotations_dict: Dict[str, IntervalAnnotation] = {
    GCContentAnnotation.get_key(): GCContentAnnotation
}

interval_annotations_dtypes: Dict[str, object] = {
    GCContentAnnotation.get_key(): float
}
