from typing import Optional

from .. import config
from ..io import io_consts


class IntegerCopyNumberSegment:
    """Represents a constant copy-number genomic interval ("segment") and stores various quality metrics.

    Note:
        Quality metric attributes are considered optional.

    Attributes:
        contig: contig name
        start: segment start genomic position
        end: segment end genomic position
        num_points: number of points (e.g. bins or intervals) contained in the segment
        call_copy_number: segment most-likely copy-number call
        baseline_copy_number: baseline copy-number for the contig
        quality_some_called: complementary phred-scaled probably that one or more of the points in the segment
            agree with the segment copy-number call (normalized by `num_points`)
        quality_all_called: complementary phred-scaled probability that all of points in the segment
            agree with the segment copy-number call
        quality_start: complementary phred-scaled probability that the leftmost point of the segment
            is a change-point to the segment copy-number call
        quality_end: complementary phred-scaled probability that the rightmost point of the segment
            is a change-point to the segment copy-number call
    """
    def __init__(self, contig: str, start: int, end: int,
                 num_points: int,
                 call_copy_number: int,
                 baseline_copy_number: int):
        self.contig: str = contig
        self.start: int = start
        self.end: int = end
        self.num_points: int = num_points
        self.call_copy_number: int = call_copy_number
        self.baseline_copy_number = baseline_copy_number
        self.quality_some_called: Optional[float] = None
        self.quality_all_called: Optional[float] = None
        self.quality_start: Optional[float] = None
        self.quality_end: Optional[float] = None

    @staticmethod
    def get_header_column_string():
        return '\t'.join([io_consts.contig_column_name,
                          io_consts.start_column_name,
                          io_consts.end_column_name,
                          io_consts.num_points_column_name,
                          io_consts.call_copy_number_column_name,
                          io_consts.baseline_copy_number_column_name,
                          io_consts.quality_some_called_column_name,
                          io_consts.quality_all_called_column_name,
                          io_consts.quality_start_column_name,
                          io_consts.quality_end_column_name])

    @staticmethod
    def _repr_quality(quality):
        return '{0:.{1}f}'.format(quality, config.phred_decimals) if quality is not None else '.'

    def __repr__(self):
        return '\t'.join([self.contig,
                          repr(self.start),
                          repr(self.end),
                          repr(self.num_points),
                          repr(self.call_copy_number),
                          repr(self.baseline_copy_number),
                          self._repr_quality(self.quality_some_called),
                          self._repr_quality(self.quality_all_called),
                          self._repr_quality(self.quality_start),
                          self._repr_quality(self.quality_end)])
