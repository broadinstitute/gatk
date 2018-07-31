import numpy as np
import logging
from typing import List, Set

from ..structs.interval import Interval
from .. import types

__all__ = ['IntervalListMask']

_logger = logging.getLogger(__name__)


class IntervalListMask:
    """This class implements different filters ("masks") on interval lists."""
    def __init__(self, interval_list: List[Interval]):
        self.num_intervals = len(interval_list)
        self.interval_list = interval_list
        self.drop_t = np.zeros((self.num_intervals,), dtype=bool)
        self.drop_reason_t: List[Set[str]] = [set() for _ in range(self.num_intervals)]

    def _assert_mask_compatibility_with_read_count_array(self, n_st: np.ndarray):
        assert n_st.shape[1] == self.num_intervals, \
            "Number of intervals in the mask instance ({0}) is not compatible with the shape of the provided " \
            "read count array (shape = {1})".format(self.num_intervals, n_st.shape)

    @staticmethod
    def _assert_read_count_int_dtype(n_st: np.ndarray):
        assert n_st.dtype in types.int_dtypes or n_st.dtype in types.uint_dtypes, \
            "Can not reliably detect cohort-wide uncovered intervals with the dtype of the given " \
            "read counts array ({0})".format(n_st.dtype)

    def get_masked_view(self, n_st: np.ndarray):
        """Applies the mask on a given interval list and read count array.

        Args:
            n_st: a read count matrix

        Returns:
            a view of the provided n_st,
            a new list containing references to the provided interval list
        """
        self._assert_mask_compatibility_with_read_count_array(n_st)
        kept_intervals_indices = [ti for ti in range(len(self.interval_list)) if not self.drop_t[ti]]
        num_dropped_intervals = self.num_intervals - len(kept_intervals_indices)
        kept_interval_list = [self.interval_list[ti] for ti in kept_intervals_indices]
        kept_n_st = n_st[:, kept_intervals_indices]
        if num_dropped_intervals > 0:
            dropped_fraction = num_dropped_intervals / self.num_intervals
            _logger.warning("Some intervals were dropped. Dropped fraction: {0:2.2}".format(dropped_fraction))
        return kept_n_st, kept_interval_list

    def keep_only_given_contigs(self, contigs_to_keep: Set[str]):
        """Only keep intervals in the given contig and mask the rest.

        Args:
            contigs_to_keep: set of contigs to keep

        Returns:
            None
        """
        inactive_interval_indices = [interval.contig not in contigs_to_keep for interval in self.interval_list]
        self.drop_t[inactive_interval_indices] = True
        for ti in inactive_interval_indices:
            self.drop_reason_t[ti].add("contig marked to be dropped")

    def drop_blacklisted_intervals(self, blacklisted_intervals: List[Interval]):
        """Mask any interval that overlaps with a list of blacklisted intervals.

        Args:
            blacklisted_intervals: list of blacklisted intervals.

        Returns:
            None
        """
        for ti, interval in enumerate(self.interval_list):
            if any([interval.overlaps_with(interval) for interval in blacklisted_intervals]):
                self.drop_t[ti] = True
                self.drop_reason_t[ti].add("blacklisted")

    def drop_cohort_wide_uncovered_intervals(self, n_st: np.ndarray):
        """Mask any interval that has zero coverage for all samples in a given read count matrix.

        Args:
            n_st: read count matrix

        Returns:
            None
        """
        self._assert_mask_compatibility_with_read_count_array(n_st)
        self._assert_read_count_int_dtype(n_st)
        for ti in range(self.num_intervals):
            if all(n_st[:, ti] == 0):
                self.drop_t[ti] = True
                self.drop_reason_t[ti].add("cohort-wide uncovered interval")
