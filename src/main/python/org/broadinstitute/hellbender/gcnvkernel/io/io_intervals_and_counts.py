import logging
import numpy as np
import pandas as pd
from typing import Optional, List, Tuple, Set

from . import io_commons
from . import io_consts
from .. import types
from ..structs.interval import Interval, interval_annotations_dtypes, interval_annotations_dict

_logger = logging.getLogger(__name__)

interval_dtypes_dict = {
    io_consts.contig_column_name: np.str,
    io_consts.start_column_name: types.med_uint,
    io_consts.end_column_name: types.med_uint
}

read_count_dtypes_dict = {
    **interval_dtypes_dict,
    io_consts.count_column_name: types.med_uint
}


def load_read_counts_tsv_file(read_counts_tsv_file: str,
                              max_rows: Optional[int] = None,
                              return_interval_list: bool = False,
                              comment=io_consts.default_comment_char,
                              delimiter=io_consts.default_delimiter_char) \
        -> Tuple[str, np.ndarray, Optional[List[Interval]]]:
    """Loads a read counts .tsv file.

    Args:
        read_counts_tsv_file: input read counts .tsv file
        max_rows: (optional) maximum number of rows to process
        return_interval_list: if true, an interval list will also be generated and returned
        delimiter: delimiter character
        comment: comment character

    Returns:
        sample name, counts, (and optionally a list of intervals if `return_interval_list` == True)
    """
    sample_name = io_commons.extract_sample_name_from_header(read_counts_tsv_file)
    counts_pd = pd.read_csv(read_counts_tsv_file, delimiter=delimiter, comment=comment, nrows=max_rows,
                            dtype={**read_count_dtypes_dict})
    if return_interval_list:
        interval_list_pd = counts_pd[list(interval_dtypes_dict.keys())]
        interval_list = _convert_interval_list_pandas_to_gcnv_interval_list(interval_list_pd, read_counts_tsv_file)
        return sample_name, counts_pd[io_consts.count_column_name].as_matrix(), interval_list
    else:
        return sample_name, counts_pd[io_consts.count_column_name].as_matrix(), None


def load_interval_list_tsv_file(interval_list_tsv_file: str,
                                comment=io_consts.default_comment_char,
                                delimiter=io_consts.default_delimiter_char) -> List[Interval]:
    """Loads an interval list .tsv file.
    Args:
        interval_list_tsv_file: input interval list .tsv file
        delimiter: delimiter character
        comment: comment character

    Returns:
        interval list
    """
    interval_list_pd = pd.read_csv(interval_list_tsv_file, delimiter=delimiter, comment=comment,
                                   dtype={**interval_dtypes_dict, **interval_annotations_dtypes})
    return _convert_interval_list_pandas_to_gcnv_interval_list(interval_list_pd, interval_list_tsv_file)


def extract_sam_header_from_file(input_file: str):
    """Extract SAM header from a file.

    Notes:
        Only contiguous SAM header lines (starting with '@') are considered. The parsing of the input file
        stops as soon as a line starting with any other character is reached.

    Returns:
        a list of str
    """
    sam_header_list: List[str] = list()
    with open(input_file, 'r') as f:
        for line in f:
            stripped_line = line.strip()
            if len(stripped_line) == 0:
                continue
            elif stripped_line[0] == '@':
                sam_header_list.append(stripped_line)
            else:
                break
    return sam_header_list


def load_counts_in_the_modeling_zone(read_count_file_list: List[str],
                                     modeling_interval_list: List[Interval]) -> Tuple[List[str], np.ndarray]:
    """Loads read counts for a given cohort corresponding to a provided list of intervals.

    Args:
        read_count_file_list: list of read counts .tsv files
        modeling_interval_list: requested list of intervals

    Raises:
        AssertionError: if some of the intervals in `modeling_interval_list` are absent in the
        provided read counts .tsv file

    Note:
        it is assumed that all read counts have the SAME intervals.
        this assumption is not asserted for speed.

    Returns:
        list of sample names, 2-dim (sample x interval) ndarray of read counts
    """
    num_intervals = len(modeling_interval_list)
    num_samples = len(read_count_file_list)
    assert num_samples > 0
    assert num_intervals > 0

    sample_names: List[str] = []
    n_st = np.zeros((num_samples, num_intervals), dtype=types.med_uint)
    master_interval_list = None
    interval_to_index_map = None
    for si, read_count_file in enumerate(read_count_file_list):
        if master_interval_list is None:  # load intervals from the first read counts table
            sample_name, n_t, master_interval_list = load_read_counts_tsv_file(
                read_count_file, return_interval_list=True)
            interval_to_index_map = {interval: ti for ti, interval in enumerate(master_interval_list)}
            assert all([interval in interval_to_index_map for interval in modeling_interval_list]), \
                "Some of the modeling intervals are absent in the provided read counts .tsv file"
        else:  # do not load intervals again for speed, assume it is the same as the first sample
            sample_name, n_t, _ = load_read_counts_tsv_file(read_count_file, return_interval_list=False)
        # subset the counts in the order dictated by modeling_interval_list
        n_st[si, :] = np.asarray([n_t[interval_to_index_map[interval]]
                                  for interval in modeling_interval_list], dtype=types.med_uint)
        sample_names.append(sample_name)
    return sample_names, n_st


def _convert_interval_list_pandas_to_gcnv_interval_list(interval_list_pd: pd.DataFrame,
                                                        input_tsv_file: str) -> List[Interval]:
    """Converts a pandas dataframe of intervals to list(Interval). Annotations will be parsed
    and added to the intervals as well.

    Args:
        interval_list_pd: interval list as a pandas dataframe
        input_tsv_file: path to the .tsv file associated to the dataframe
            (only used to generate exception messages)

    Returns:
        a list of intervals
    """
    interval_list: List[Interval] = list()
    columns = {str(x) for x in interval_list_pd.columns.values}
    io_commons.assert_mandatory_columns(set(interval_dtypes_dict.keys()), columns, input_tsv_file)
    for contig, start, end in zip(interval_list_pd[io_consts.contig_column_name],
                                  interval_list_pd[io_consts.start_column_name],
                                  interval_list_pd[io_consts.end_column_name]):
        interval = Interval(contig, start, end)
        interval_list.append(interval)

    found_annotation_keys: Set[str] = columns.intersection(interval_annotations_dict.keys())
    if len(found_annotation_keys) > 0:
        _logger.info("The given interval list provides the following interval annotations: "
                     + str(found_annotation_keys))
        for annotation_key in found_annotation_keys:
            bad_annotations_found = False
            for ti, raw_value in enumerate(interval_list_pd[annotation_key]):
                try:
                    annotation = interval_annotations_dict[annotation_key](raw_value)
                    interval_list[ti].add_annotation(annotation_key, annotation)
                except ValueError:
                    bad_annotations_found = True
            if bad_annotations_found:
                _logger.warning("Some of the annotations for {0} contained bad values and were ignored".format(
                    annotation_key))

    return interval_list


def write_interval_list_to_tsv_file(output_file: str, interval_list: List[Interval],
                                    sam_header_lines: Optional[List[str]] = None):
    """Write a list of interval list to .tsv file.

    Note:
        If all intervals have an annotation, that annotation will be written to the .tsv files.
        If only some intervals have an annotation, that annotation will be ignored and a warning
        will be logged.

    Args:
        output_file: output .tsv file
        interval_list: list of intervals to write to .tsv file
        sam_header_lines: (optional) SAM header lines

    Returns:
        None
    """
    assert len(interval_list) > 0, "Can not write an empty interval list to disk"
    annotation_found_keys: Set[str] = set()
    for interval in interval_list:
        for key in interval.annotations.keys():
            annotation_found_keys.add(key)
    mutual_annotation_key_list: List[str] = []
    for key in annotation_found_keys:
        if all(key in interval.annotations.keys() for interval in interval_list):
            mutual_annotation_key_list.append(key)
        else:
            _logger.warning("Only a subset of intervals contain annotation \"{0}\"; "
                            "Cannot write this annotation to a .tsv file; "
                            "Neglecting \"{0}\" and proceeding...".format(key))
    with open(output_file, 'w') as out:
        if sam_header_lines is not None:
            for sam_header_line in sam_header_lines:
                out.write(sam_header_line + '\n')
        header = '\t'.join([io_consts.contig_column_name,
                            io_consts.start_column_name,
                            io_consts.end_column_name]
                           + mutual_annotation_key_list)
        out.write(header + '\n')
        for interval in interval_list:
            row = '\t'.join([interval.contig, repr(interval.start), repr(interval.end)] +
                            [str(interval.annotations[key]) for key in mutual_annotation_key_list])
            out.write(row + '\n')
