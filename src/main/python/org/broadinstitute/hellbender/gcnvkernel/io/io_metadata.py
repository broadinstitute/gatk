import logging
import numpy as np
import os
import pandas as pd
from typing import List
from collections import OrderedDict

from . import io_commons
from . import io_consts
from .. import types
from ..structs.metadata import SampleReadDepthMetadata, SamplePloidyMetadata, SampleCoverageMetadata, \
    SampleMetadataCollection

_logger = logging.getLogger(__name__)


def read_sample_coverage_metadata(sample_metadata_collection: SampleMetadataCollection,
                                  input_files: List[str],
                                  comment=io_consts.default_comment_char,
                                  delimiter=io_consts.default_delimiter_char) -> List[str]:
    """Reads sample coverage metadata from .tsv files and adds them to `sample_metadata_collection`.

    Args:
        sample_metadata_collection: collection to which the coverage metadata is to be added
        input_files: input sample coverage metadata .tsv file
        comment: comment character
        delimiter: delimiter character

    Returns:
        list of samples in the same order as encountered in `input_files`
    """
    sample_names = []
    max_count = None
    contig_list = []
    for sample_index, input_file in enumerate(input_files):
        coverage_metadata_pd = pd.read_csv(input_file, delimiter=delimiter, comment=comment)
        found_columns_list = [str(column) for column in coverage_metadata_pd.columns.values]
        io_commons.assert_mandatory_columns({io_consts.contig_column_name}, set(found_columns_list), input_file)
        count_columns = found_columns_list.copy()
        count_columns.remove(io_consts.contig_column_name)
        if max_count is None:
            max_count = len(count_columns) - 1
            contig_list = coverage_metadata_pd[io_consts.contig_column_name].tolist()
        else:
            assert len(count_columns) - 1 == max_count, \
                "Maximum count in per-contig count distribution file \"{0}\" " \
                "does not match that in other files.".format(input_file)
            assert coverage_metadata_pd[io_consts.contig_column_name].tolist() == contig_list, \
                "Contigs {0} in per-contig count distribution file \"{1}\" " \
                "do not match those in other files.".format(coverage_metadata_pd[io_consts.contig_column_name].tolist(), input_file)
        sample_name = io_commons.extract_sample_name_from_header(input_file)
        sample_names.append(sample_name)
        contig_hist_m = OrderedDict((contig, np.asarray(coverage_metadata_pd.loc[contig_index, count_columns],
                                                        dtype=types.med_uint))
                                    for contig_index, contig in enumerate(contig_list))
        sample_metadata_collection.add_sample_coverage_metadata(SampleCoverageMetadata(
            sample_name, contig_hist_m))

    return sample_names


def update_sample_metadata_collection_from_ploidy_determination_calls(
        sample_metadata_collection: SampleMetadataCollection,
        input_calls_path: str,
        comment=io_consts.default_comment_char,
        delimiter=io_consts.default_delimiter_char):
    """Reads the output of contig ploidy determination tool and updates the given instance of
    `SampleMetadataCollection` for read depth and ploidy metadata.

    Args:
        sample_metadata_collection: the instance of `SampleMetadataCollection` to be updated
        input_calls_path: posterior output path of contig ploidy determination tool
        comment: comment character
        delimiter: delimiter character

    Returns:
        None
    """

    def get_sample_read_depth_metadata(input_path: str) -> SampleReadDepthMetadata:
        sample_read_depth_file = os.path.join(input_path, io_consts.default_sample_read_depth_tsv_filename)
        assert os.path.exists(sample_read_depth_file), \
            "Sample read depth could not be found in the contig ploidy results " \
            "located at \"{0}\"".format(input_path)

        _sample_name = io_commons.extract_sample_name_from_header(sample_read_depth_file)

        sample_read_depth_pd = pd.read_csv(sample_read_depth_file, delimiter=delimiter, comment=comment)

        io_commons.assert_mandatory_columns(
            SampleReadDepthMetadata.mandatory_tsv_columns,
            {str(column) for column in sample_read_depth_pd.columns.values},
            sample_read_depth_file)

        global_read_depth = sample_read_depth_pd[io_consts.global_read_depth_column_name].values[0]
        average_ploidy = sample_read_depth_pd[io_consts.average_ploidy_column_name].values[0]

        return SampleReadDepthMetadata(_sample_name, global_read_depth, average_ploidy)

    def get_sample_ploidy_metadata(input_path: str) -> SamplePloidyMetadata:
        sample_ploidy_file = os.path.join(input_path, io_consts.default_sample_contig_ploidy_tsv_filename)
        assert os.path.exists(sample_ploidy_file), \
            "Sample ploidy results could not be found in the contig ploidy results " \
            "located at \"{0}\"".format(input_path)

        _sample_name = io_commons.extract_sample_name_from_header(sample_ploidy_file)

        sample_ploidy_pd = pd.read_csv(sample_ploidy_file, delimiter=delimiter, comment=comment)

        io_commons.assert_mandatory_columns(
            SamplePloidyMetadata.mandatory_tsv_columns,
            {str(column) for column in sample_ploidy_pd.columns.values},
            sample_ploidy_file)

        contig_list = [str(x) for x in sample_ploidy_pd[io_consts.contig_column_name].values]
        ploidy_list = [int(x) for x in sample_ploidy_pd[io_consts.ploidy_column_name].values]
        ploidy_gq_list = [float(x) for x in sample_ploidy_pd[io_consts.ploidy_gq_column_name].values]

        return SamplePloidyMetadata(_sample_name,
                                    np.asarray(ploidy_list, dtype=types.small_uint),
                                    np.asarray(ploidy_gq_list, dtype=types.floatX),
                                    contig_list)

    _logger.info("Loading germline contig ploidy and global read depth metadata...")
    assert os.path.exists(input_calls_path) and os.path.isdir(input_calls_path), \
        "The provided path to ploidy determination results \"{0}\" is not a directory".format(input_calls_path)

    subdirs = os.listdir(input_calls_path)
    for subdir in subdirs:
        if subdir.find(io_consts.sample_folder_prefix) >= 0:
            sample_ploidy_results_dir = os.path.join(input_calls_path, subdir)

            sample_name = io_commons.get_sample_name_from_txt_file(sample_ploidy_results_dir)
            sample_read_depth_metadata = get_sample_read_depth_metadata(sample_ploidy_results_dir)
            sample_ploidy_metadata = get_sample_ploidy_metadata(sample_ploidy_results_dir)

            assert (sample_read_depth_metadata.sample_name == sample_name and
                    sample_ploidy_metadata.sample_name == sample_name), \
                "Inconsistency detected in the ploidy determination results in {0}: sample name in the .txt " \
                "file does not match with sample name in the posterior headers".format(sample_ploidy_results_dir)

            sample_metadata_collection.add_sample_read_depth_metadata(sample_read_depth_metadata)
            sample_metadata_collection.add_sample_ploidy_metadata(sample_ploidy_metadata)
