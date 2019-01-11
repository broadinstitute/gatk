import csv
import logging
import numpy as np
import os
import pandas as pd
from typing import List

from . import io_commons
from . import io_consts
from .. import types
from ..structs.metadata import SampleReadDepthMetadata, SamplePloidyMetadata, SampleCoverageMetadata, \
    SampleMetadataCollection

_logger = logging.getLogger(__name__)


def write_sample_coverage_metadata(sample_metadata_collection: SampleMetadataCollection,
                                   sample_names: List[str],
                                   output_file: str):
    """Write coverage metadata for all samples in a given `SampleMetadataCollection` to a single .tsv file
    in the same order as `sample_names`.

    Args:
        sample_metadata_collection: an instance of `SampleMetadataCollection`
        sample_names: list of samples to process
        output_file: output .tsv file

    Raises:
        AssertionError: if some of the samples do not have `SampleCoverageMetadata` annotation

    Returns:
        None
    """
    assert len(sample_names) > 0
    assert sample_metadata_collection.all_samples_have_coverage_metadata(sample_names)
    contig_list = sample_metadata_collection.sample_coverage_metadata_dict[sample_names[0]].contig_list
    for sample_name in sample_names:
        assert sample_metadata_collection.sample_coverage_metadata_dict[sample_name].contig_list == contig_list
    parent_path = os.path.dirname(output_file)
    io_commons.assert_output_path_writable(parent_path)
    with open(output_file, 'w') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        header = [io_consts.sample_name_column_name] + [contig for contig in contig_list]
        writer.writerow(header)
        for sample_name in sample_names:
            sample_coverage_metadata = sample_metadata_collection.get_sample_coverage_metadata(sample_name)
            row = ([sample_name] + [repr(sample_coverage_metadata.n_j[j]) for j in range(len(contig_list))])
            writer.writerow(row)


def read_sample_coverage_metadata(sample_metadata_collection: SampleMetadataCollection,
                                  input_file: str,
                                  comment=io_consts.default_comment_char,
                                  delimiter=io_consts.default_delimiter_char) -> List[str]:
    """Reads sample coverage metadata from a .tsv file and adds them to `sample_metadata_collection`.

    Args:
        sample_metadata_collection: collection to which the coverage metadata is to be added
        input_file: input sample coverage metadata .tsv file
        comment: comment character
        delimiter: delimiter character

    Returns:
        list of samples in the same order as encountered in `input_file`
    """
    coverage_metadata_pd = pd.read_csv(input_file, delimiter=delimiter, comment=comment)
    found_columns_list = [str(column) for column in coverage_metadata_pd.columns.values]
    io_commons.assert_mandatory_columns({io_consts.sample_name_column_name}, set(found_columns_list), input_file)
    contig_list = found_columns_list.copy()
    contig_list.remove(io_consts.sample_name_column_name)
    num_contigs = len(contig_list)
    sample_names = []
    for tup in zip(coverage_metadata_pd[io_consts.sample_name_column_name],
                   *(coverage_metadata_pd[contig] for contig in contig_list)):
        sample_name = str(tup[0])
        n_j = np.asarray([int(tup[k + 1]) for k in range(num_contigs)], dtype=types.big_uint)
        sample_metadata_collection.add_sample_coverage_metadata(SampleCoverageMetadata(
            sample_name, n_j, contig_list))
        sample_names.append(sample_name)

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
