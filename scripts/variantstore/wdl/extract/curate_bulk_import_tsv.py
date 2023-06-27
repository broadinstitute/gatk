# -*- coding: utf-8 -*-
import argparse
import numpy as np
from typing import Union


def curate_bulk_import_data(sample_map_to_be_loaded_file_name: str,
                            bulk_import_input_file_name: str,
                            bulk_import_output_file_name: Union[str, None] = None) -> Union[None, list[list[str]]]:

    with open(sample_map_to_be_loaded_file_name, 'r') as sample_map_to_be_loaded_file, \
            open(bulk_import_input_file_name, 'r') as bulk_import_input_file:

        # The sample map has a header so skiprows=1.
        sample_map = np.loadtxt(sample_map_to_be_loaded_file, dtype=str, delimiter=",", ndmin=2, skiprows=1)
        # The sample map format is CSV with fields sample_id, sample_name. We want only the sample_name:
        samples_to_load = [arr[1] for arr in sample_map]

        bulk_import_data = np.loadtxt(bulk_import_input_file, delimiter="\t", dtype=str, ndmin=2).tolist()

        # Delete back to front to not throw off indices:
        for i in reversed(range(len(bulk_import_data))):
            # The 0th element in the bulk import data is the sample name.
            sample = bulk_import_data[i][0]
            if sample not in samples_to_load:
                #  If this isn't in the array of samples to load, delete the whole bulk import row.
                del(bulk_import_data[i])

        if bulk_import_output_file_name:
            # If an output file is specified write the curated bulk import data to that.
            with open(bulk_import_output_file_name, 'w') as output:
                np.savetxt(output, bulk_import_data, fmt='%s')
        else:
            # Otherwise return the curated bulk import data as a list[list[str]].
            return bulk_import_data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Curate bulk import TSV to remove already-loaded samples.')

    parser.add_argument('--samples-to-load', type=str, help='Path to a two column CSV *with a header* containing '
                                                            'sample_id, sample_name to load',
                        required=True)
    parser.add_argument('--bulk-import-input-tsv', type=str,
                        help='Path to input uncurated bulk import TSV: sample names, VCF paths, VCF index paths',
                        required=True)
    parser.add_argument('--bulk-import-output-tsv', type=str,
                        help='Path to output curated bulk import TSV: sample names, VCF paths, VCF index paths',
                        required=True)
    args = parser.parse_args()

    curate_bulk_import_data(args.samples_to_load,
                            args.bulk_import_input_tsv,
                            args.bulk_import_output_tsv)
