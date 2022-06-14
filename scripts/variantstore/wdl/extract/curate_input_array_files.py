# -*- coding: utf-8 -*-
import numpy as np
from contextlib import contextmanager
import argparse

SAMPLE_MAP_TO_BE_LOADED_FILE_SUFFIX = "samples_to_be_loaded_map_file"
SAMPLE_NAME_FILE_SUFFIX = "sample_name_list_file"
VCF_FILE_SUFFIX = "vcf_list_file"
VCF_INDEX_FILE_SUFFIX = "vcf_index_list_file"

@contextmanager
def handle_file_error(file_name):
    try:
        yield
    except:
        print(f"ERROR: required file named '{file_name}' does not exist.")


def curate_input_arrays(sample_map_to_be_loaded_file_name,
                        sample_name_list_file_name,
                        vcf_list_file_name,
                        vcf_index_list_file_name,
                        output_files):
    sample_map_to_be_loaded_array = vcf_array = vcf_indexes_array = sample_names_array = []
    with handle_file_error(sample_map_to_be_loaded_file_name):
        sample_map_to_be_loaded_array = np.loadtxt(sample_map_to_be_loaded_file_name, dtype=str, delimiter=",")
    with handle_file_error(vcf_list_file_name):
        vcf_array = np.loadtxt(vcf_list_file_name, dtype=str)
    with handle_file_error(vcf_index_list_file_name):
        vcf_indexes_array = np.loadtxt(vcf_index_list_file_name, dtype=str)
    with handle_file_error(sample_name_list_file_name):
        sample_names_array = np.loadtxt(sample_name_list_file_name, dtype=str)
    rows_to_delete = []

    # use input_sample_names_array to figure out which index "rows" to delete
    for i in range(len(sample_names_array)):
        if sample_names_array[i] not in sample_map_to_be_loaded_array:
            rows_to_delete.append(i)

    # re-create input arrays using array of "rows" to delete
    vcf_array = [vcf_array[i] for i in range(len(vcf_array)) if i not in rows_to_delete]
    vcf_indexes_array = [vcf_indexes_array[i] for i in range(len(vcf_indexes_array)) if
                         i not in rows_to_delete]
    sample_names_array = [sample_names_array[i] for i in range(len(sample_names_array)) if
                          i not in rows_to_delete]

    if output_files:
        print(f"Creating 'output_{SAMPLE_NAME_FILE_SUFFIX}', 'output_{VCF_FILE_SUFFIX}' and 'output_{VCF_INDEX_FILE_SUFFIX}'.")
        np.savetxt(f"output_{SAMPLE_NAME_FILE_SUFFIX}", sample_names_array, fmt='%s')
        np.savetxt(f"output_{VCF_FILE_SUFFIX}", vcf_array, fmt='%s')
        np.savetxt(f"output_{VCF_INDEX_FILE_SUFFIX}", vcf_indexes_array, fmt='%s')
    else:
        d = dict();
        d['sample_names_array'] = sample_names_array
        d['vcf_array'] = vcf_array
        d['vcf_indexes_array'] = vcf_indexes_array
        return d

if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Curate GvsImportGenomes arrays to remove duplicate samples')

    parser.add_argument('--sample_map_to_be_loaded_file_name',type=str, help='name of sample_map file', required=False, default=f"input_{SAMPLE_MAP_TO_BE_LOADED_FILE_SUFFIX}")
    parser.add_argument('--sample_name_list_file_name',type=str, help='name of sample name list file', required=False, default=f"input_{SAMPLE_NAME_FILE_SUFFIX}")
    parser.add_argument('--vcf_list_file_name',type=str, help='name of VCF list file', required=False, default=f"input_{VCF_FILE_SUFFIX}")
    parser.add_argument('--vcf_index_list_file_name',type=str, help='name of VCF index list file', required=False, default=f"input_{VCF_INDEX_FILE_SUFFIX}")
    parser.add_argument('--output_files',type=bool, help='true (default): outputs are files; false: outputs are arrays', required=False, default=True)
    args = parser.parse_args()

    curate_input_arrays(args.sample_map_to_be_loaded_file_name,
                        args.sample_name_list_file_name,
                        args.vcf_list_file_name,
                        args.vcf_index_list_file_name,
                        args.output_files)
