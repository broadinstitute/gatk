# -*- coding: utf-8 -*-
import argparse
import csv
import numpy as np

SAMPLE_MAP_TO_BE_LOADED_FILE_SUFFIX = "samples_to_be_loaded_map_file"
SAMPLE_NAME_FILE_SUFFIX = "sample_name_list_file"
VCF_FILE_SUFFIX = "vcf_list_file"
VCF_INDEX_FILE_SUFFIX = "vcf_index_list_file"


def read_single_column_file(file_name):
    ret = []
    with open(file_name, 'r') as fd:
        reader = csv.reader(fd)
        for row in reader:
            ret.append(row[0])
    return ret

def curate_input_arrays(sample_map_to_be_loaded_file_name,
                        sample_name_list_file_name,
                        vcf_list_file_name,
                        vcf_index_list_file_name,
                        output_files):

    sample_map = set()

    with open(sample_map_to_be_loaded_file_name, 'r') as fd:
        reader = csv.reader(fd)
        for row in reader:
            sample_map.add(row[1])
    vcfs = read_single_column_file(vcf_list_file_name)
    vcf_indexes = read_single_column_file(vcf_index_list_file_name)
    sample_names = read_single_column_file(sample_name_list_file_name)

    # Work back to front to not mess up indexes if deleting.
    for i in reversed(range(len(sample_names))):
        if sample_names[i] not in sample_map:
            del vcfs[i]
            del vcf_indexes[i]
            del sample_names[i]

    if output_files:
        print(f"Creating 'output_{SAMPLE_NAME_FILE_SUFFIX}', 'output_{VCF_FILE_SUFFIX}' and 'output_{VCF_INDEX_FILE_SUFFIX}'.")
        np.savetxt(f"output_{SAMPLE_NAME_FILE_SUFFIX}", sample_names, fmt='%s')
        np.savetxt(f"output_{VCF_FILE_SUFFIX}", vcfs, fmt='%s')
        np.savetxt(f"output_{VCF_INDEX_FILE_SUFFIX}", vcf_indexes, fmt='%s')
    else:
        d = dict();
        d['sample_names_array'] = sample_names
        d['vcf_array'] = vcfs
        d['vcf_indexes_array'] = vcf_indexes
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