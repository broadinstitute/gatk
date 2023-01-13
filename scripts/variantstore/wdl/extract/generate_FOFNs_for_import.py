import os
import argparse

from terra_notebook_utils import table
from terra_notebook_utils import workspace

from google.cloud import storage
from pathlib import Path

import utils

# just keep these around as constants for now
vcf_files_name = "vcf_files.txt"
vcf_index_files_name = "vcf_index_files.txt"
sample_names_file_name = "sample_names.txt"


def generate_FOFNs_from_data_table(data_table_name, vcf_files_column_name, vcf_index_files_column_name):

    vcf_files = open(vcf_files_name, "w")
    vcf_index_files = open(vcf_index_files_name, "w")
    sample_names_file = open(sample_names_file_name, "w")

    # Replace this with the actual data in the tables later
    for row in table.list_rows(data_table_name):
        sample_names_file.write(f"{row.name}\n")
        vcf_files.write(f'{row.attributes[vcf_files_column_name]}\n')
        vcf_index_files.write(f'{row.attributes[vcf_index_files_column_name]}\n')

    vcf_files.close()
    vcf_index_files.close()
    sample_names_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Turn the data tables into 3 FOFNs, one for sample names, vcf paths, and vcf index paths')

    parser.add_argument('--data_table_name', type=str,
                        help='The name of the table in your workspace that holds your sample data', default='sample')
    parser.add_argument('--vcf_files_column_name', type=str, help='name of the column that holds the paths to the vcf files',
                        required=True)
    parser.add_argument('--vcf_index_files_column_name', type=str, help='name of the column that holds the paths to the vcf index files', required=True)

    args = parser.parse_args()

    generate_FOFNs_from_data_table(args.data_table_name,
                              args.vcf_files_column_name,
                              args.vcf_index_files_column_name)
