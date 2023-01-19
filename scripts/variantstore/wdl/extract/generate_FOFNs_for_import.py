import os
import argparse

from terra_notebook_utils import table
from terra_notebook_utils import workspace

import time

# just keep these around as constants for now
vcf_files_name = "vcf_files.txt"
vcf_index_files_name = "vcf_index_files.txt"
sample_names_file_name = "sample_names.txt"
error_file_name = "errors.txt"


# make a default, but allow a user to overwrite it just in case someone wants to tweak this while tuning
attempts_between_pauses=500

def generate_FOFNs_from_data_table(data_table_name, sample_id_column_name, vcf_files_column_name, vcf_index_files_column_name):

    vcf_files = open(vcf_files_name, "w")
    vcf_index_files = open(vcf_index_files_name, "w")
    sample_names_file = open(sample_names_file_name, "w")
    error_file = open(error_file_name, "w")

    # Replace this with the actual data in the tables later
    for row in table.list_rows(data_table_name):
        try:
            current_sample_name = row.attributes[sample_id_column_name]
            current_vcf_file = row.attributes[vcf_files_column_name]
            current_vcf_index_file = row.attributes[vcf_index_files_column_name]

            sample_names_file.write(f'{current_sample_name}\n')
            vcf_files.write(f'{current_vcf_file}\n')
            vcf_index_files.write(f'{current_vcf_index_file}\n')
        except KeyError:
            error_file.write(f'Row "{row.name}" skipped: missing columns\n')

        count += 1
        if count >= attempts_between_pauses:
            print("sleeping between requests...")
            time.sleep(1)
            count = 0

    vcf_files.close()
    vcf_index_files.close()
    sample_names_file.close()
    error_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Turn the data tables into 3 FOFNs, one for sample names, vcf paths, and vcf index paths')

    parser.add_argument('--data_table_name', type=str,
                        help='The name of the table in your workspace that holds your sample data', default='sample')
    parser.add_argument('--sample_id_column_name', type=str, help='name of the column that holds the external sample names',
                            required=True)
    parser.add_argument('--vcf_files_column_name', type=str, help='name of the column that holds the paths to the vcf files',
                        required=True)
    parser.add_argument('--vcf_index_files_column_name', type=str, help='name of the column that holds the paths to the vcf index files', required=True)

    parser.add_argument('--attempts_between_pauses', type=int,
                        help='The number of rows in the db that are processed before we pause', default=500)

    args = parser.parse_args()

    # allow this to be overridden, but default it to 500
    if "attempts_between_pauses" in args:
        attempts_between_pauses = args.attempts_between_pauses


    generate_FOFNs_from_data_table(args.data_table_name,
                              args.sample_id_column_name,
                              args.vcf_files_column_name,
                              args.vcf_index_files_column_name)
