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
    with open(vcf_files_name, "w") as vcf_files, open(vcf_index_files_name, "w") as vcf_index_files, open(sample_names_file_name, "w") as sample_names_file, open(error_file_name, "w") as error_file
        count = 0
        processed_entities = 0
        for row in table.list_rows(data_table_name):
            try:
                # The table data model is weird, as each row has a "name" property followed by a hash of attributes.  But it
                # is all displayed as though they are equivalent columns, and the name field is mapped to an imaginary
                # column with the name table_name_id.  Map "name" to that imaginary attribute here so none of the lookups
                # below fail if they reference that call.  This is safe, as this Row object is a read only copy of the actual
                # table content, and it is disposed
                row.attributes[f"{data_table_name}_id"] = row.name

                current_sample_name = row.attributes[sample_id_column_name]
                current_vcf_file = row.attributes[vcf_files_column_name]
                current_vcf_index_file = row.attributes[vcf_index_files_column_name]

                sample_names_file.write(f'{current_sample_name}\n')
                vcf_files.write(f'{current_vcf_file}\n')
                vcf_index_files.write(f'{current_vcf_index_file}\n')
            except KeyError:
                error_file.write(f'Row "{row.name}" skipped: missing columns\n')

            count += 1
            processed_entities += 1
            if count >= attempts_between_pauses:
                print(f"sleeping between requests ({processed_entities} processed)...")
                time.sleep(1)
                count = 0

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
