import os
import argparse

from terra_notebook_utils import table
from terra_notebook_utils import workspace

import time



def get_entities_in_set(data_table_name, sample_set_name):
    set_of_entities = set()
    sample_set_list = list(table.list_rows(f"{data_table_name}_set"))
    for sample_set in sample_set_list:
        if sample_set.name == sample_set_name: ## else it's some other sample set and we dont care about it
            for sample in sample_set.attributes['samples']: ## TODO why is this samples in Terra? Is this some Terra hardcoding?!?!
                set_of_entities.add(sample['entityName'])

    return set_of_entities


## Note that there are two different possible <entity>_id options:
## There is the id col that is used by the <entity>_set to differntiate what is and isn't in the set / inclusion list
## The user may want to use a _different_ column other than the <entity>_id col for sample level info. E.g. the "research_id" col


def generate_FOFNs_from_data_table_with_sample_set(data_table_name, sample_id_column_name, vcf_files_column_name, vcf_index_files_column_name, set_of_entities):
    with open(args.vcf_files_name, "w") as vcf_files, open(args.vcf_index_files_name, "w") as vcf_index_files, open(args.sample_names_file_name, "w") as sample_names_file, open(args.error_file_name, "w") as error_file:
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
                if current_sample_name in set_of_entities:
                    current_vcf_file = row.attributes[vcf_files_column_name]
                    current_vcf_index_file = row.attributes[vcf_index_files_column_name]
                    sample_names_file.write(f'{current_sample_name}\n')
                    vcf_files.write(f'{current_vcf_file}\n')
                    vcf_index_files.write(f'{current_vcf_index_file}\n')

            except KeyError:
                error_file.write(f'Row "{row.name}" skipped: missing columns\n')

            count += 1
            processed_entities += 1 ## TODO is this just for the sample set we are looking at!??!!
            if count >= attempts_between_pauses:
                print(f"sleeping between requests ({processed_entities} processed)...")
                time.sleep(1)
                count = 0

def generate_FOFNs_from_data_table(data_table_name, sample_id_column_name, vcf_files_column_name, vcf_index_files_column_name):
    with open(vcf_files_name, "w") as vcf_files, open(vcf_index_files_name, "w") as vcf_index_files, open(sample_names_file_name, "w") as sample_names_file, open(error_file_name, "w") as error_file:
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

    parser.add_argument('--entity_set_name', type=str, help='The name of the entity / sample set to use', default=None)
    parser.add_argument('--user_defined_sample_id_col_name', type=str, help='The name of the sample set to use', default=None)

    parser.add_argument('--sample_names_file_name', type=str, help='The name of the sample set to use', default="sample_names.txt")
    parser.add_argument('--vcf_files_name', type=str, help='The name of the sample set to use', default="vcf_files.txt")
    parser.add_argument('--vcf_index_files_name', type=str, help='The name of the sample set to use', default="vcf_index_files.txt")
    parser.add_argument('--error_file_name', type=str, help='The name of the sample set to use', default="errors.txt")

    args = parser.parse_args()


    # allow this to be overridden, but default it to 500
    if "attempts_between_pauses" in args:
        attempts_between_pauses = args.attempts_between_pauses

    # check for selected sample sets
    if "sample_set_name" in args:
        sample_set_name = args.entity_set_name


    if sample_set_name:
        ## Note: check to see if there is a suer defined id column
        if "user_defined_sample_id_col_name" in args:
            sample_id_column_name = args.user_defined_sample_id_col_name
        else:
            sample_id_column_name = args.sample_id_column_name

        set_of_entities = get_entities_in_set(
            args.data_table_name,
            sample_set_name)
        generate_FOFNs_from_data_table_with_sample_set(
            args.data_table_name,
            sample_id_column_name,
            args.vcf_files_column_name,
            args.vcf_index_files_column_name,
            set_of_entities)
    else:
        generate_FOFNs_from_data_table(
            args.data_table_name,
            args.sample_id_column_name,
            args.vcf_files_column_name,
            args.vcf_index_files_column_name)
