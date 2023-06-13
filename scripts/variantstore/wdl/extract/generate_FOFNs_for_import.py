import argparse
import time

from terra_notebook_utils import table


def get_entities_in_set(data_table_name, sample_set_name):
    set_of_entities = set()
    sample_set_list = list(table.list_rows(f"{data_table_name}_set"))
    for sample_set in sample_set_list:
        if sample_set.name == sample_set_name: ## else it's some other sample set and we dont care about it
            for sample in sample_set.attributes[f"{data_table_name}s"]:
                set_of_entities.add(sample['entityName'])

    return set_of_entities


## Note that there are two different possible <entity>_id options:
## There is the id column that is used by the <entity>_set to differentiate what is and isn't in the set / inclusion list
## The user may want to use a _different_ column other than the <entity>_id column for sample level info. E.g. the "research_id" col
## The user passes this in as the sample_id_column_name, which is set to be the <entity>_id if it is missing


def generate_FOFNs_from_data_table_with_sample_set(
        data_table_name, sample_id_column_name,
        vcf_files_name, vcf_index_files_name,
        vcf_files_column_name, vcf_index_files_column_name, sample_names_file_name, error_file_name,
        set_of_entities):
    with open(vcf_files_name, "w") as vcf_files, open(vcf_index_files_name, "w") as vcf_index_files, open(sample_names_file_name, "w") as sample_names_file, open(error_file_name, "w") as error_file:
        count = 0
        processed_entities = 0
        # Cycle through each row / sample in the data table.
        # If there is an inclusion list (a sample_set), check to see if the sample is in it.
        # Then write what is needed for the FOFN.
        for row in table.list_rows(data_table_name):
            try:
                # The table data model is weird, as each row has a "name" property followed by a hash of attributes.  But it
                # is all displayed as though they are equivalent columns, and the name field is mapped to an imaginary
                # column with the name table_name_id.  Map "name" to that imaginary attribute here so none of the lookups
                # below fail if they reference that call.  This is safe, as this Row object is a read only copy of the actual
                # table content, and it is disposed
                row.attributes[f"{data_table_name}_id"] = row.name  ## NOTE: this is the <entity>_id value and will be used by the sample_set

                current_sample_id = row.attributes[f"{data_table_name}_id"] ## NOTE: this is the sample id based on the <entity>_id
                current_sample_name = row.attributes[sample_id_column_name] ## NOTE: this is the sample name based on the user defined sample id col name

                if set_of_entities:
                    if current_sample_id in set_of_entities:
                        current_vcf_file = row.attributes[vcf_files_column_name]
                        current_vcf_index_file = row.attributes[vcf_index_files_column_name]
                        sample_names_file.write(f'{current_sample_name}\n')
                        vcf_files.write(f'{current_vcf_file}\n')
                        vcf_index_files.write(f'{current_vcf_index_file}\n')
                else:
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
    parser.add_argument('--sample_id_column_name', type=str,
                        help='The name of the column that holds the external sample names')
    parser.add_argument('--vcf_files_column_name', type=str,
                        help='The name of the column that holds the paths to the vcf files',
                        required=True)
    parser.add_argument('--vcf_index_files_column_name', type=str,
                        help='The name of the column that holds the paths to the vcf index files',
                        required=True)
    parser.add_argument('--attempts_between_pauses', type=int,
                        help='The number of rows in the db that are processed before we pause', default=500)
    parser.add_argument('--sample_set_name', type=str,
                        help='The name of the entity / sample set to use', default=None)
    parser.add_argument('--sample_names_file_name', type=str,
                        help='The text file we save the sample column name to', default="sample_names.txt")
    parser.add_argument('--vcf_files_name', type=str,
                        help='The text file we save the vcf column name to', default="vcf_files.txt")
    parser.add_argument('--vcf_index_files_name', type=str,
                        help='The text file we save the vcf index column name to', default="vcf_index_files.txt")
    parser.add_argument('--error_file_name', type=str,
                        help='The text file we save error msgs to', default="errors.txt")

    args = parser.parse_args()


    # allow this to be overridden, but default it to 500
    if "attempts_between_pauses" in args:
        attempts_between_pauses = args.attempts_between_pauses

    # check for selected sample sets
    if "sample_set_name" in args:
        sample_set_name = args.sample_set_name

    ## Note: check to see if there is a user defined sample name column
    if "sample_id_column_name" in args:
        sample_id_column_name = args.sample_id_column_name
    else:
        sample_id_column_name = args.data_table_name + "_id" ## <-- NOTE: This makes the assumption that the default setting for the GVS sample name is <entity>_id if not specified



    if sample_set_name:
        set_of_entities = get_entities_in_set(
            args.data_table_name,
            sample_set_name)

    else:
        set_of_entities = None # set_of_entities is now nothing

    generate_FOFNs_from_data_table_with_sample_set(
        args.data_table_name,
        sample_id_column_name,
        args.vcf_files_name,
        args.vcf_index_files_name,
        args.vcf_files_column_name,
        args.vcf_index_files_column_name,
        args.sample_names_file_name,
        args.error_file_name,
        set_of_entities)
