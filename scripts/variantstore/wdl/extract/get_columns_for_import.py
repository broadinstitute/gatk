import argparse
import re
from terra_notebook_utils import table


# The goal of this code is to validate and determine the 5 values for:
# 1. The entity type (default: sample)
# 2. The entity type id column name. This is never anything but #1 + _id (default: sample_id)
# 3. The sample name column name (default: sample_id if it exists, otherwise <entity>_id)
# 4. The input GVCFs path column name. This has no explicit default, but uses heuristics above to be determined
# 5. The input GVCF index files path column name. This has no explicit default, but uses heuristics above to be determined

# The method get_entity_data returns the values for the first two
# The method get_column_data tells us how much data we want to survey
# The method get_column_values looks in the columns and gives us the best guess for the fourth and fifth values

### NOTE: FOR NOW WE JUST ALLOW "SAMPLE" TO BE THE DEFAULT ENTITY TYPE AND THIS CODE JUST VALIDATES THAT



maxNumSamples = 50


def get_entity_data(user_defined_entity, entity_set):
    tables_list=list(table.list_tables())
    ## When the user has supplied a defined entity
    if user_defined_entity:
        # validate that the user defined entity is present in the workspace
        if user_defined_entity in tables_list:
            entity = user_defined_entity
            # if there is a user specified <entity>_set, validate that it is of the expected entity type
            if entity_set:
                entity_set_list=[]
                [entity_set_list.append(row.name) for row in table.list_rows(entity+"_set")]
                if entity_set not in entity_set_list:
                    print("error in user defined set")

        else:
            print("error in user defined entity type")

    # In the case that there is only one or two possible entities, they can easily be derived
    elif len(tables_list) == 1:
        entity = tables_list[0]
    elif len(tables_list) == 2:
        if tables_list[0]+"_set" ==tables_list[1]:
            entity = tables_list[0]
        elif tables_list[1]+"_set" ==tables_list[0]:
            entity = tables_list[1]


    ## If the user has selected an entity_set name (but not defined an entity), grab all tables that end in "_set"
    # and validate that the named <entity>_set is present in at least one entity type, and if it's in only one, set that as the entity type
    elif entity_set:
        entity_table_that_ends_in_set=[]
        for entity_table in tables_list:
            if entity_table[-4:]=="_set":
                entity_table_that_ends_in_set.append(entity_table)

        #check that the named <entity>_set is present as at least one set
        presence=[]
        for set_table in entity_table_that_ends_in_set:
            for row in table.list_rows(set_table):
                if row.name == entity_set:
                    presence.append(set_table[0:-4])
            if len(presence) > 1:
                break
        if len(presence) == 1:
            entity = presence[0]

    elif len(tables_list) == 1:
        entity = tables_list[0]
    elif len(tables_list) == 2:
        if tables_list[0]+"_set" ==tables_list[1]:
            entity = tables_list[0]
        elif tables_list[1]+"_set" ==tables_list[0]:
            entity = tables_list[1]
    else:
        ## TODO: move this to the top or do it throughout the gates--not sure which is better. Currently it is set immediately by the default value in the argparse arg
        print("default set to entity type sample")
        entity = default_entity

    if entity not in tables_list:
        print("error in possible entities")

    return entity


def get_column_data(entity_type):
    # We need to identify 3 things
    # 1. entity id field (just <entity_type>_id)
    # 2. vcf column name
    # 3. vcf index column name

    # We only sample a certain number of columns
    numSamples = 0
    columnSamples = {}

    table_name = entity_type # e.g.: "sample"
    # list_rows is a generator and it uses paging behind the scenes to make
    # this call much more efficient
    for row in table.list_rows(table_name):
        # the id field is special, so handle it like it is...
        if f"{table_name}_id" in columnSamples:
            existingList = columnSamples[f"{table_name}_id"]
            existingList.append(row.name)
        else:
            newList = [row.name];
            columnSamples[f"{table_name}_id"] = newList

        # handle the more general attributes
        for columnName in row.attributes:
            if columnName in columnSamples:
                existingList = columnSamples[columnName]
                existingList.append(row.name)
            else:
                columnSamples[columnName] = [row.attributes[columnName]]

        # done iterating columns
        numSamples += 1
        if numSamples == maxNumSamples:
            break
    return (columnSamples, numSamples)

def get_column_values(columnSamples, numSamples, user_defined_vcf, user_defined_index):
    if user_defined_vcf:
        if user_defined_vcf not in columnSamples: # validate that the user defined vcf column is present in the table
            raise ValueError("vcf column not in table")
        if user_defined_index:
            if user_defined_index not in columnSamples: # validate that the user defined index column is present in the table
                raise ValueError("vcf index column not in table")
            else:
                return (user_defined_vcf, user_defined_index)
        elif user_defined_vcf + "_index" in columnSamples: # If no index col was provided, but there is one that matches a provided vcf column, we go with that
            return (user_defined_vcf, user_defined_vcf + "_index")
        else:
            raise ValueError(f"No index column specified for: {user_defined_vcf}")
    elif user_defined_index:
        raise ValueError(f"No GVCF column specified for: {user_defined_index}")

    # time to start gathering some potential rows.
    # how many samples did we actually take?
    numSampledRows = numSamples
    cutoffPoint = numSampledRows * 0.95
    print(f"Sampled {numSampledRows} rows total. Throwing away any under {cutoffPoint}")

    # match column names that end with "vcf"
    ends_in_vcf_pattern = "^.*vcf$"
    ends_in_vcf = set()

    # match column names that end with "vcf_index"
    ends_in_vcf_index_pattern = "^.*vcf_index$"
    ends_in_vcf_index = set()

    # match column names that contain the word "reblocked"
    contains_reblocked_pattern = ".*reblocked.*"
    contains_reblocked = set()

    # match path that end with ".vcf.gz"
    path_ends_in_vcf_gz_pattern = "^.*\.vcf\.gz$"
    path_ends_in_vcf_gz = set()

    # match path that end with ".vcf.gz.tbi"
    path_ends_in_vcf_gz_tbi_pattern = "^.*\.vcf\.gz\.tbi$"
    path_ends_in_vcf_gz_tbi = set()

    # start sorting the columns that we've seen into buckets to be compared later

    for key in columnSamples:
        print(f"Found key: {key} with {len(columnSamples[key])} entries")
        samplesData = columnSamples[key]
        if len(samplesData) < cutoffPoint:
            # ignoring this completely
            continue

        if type(samplesData[0]) != str:
            # ignoring non-strings completely.
            # We can grab the 0th index because the way that the values come in, there are no missing vals
            continue

        # ends in vcf?
        result = re.search(ends_in_vcf_pattern, key)
        if result:
            ends_in_vcf.add(key)

        # ends in vcf_index?
        result = re.search(ends_in_vcf_index_pattern, key)
        if result:
            ends_in_vcf_index.add(key)

        # contains the word reblocked?
        result = re.search(contains_reblocked_pattern, key)
        if result:
            contains_reblocked.add(key)

        # has a path that ends in vcf.gz
        result = re.search(path_ends_in_vcf_gz_pattern, samplesData[0])
        if result:
            path_ends_in_vcf_gz.add(key)

        # has a path that ends in vcf.gz.tbi
        result = re.search(path_ends_in_vcf_gz_tbi_pattern, samplesData[0])
        if result:
            path_ends_in_vcf_gz_tbi.add(key)


    final_vcf_column = ""
    final_vcf_index_column = ""

    found_vcf_column = False

    # this is for returning debug info to the user, showing the multiple columns that matched
    matching_vcf_columns = set()

    # super simple heuristic: Are there columns that end in vcf?
    for col in ends_in_vcf:
        # ...and have an analogue that ends in vcf_index?
        index_column = f"{col}_index"
        if index_column in ends_in_vcf_index:
            # woohoo!
            final_vcf_column = col
            final_vcf_index_column = index_column
            found_vcf_column = True

            matching_vcf_columns.add(col)

    # check the contents of the columns if we haven't found a single one yet
    if len(matching_vcf_columns) > 1:
        # reset instead of creating yet another boolean to track this
        found_vcf_column = False

        # We have multiple potential columns to go with.  Two ways to trim them down
        # 1: See if any specifically have "reblocked" in their names.  If that gets us down to 1, go with it
        also_contains_reblocked = matching_vcf_columns.intersection(contains_reblocked)
        if len(also_contains_reblocked) == 1:
            # They likely just had a raw AND reblocked columns in the data. Pick the reblocked one
            found_vcf_column = True
            final_vcf_column = also_contains_reblocked.pop()
            final_vcf_index_column = f"{final_vcf_column}_index"

    content_matching_vcfs = set()
    # We STILL weren't able to uniquely determine the correct columns for the vcf and index files going by column names?
    if not found_vcf_column:
        # Check the contents of the columns: the duck algorithm.  If its contents LOOK like vcfs and indexes, go from there
        for col in path_ends_in_vcf_gz:
            # ...and has an analogue that looks like an index file?
            index_column = f"{column}_index"
            # is this the correct logic?  Or do we just want to look for ANY singular column that has contents that
            # look liks a vcf and ANY other singular column with contents that look like an index file?  Stick with
            # enforcing a naming convention for now...
            if index_column in path_ends_in_vcf_gz_tbi:
                # woohoo!
                final_vcf_column = column
                final_vcf_index_column = index_column
                found_vcf_column = True

                content_matching_vcfs.add(col)

    if len(content_matching_vcfs) > 1:
        # This means we went through the second loop above and still weren't able to uniquely identify vcf and index columns
        # Try further heuristics or throw our hands up at this point
        found_vcf_column = False


    if not found_vcf_column:
        print("Unable to uniquely identify columns for the vcf and vcf index files")
        error_string = "No appropriate columns were found for the vcf and vcf index file"
        if len(matching_vcf_columns) > 1:
            error_string = f"Multiple columns that look like the right naming convention: {matching_vcf_columns}"
        if len(content_matching_vcfs) > 1:
            error_string = f"Multiple columns that look like have contents that look like zipped vcfs: {content_matching_vcfs}"

        # debug vomit
        print("\n\nDebug info:")
        print(f"ends_in_vcf: {ends_in_vcf}")
        print(f"ends_in_vcf_index: {ends_in_vcf_index}")
        print(f"contains_reblocked: {contains_reblocked}")
        print(f"path_ends_in_vcf_gz: {path_ends_in_vcf_gz}")
        print(f"path_ends_in_vcf_gz_tbi: {path_ends_in_vcf_gz_tbi}")

        raise ValueError(error_string)
    else:
        print("\nDerived values:");
        print(f"final_vcf_column: {final_vcf_column}")
        print(f"final_vcf_index_column: {final_vcf_index_column}")

    return (final_vcf_column, final_vcf_index_column)

def write_column_names(final_vcf_column, final_vcf_index_column, vcf_output, vcf_index_output):
    with open(vcf_output, "w") as vcf_output, open(vcf_index_output, "w") as vcf_index_output:
        vcf_output.write(f'{final_vcf_column}\n')
        vcf_index_output.write(f'{final_vcf_index_column}\n')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Get workspace information')

    parser.add_argument('--entity_type', type=str,
                        help='The name of the entity being ingested--we default to sample')

    parser.add_argument('--entity_set_name', type=str,
                        help='Optional set of entities--often called a sample_set. Having a sample_set is recommended for loading large amounts of data',
                        required=False)

    parser.add_argument('--user_defined_sample_id', type=str,
                        help='The column that the user would like to use as the sample name in GVS and in the final extract',
                        required=False)

    parser.add_argument('--vcf_output', type=str,
                        help='The location to write the suggested vcf column name',
                        required=False)

    parser.add_argument('--vcf_index_output', type=str,
                        help='The location to write the suggested index column name',
                        required=False)

    parser.add_argument('--user_defined_vcf', type=str,
                        help='The column that the user has specified as the vcf column name',
                        required=False)

    parser.add_argument('--user_defined_index', type=str,
                        help='The column that the user has specified as the index column name',
                        required=False)



    args = parser.parse_args()

    # allow this to be overridden, but default it to 500
    if "attempts_between_pauses" in args:
        attempts_between_pauses = args.attempts_between_pauses

    # if the user is using sample_sets, then that will potentially help us get the first (and thus also the second) or validate a user defined entity_type
    # 1. The entity type
    entity_type = ""
    if "entity_set_name" in args:
        entity_type = get_entity_data(args.entity_type, args.entity_set_name)
    # if there is no sample_set, we must rely on the user defining the entity type
    elif "entity_type" in args:
        entity_type = args.entity_type
    # or let it default to "sample"
    else:
        entity_type = "sample"

    # 2. The entity type id column name
    entity_id = entity_type + "_id"

    # 3. The sample name column name
    sample_id = ""
    if "user_defined_entity_id" in args:
        sample_id = args.user_defined_entity_id
    else:
        sample_id = entity_id

    # validate the entity information
    (columnSamples, numSamples) = get_column_data(entity_type)

    # 4 and # 5: The input GVCFs path column name and the input GVCF index path column name
    (vcf_column, vcf_column_index) = get_column_values(columnSamples, numSamples, args.user_defined_vcf, args.user_defined_index)

    # Write the derived values to their corresponding locations -- we just pass the vcf paths for now, as at this point, everything else is just for validation rather than discovery
    write_column_names(vcf_column, vcf_column_index, args.vcf_output, args.vcf_index_output)
