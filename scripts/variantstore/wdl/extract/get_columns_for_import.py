import requests
import argparse

from terra_notebook_utils import table
import re


maxNumSamples = 50

def get_column_data(workspace_id):
    # We need to identify 3 things
    # 1. Sample id field
    # 2. vcf column name
    # 3. vcf index column name

    # We only sample a certain number of columns
    numSamples = 0
    columnSamples = {}

    table_name = "sample"
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

def get_column_values(columnSamples, numSamples):
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

    # this is for returning debug info to the user, showing the multiple columns that matched?
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
        also_contains_reblocked = matching_vcf_columns.intersect(contains_reblocked)
        if len(also_contains_reblocked) == 1:
            # They likely just had a raw AND reblocked columns in the data. Pick the reblocked one
            found_vcf_column = True
            final_vcf_column = also_contains_reblocked.pop()
            final_vcf_index_column = f"{final_vcf_column}_index"

    content_matching_vcfs = set()
    # We STILL weren't able to uniquely determine the correct columns for the vcf and index files going by column names?
    if not found_vcf_column:
        # Check the contents of the columns: the duck algorithm.  If its contents LOOK like like vcfs and indexes, go from there
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
        if len(matching_vcf_columns) > 1:
            print(f"Multiple columns that look like the right naming convention: {matching_vcf_columns}")
        if len(content_matching_vcfs) > 1:
            print(f"Multiple columns that look like have contents that look like zipped vcfs: {content_matching_vcfs}")
    else:
        print("\nDerived values:");
        print(f"final_vcf_column: {final_vcf_column}")
        print(f"final_vcf_index_column: {final_vcf_index_column}")

    # debug vomit
    print("\n\nDebug info:")
    print(f"ends_in_vcf: {ends_in_vcf}")
    print(f"ends_in_vcf_index: {ends_in_vcf_index}")
    print(f"contains_reblocked: {contains_reblocked}")

    print(f"path_ends_in_vcf_gz: {path_ends_in_vcf_gz}")
    print(f"path_ends_in_vcf_gz_tbi: {path_ends_in_vcf_gz_tbi}")

    column_names = {
        'vcf_column' : final_vcf_column,
        'vcf_column_index': final_vcf_index_column
    }

    return column_names


def write_column_names(column_names, vcf_output, vcf_index_output):
    final_vcf_column = column_names['vcf_column']
    final_vcf_index_column = column_names['vcf_column_index']
    with open(vcf_output, "w") as vcf_output, open(vcf_index_output, "w") as vcf_index_output:
        vcf_output.write(f'{final_vcf_column}\n')
        vcf_index_output.write(f'{final_vcf_index_column}\n')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Get workspace information')

    parser.add_argument('--workspace_id', type=str,
                        help='The ID of your workspace that holds your sample data',
                        required=True)

    parser.add_argument('--workspace_name', type=str,
                        help='The name of your workspace that holds your sample data',
                        required=False)

    parser.add_argument('--vcf_output', type=str,
                        help='The location to write the suggested vcf col name',
                        required=False)

    parser.add_argument('--vcf_index_output', type=str,
                        help='The location to write the suggested index col name',
                        required=False)

    parser.add_argument('--attempts_between_pauses', type=int,
                        help='The number of rows in the db that are processed before we pause', default=500)

    args = parser.parse_args()

    # allow this to be overridden, but default it to 500
    if "attempts_between_pauses" in args:
        attempts_between_pauses = args.attempts_between_pauses


    (columnSamples, numSamples) = get_column_data(args.workspace_id)
    column_names = get_column_values(columnSamples, numSamples)
    write_column_names(column_names, args.vcf_output, args.vcf_index_output)
