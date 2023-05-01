import json
import requests

from terra_notebook_utils import table
from terra_notebook_utils import terra_auth
from terra_notebook_utils import workspace
import re

# make a default, but allow a user to overwrite it


def get_workspace_name(workspace_id):
    token=terra_auth.get_terra_access_token()
    # grab the workspace information from rawls
    rawls = 'https://rawls.dsde-prod.broadinstitute.org/api/workspaces/id/{}?fields=workspace.namespace,workspace.googleProject'.format(workspace_id)
    head = {'Authorization': 'token {}'.format(token)}
    response = requests.get(rawls, headers=head)
    # then extract the googleProject info
    print(response)
    proj_id=response.workspace.googleProject
    workspace_name=response.workspace.namespace
    return workspace_name


def get_sample_sets(workspace_namespace, workspace_name):
    response = requests.get('https://rawls.dsde-prod.broadinstitute.org/api/workspaces/{workspace_namespace}/{workspace_name}/entities?useCache=true')
    sample_sets = response.sample_set
    return sample_sets


def get_column_names(workspace_id, workspace_name):


def get_column_values(workspace_id, workspace_name):
    # We need to identify 3 things
    # 1. Sample id field
    # 2. vcf column name
    # 3. vcf index column name

    # We only sample a certain number of rows for each column
    numSamples = 50
    columnSamples = {}

    table_name = "sample" ## TODO is this an assumption that we are making? If the name of the table is sample, then wont the id be sample_id unless override?
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
        numSamples -= 1
        if numSamples == 0:
            break


    # time to start gathering some potential rows.
    # how many samples did we actually take?
    numSampledRows = 50 - numSamples
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



    print(f"ends_in_vcf: {ends_in_vcf}")
    print(f"ends_in_vcf_index: {ends_in_vcf_index}")
    print(f"contains_reblocked: {contains_reblocked}")

    print(f"path_ends_in_vcf_gz: {path_ends_in_vcf_gz}")
    print(f"path_ends_in_vcf_gz_tbi: {path_ends_in_vcf_gz_tbi}")



    # super simple heuristic: Is there a single entry that ends in vcf
    #for col in ends_in_vcf:


    # and has an analogue


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Get workspace information')

    parser.add_argument('--workspace_id', type=str,
                        help='The ID of your workspace that holds your sample data',
                        required=True)

    #parser.add_argument('--attempts_between_pauses', type=int,
     #                   help='The number of rows in the db that are processed before we pause', default=500)

    args = parser.parse_args()

    # allow this to be overridden, but default it to 500
    if "attempts_between_pauses" in args:
        attempts_between_pauses = args.attempts_between_pauses


    workspace_name = get_workspace_name(args.workspace_id)
    # get_column_names(args.workspace_id, workspace_name)
