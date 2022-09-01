import argparse
import json
import re
from collections import defaultdict


def generate_avro_dict(avro_prefix, gcs_listing):
    """
    Generate a dictionary of the Avro arguments for the `hl.import_gvs` invocation. The keys will match the formal
    parameters of the `hl.import_gvs` function and should include:

    * vets (list of lists, one outer list per GVS superpartition of 4000 samples max)
    * refs (list of lists, one outer list per GVS superpartition of 4000 samples max)
    * sample_mapping_data (list)
    * site_filtering_data (list)
    * vqsr_filtering_data (list)
    * vqsr_tranche_data (list)
    """

    avro_file_arguments = defaultdict(list)
    superpartitioned_keys = {'vets', 'refs'}
    slashed_avro_prefix = avro_prefix + "/" if not avro_prefix[-1] == '/' else avro_prefix

    for full_path in gcs_listing.readlines():
        full_path = full_path.strip()
        if not full_path.endswith(".avro"):
            continue
        relative_path = full_path[len(slashed_avro_prefix):]
        parts = relative_path.split('/')
        key = parts[0]

        if key in superpartitioned_keys:
            # Get the zero based index from the `vet_001` component
            index = int(parts[1].split('_')[-1]) - 1
            if len(avro_file_arguments[key]) == index:
                avro_file_arguments[key].append([])
            avro_file_arguments[key][index].append(full_path)
        else:
            avro_file_arguments[key].append(full_path)

    return dict(avro_file_arguments)


def generate_avro_args(avro_dict):
    """
    Stringifies an Avro arguments dictionary such that it can be interpolated into a `hl.gvs_import` invocation to
    encompass all required Avro parameters.
    """
    avro_args = []
    for k, v in avro_dict.items():
        s = f'{k}={json.dumps(v, indent=4)}'
        avro_args.append(s)
    return ',\n'.join(avro_args)


def generate_gvs_import_script(avro_prefix, gcs_listing, vds_output_path, temp_dir):
    """
    Generate a Python script that when executed will:

    * import GVS Avro files into a Hail VDS
    * export from the Hail VDS to VCF

    """
    avro_dict = generate_avro_dict(avro_prefix, gcs_listing)
    avro_args = generate_avro_args(avro_dict)
    indented_avro_args = re.sub('\n', '\n    ', avro_args)
    hail_script = f"""
# * Introduction
#
# The following instructions can be used from the terminal of a Terra notebook to import GVS QuickStart Avro files
# and generate a VDS.
#
# * Hail installation
#
# Copy the appropriate Hail wheel locally and install:
#
# gsutil -m cp 'gs://gvs-internal-scratch/hail-wheels/<date>-<short git hash>/hail-<Hail version>-py3-none-any.whl' .
# pip install --force-reinstall hail-<Hail version>-py3-none-any.whl
#
# * Additional non-cluster configuration:
#
# If running with a Spark cluster, setup should be complete. If running with a non-cluster environment,
# the following lines are required in the terminal:
# 
# export PYSPARK_SUBMIT_ARGS='--driver-memory 16g --executor-memory 16g pyspark-shell'
# gcloud auth application-default login
# curl -sSL https://broad.io/install-gcs-connector | python3
#
# * Running the GVS to Hail VDS conversion script:
# 
# `python hail_gvs_import.py`

import hail as hl

rg38 = hl.get_reference('GRCh38')
rg38.add_sequence('gs://hail-common/references/Homo_sapiens_assembly38.fasta.gz',
                  'gs://hail-common/references/Homo_sapiens_assembly38.fasta.fai')
                  
vds_output_path="{vds_output_path}"

hl.import_gvs(
    {indented_avro_args},
    final_path=vds_output_path,
    tmp_dir="{temp_dir}",
    reference_genome=rg38,
)

"""
    return hail_script


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Generate Hail gvs_import invocation script')

    parser.add_argument('--avro_prefix', type=str, help='Avro prefix', required=True)
    parser.add_argument('--avro_listing_file', type=str, help='File containing a recursive listing under `avro_prefix',
                        required=True)
    parser.add_argument('--vds_output_path', type=str, help='GCS location for VDS output', required=True)
    parser.add_argument('--gcs_temporary_path', type=str, help='GCS location under which to create temporary files',
                        required=True)

    args = parser.parse_args()

    with open(args.avro_listing_file, 'r') as listing:
        script = generate_gvs_import_script(args.avro_prefix,
                                            listing,
                                            args.vds_output_path,
                                            args.gcs_temporary_path)
        print(script)
