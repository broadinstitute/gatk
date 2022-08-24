import argparse
import json
import re
from collections import defaultdict


def generate_avro_dict(avro_prefix, gcs_listing):
    avro_file_arguments = defaultdict(list)
    super_partitioned_keys = {'vets', 'refs'}
    slashed_avro_prefix = avro_prefix + "/" if not avro_prefix[-1] == '/' else avro_prefix

    for full_path in gcs_listing.readlines():
        full_path = full_path.strip()
        if not full_path.endswith(".avro"):
            continue
        relative_path = full_path[len(slashed_avro_prefix):]
        parts = relative_path.split('/')
        key = parts[0]

        if key in super_partitioned_keys:
            # Get the zero based index from the `vet_001` component
            index = int(parts[1].split('_')[-1]) - 1
            if len(avro_file_arguments[key]) == index:
                avro_file_arguments[key].append([])
            avro_file_arguments[key][index].append(full_path)
        else:
            avro_file_arguments[key].append(full_path)

    return dict(avro_file_arguments)


def generate_avro_args(avro_dict):
    args = []
    for k, v in avro_dict.items():
        s = f'{k}={json.dumps(v, indent=4)}'
        args.append(s)
    return ',\n'.join(args)


def generate_gvs_import_script(avro_prefix, listing, vds_output_path, vcf_output_path, temp_dir):
    avro_dict = generate_avro_dict(avro_prefix, listing)
    avro_args = generate_avro_args(avro_dict)
    indented_avro_args = re.sub('\n', '\n    ', avro_args)
    script = f"""
# The following instructions can be used from the terminal of a Terra notebook to import GVS QuickStart Avro files
# and generate a VDS.
#
# Copy the appropriate Hail wheel locally first:
#
# gsutil -m cp 'gs://gvs-internal-scratch/hail-wheels/2022-08-18-01f7b77ebbcc/hail-0.2.97-py3-none-any.whl' .
#
# If running locally (non-Spark cluster) set this in the environment before launching Python:
# export PYSPARK_SUBMIT_ARGS='--driver-memory 16g --executor-memory 16g pyspark-shell'
#
# Hail wants Java 8, Java 11+ will not do. Make sure you have a Java 8 in your path with `java -version`.
#
# pip install hail-0.2.97-py3-none-any.whl
# gcloud auth application-default login
# curl -sSL https://broad.io/install-gcs-connector | python3
#

import hail as hl

rg38 = hl.get_reference('GRCh38')
rg38.add_sequence('gs://hail-common/references/Homo_sapiens_assembly38.fasta.gz',
                  'gs://hail-common/references/Homo_sapiens_assembly38.fasta.fai')

hl.import_gvs(
    {indented_avro_args},
    final_path="{vds_output_path}",
    tmp_dir="{temp_dir}",
    reference_genome=rg38,
)

vds = hl.vds.read_vds('{vds_output_path}')

mt = hl.vds.to_dense_mt(vds)
fail_case = 'FAIL'
mt = mt.annotate_entries(FT=hl.if_else(mt.FT, 'PASS', fail_case))
hl.export_vcf(mt, '{vcf_output_path}')
"""
    return script


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Generate Hail gvs_import invocation script')

    parser.add_argument('--avro_prefix', type=str, help='Avro prefix', required=True)
    parser.add_argument('--avro_listing_file', type=str, help='File containing a recursive listing under `avro_prefix',
                        required=True)
    parser.add_argument('--vds_output_path', type=str, help='GCS location for VDS output', required=True)
    parser.add_argument('--vcf_output_path', type=str, help='GCS location for VCF output generated from VDS',
                        required=True)
    parser.add_argument('--gcs_temporary_path', type=str, help='GCS location under which to create temporary files',
                        required=True)

    args = parser.parse_args()

    with open(args.avro_listing_file, 'r') as listing:
        script = generate_gvs_import_script(args.avro_prefix,
                                            listing,
                                            args.vds_output_path,
                                            args.vcf_output_path,
                                            args.gcs_temporary_path)
        print(script)
