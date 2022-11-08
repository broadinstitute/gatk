"""
The following instructions can be used from the terminal of a Terra notebook to import GVS QuickStart Avro files
and generate a VDS.

* Hail installation

Copy the appropriate Hail wheel locally and install:

```
gsutil -m cp 'gs://gvs-internal-scratch/hail-wheels/<date>-<short git hash>/hail-<Hail version>-py3-none-any.whl' .
pip install --force-reinstall hail-<Hail version>-py3-none-any.whl
```
"""


from google.cloud import storage

import argparse
import os
import re

gcs_re = re.compile("^gs://(?P<bucket_name>[^/]+)/(?P<object_prefix>.*)$")


def import_gvs(argsfn, vds_path, references_path, temp_path):
    import hail as hl
    hl.init(tmp_dir=f'{temp_path}/hail_tmp_general')

    rg38 = hl.get_reference('GRCh38')
    rg38.add_sequence(f'{references_path}/Homo_sapiens_assembly38.fasta.gz',
                      f'{references_path}/Homo_sapiens_assembly38.fasta.fai')


    hl.import_gvs(
        vets=argsfn('vets'),
        refs=argsfn('refs'),
        sample_mapping=argsfn('sample_mapping'),
        site_filtering_data=argsfn('site_filtering_data'),
        vqsr_filtering_data=argsfn('vqsr_filtering_data'),
        vqsr_tranche_data=argsfn('vqsr_tranche_data'),
        reference_genome=rg38,
        final_path=vds_path,
        tmp_dir=f'{temp_path}/hail_tmp_import_gvs',
    )


def gcs_generate_avro_args(bucket, blob_prefix, key):
    """
    Generate a list of the Avro arguments for the `hl.import_gvs` invocation for the specified key. The datatype should
    match the parameters of the `hl.import_gvs` function:

    * vets (list of lists, one outer list per GVS superpartition of 4000 samples max)
    * refs (list of lists, one outer list per GVS superpartition of 4000 samples max)
    * sample_mapping (list)
    * site_filtering_data (list)
    * vqsr_filtering_data (list)
    * vqsr_tranche_data (list)
    """

    keyed_prefix = f"{blob_prefix}/{key}/"

    def superpartitioned_handler(blob_name):
        relative_path = blob_name[len(keyed_prefix):]
        parts = relative_path.split('/')

        index = int(parts[0].split('_')[-1]) - 1
        if len(ret) == index:
            ret.append([])
        ret[index].append(f'gs://{bucket.name}/{blob_name}')

    def regular_handler(blob_name):
        ret.append(f'gs://{bucket.name}/{blob_name}')

    superpartitioned_keys = {'vets', 'refs'}
    entry_handler = superpartitioned_handler if key in superpartitioned_keys else regular_handler

    ret = []

    # `list_blobs` paginates under the covers, explicit pagination not required regardless of the number of Avro files.
    # https://stackoverflow.com/a/43646557
    count = 0
    log_interval = 1000
    for blob in bucket.list_blobs(prefix=keyed_prefix):
        count = count + 1
        if count % log_interval == 0:
            print(f"Processed {count} {key} blobs...")

        if not blob.name.endswith(".avro"):
            continue
        entry_handler(blob.name)

    return ret


def local_generate_avro_args(avro_prefix, key):
    def superpartitioned_handler():
        parts = root.split('/')

        index = int(parts[-1].split('_')[-1]) - 1
        if len(ret) == index:
            ret.append([])
        ret[index].append(f'{root}/{file}')

    def regular_handler():
        ret.append(f'{root}/{file}')

    superpartitioned_keys = {'vets', 'refs'}
    entry_handler = superpartitioned_handler if key in superpartitioned_keys else regular_handler

    ret = []

    for root, dir, files in os.walk(f'{avro_prefix}/{key}'):
        for file in files:
            if file.endswith('avro'):
                entry_handler()
    return ret


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Create base Hail VDS from exported GVS Avro files.')
    parser.add_argument('--avro-path', type=str, help='Path at which exported GVS Avro files are found',
                        default="@AVRO_PATH")
    parser.add_argument('--vds-path', type=str, help='Path to which the VDS should be written', default="@VDS_PATH@")
    parser.add_argument('--temp-path', type=str, help='Path to temporary directory', default="@TEMP_DIR@")
    parser.add_argument('--references-path', type=str, help='Path to references, only required for local files')

    args = parser.parse_args()

    # Remove trailing slashes if present.
    avro_path, temp_path, vds_path = [p if not p.endswith('/') else p[:-1] for p in
                                      [args.avro_path, args.temp_path, args.vds_path]]

    is_gcs = [gcs_re.match(p) for p in [avro_path, temp_path, vds_path]]
    is_not_gcs = [not g for g in is_gcs]

    if all(is_gcs):
        avro_bucket_name, avro_object_prefix = gcs_re.match(avro_path).groups()
        avro_bucket = storage.Client().get_bucket(avro_bucket_name)

        def args(key):
            return gcs_generate_avro_args(avro_bucket, avro_object_prefix, key)

        import_gvs(args, vds_path, 'gs://hail-common/references', temp_path)

    elif all(is_not_gcs):
        references_path = args.references_path
        if not references_path:
            raise ValueError(f"--references-path must be specified with local files")
        if gcs_re.match(references_path):
            raise ValueError(f"--references-path must refer to a local path")

        def args(key):
            return local_generate_avro_args(avro_path, key)

        import_gvs(args, vds_path, references_path, temp_path)
    else:
        raise ValueError("Arguments appear to be some unsavory mix of GCS and local paths, all or nothing please.")
