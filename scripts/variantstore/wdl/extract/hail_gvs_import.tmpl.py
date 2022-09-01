"""
The following instructions can be used from the terminal of a Terra notebook to import GVS QuickStart Avro files
and generate a VDS.

* Hail installation

Copy the appropriate Hail wheel locally and install:

```
gsutil -m cp 'gs://gvs-internal-scratch/hail-wheels/<date>-<short git hash>/hail-<Hail version>-py3-none-any.whl' .
pip install --force-reinstall hail-<Hail version>-py3-none-any.whl
```

* Additional non-cluster configuration:

If running with a Spark cluster, setup should be complete. If running with a non-cluster environment,
the following lines are required in the terminal:

```
export PYSPARK_SUBMIT_ARGS='--driver-memory 16g --executor-memory 16g pyspark-shell'
gcloud auth application-default login
curl -sSL https://broad.io/install-gcs-connector | python3
```

* Running the GVS to Hail VDS conversion script:

`python hail_gvs_import.tmpl.py`
"""
import os

from google.cloud import storage

import hail as hl
import re
import uuid


def generate_avro_args(bucket, object_prefix, key):
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

    object_keyed_prefix=f"{object_prefix}/{key}/"
    full_keyed_prefix = f"gs://{bucket.name}/{object_keyed_prefix}"

    def superpartitioned_handler(full_path):
        relative_path = full_path[len(full_keyed_prefix):]
        parts = relative_path.split('/')

        index = int(parts[0].split('_')[-1]) - 1
        if len(ret) == index:
            ret.append([])
        ret[index].append(full_path)

    def regular_handler(full_path):
        ret.append(full_path)

    superpartitioned_keys = {'vets', 'refs'}
    entry_handler = superpartitioned_handler if key in superpartitioned_keys else regular_handler
    bucket_str = bucket.name

    ret = []

    # `list_blobs` paginates under the covers: https://stackoverflow.com/a/43646557
    count = 0
    log_interval = 1000
    for blob in bucket.list_blobs(prefix=object_keyed_prefix):
        count = count + 1
        if count % log_interval == 0:
            print(f"Processed {count} {key} blobs...")

        if not blob.name.endswith(".avro"):
            continue
        full_path = f'gs://{bucket_str}/{blob.name}'
        entry_handler(full_path)

    return ret


def import_gvs(bucket, object_prefix, vds_output_path, tmp_dir):

    rg38 = hl.get_reference('GRCh38')
    rg38.add_sequence('gs://hail-common/references/Homo_sapiens_assembly38.fasta.gz',
                      'gs://hail-common/references/Homo_sapiens_assembly38.fasta.fai')

    def args(key):
        return generate_avro_args(bucket, object_prefix, key)

    hl.import_gvs(
        vets=args('vets'),
        refs=args('refs'),
        sample_mapping=args('sample_mapping'),
        site_filtering_data=args('site_filtering_data'),
        vqsr_filtering_data=args('vqsr_filtering_data'),
        vqsr_tranche_data=args('vqsr_tranche_data'),
        reference_genome=rg38,
        final_path=vds_output_path,
        tmp_dir=tmp_dir
    )


if __name__ == '__main__':

    avro_prefix = "$AVRO_PREFIX$"
    output_prefix = "$OUTPUT_PREFIX$"
    # Remove a trailing slash if present.
    output_prefix = output_prefix if not output_prefix.endswith("/") else output_prefix[:-1]

    gcs_re = re.compile("^gs://(?P<bucket_name>[^/]+)/(?P<object_prefix>.*)$")
    match = gcs_re.match(avro_prefix)
    if not match:
        raise ValueError(f"Avro prefix '{avro_prefix}' does not look like a GCS path")
    bucket_name, object_prefix = match.groups()

    client = storage.Client()
    bucket = client.get_bucket(bucket_name)

    # A little pseudorandomness to write without clobbering.
    rand = uuid.uuid4().hex[:8]
    temp_dir = f"{output_prefix}/temp_hail_gvs_import_{rand}"
    vds_output_path = f"{output_prefix}/final-{rand}.vds"
    print(f"Using temporary dir '{temp_dir}', final VDS to be written to '{vds_output_path}'.")

    import_gvs(bucket, object_prefix, vds_output_path, temp_dir)
