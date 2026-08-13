from __future__ import annotations

import os
import re

from google.cloud import storage

try:
    import hail as hl
except ModuleNotFoundError:
    hl = None
gcs_pattern = re.compile("^gs://(?P<bucket_name>[^/]+)/(?P<object_prefix>.*)$")

def gcs_generate_avro_args(bucket, blob_prefix, key):
    """
    Generate a list of the Avro arguments for the `import_gvs` invocation for the specified key. The datatype should
    match these parameters:

    * vets (list of lists, one outer list per GVS superpartition of 4000 samples max)
    * refs (list of lists, one outer list per GVS superpartition of 4000 samples max)
    * sample_mapping (list)
    * site_filtering_data (list)
    * vets_filtering_data (list)
    """

    keyed_prefix = f"{blob_prefix}/{key}/"

    def superpartitioned_handler(blob_name):
        relative_path = blob_name[len(keyed_prefix):]
        parts = relative_path.split('/')

        index = int(parts[0].split('_')[-1]) - 1
        if index not in seen_indexes:
            seen_indexes.add(index)
            ret.append([])
        ret[-1].append(f'gs://{bucket.name}/{blob_name}')

    def regular_handler(blob_name):
        ret.append(f'gs://{bucket.name}/{blob_name}')

    superpartitioned_keys = {'vets', 'refs'}
    entry_handler = superpartitioned_handler if key in superpartitioned_keys else regular_handler

    seen_indexes = set()
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


def load_samples_to_remove(samples_to_remove_path: str) -> hl.Table:
    """Load and validate the samples-to-remove file, returning a Hail Table keyed by ``s``.

    The file must be a CSV with a ``research_id`` header column. A :class:`ValueError` is raised
    if the file contains duplicate ``research_id`` values.

    Parameters
    ----------
    samples_to_remove_path : str
        Path to a CSV file with a ``research_id`` header, one research ID per line.

    Returns
    -------
    hl.Table
        A Hail Table keyed by ``s`` (the Hail sample-column key), containing one row per
        unique research ID to remove.
    """
    samples_to_remove_table = hl.import_table(samples_to_remove_path, delimiter=',')

    # Fail fast if the removal file contains duplicate research_ids.
    research_ids = samples_to_remove_table.select('research_id')
    n_to_remove = research_ids.count()
    n_distinct = research_ids.key_by('research_id').distinct().count()
    if n_to_remove != n_distinct:
        raise ValueError(
            f"samples_to_remove file contains {n_to_remove - n_distinct} duplicate research_id(s). "
            "Please deduplicate before proceeding."
        )

    return samples_to_remove_table.key_by(s=samples_to_remove_table.research_id)


def filter_samples_and_remove_monomorphic_rows(
        vds: hl.vds.VariantDataset,
        samples_to_remove_table: hl.Table,
) -> hl.vds.VariantDataset:
    """Remove samples from a VDS and drop variant-data rows that become monomorphic reference.

    Steps performed:

    1. Filter out the specified samples and remove dead alleles (``remove_dead_alleles=True``).
    2. Drop variant-data rows that no longer carry any alternate alleles (monomorphic reference
       rows left behind when only the removed samples carried an alt allele at those sites).

    Parameters
    ----------
    vds : hl.vds.VariantDataset
        The input VDS from which samples will be removed.
    samples_to_remove_table : hl.Table
        A Hail Table keyed by ``s`` listing the samples to remove, as returned by
        :func:`load_samples_to_remove`.

    Returns
    -------
    hl.vds.VariantDataset
        A new VDS with the specified samples removed and monomorphic reference rows dropped.
    """
    filtered_vds = hl.vds.filter_samples(
        vds,
        samples_to_remove_table,
        keep=False,
        remove_dead_alleles=True,
    )

    # Drop rows that no longer have any non-reference calls after sample removal.
    filtered_vds = hl.vds.VariantDataset(
        filtered_vds.reference_data,
        filtered_vds.variant_data.filter_rows(
            hl.agg.any(filtered_vds.variant_data.LGT.is_non_ref())
        ),
    )

    return filtered_vds


def remove_trailing_slashes(*paths):
    return [p if not p.endswith('/') else p[:-1] for p in paths]


def determine_arguments_function(args):
    """
    Determine whether this function should look for Avro files in GCS or locally.

    :return: "filesystem"-appropriate avro arguments function.
    """
    avro_path, temp_path, vds_path = remove_trailing_slashes(args.avro_path, args.temp_path, args.vds_path)

    is_gcs = [gcs_pattern.match(p) for p in [avro_path, temp_path, vds_path]]

    if all(is_gcs):
        avro_bucket_name, avro_object_prefix = gcs_pattern.match(avro_path).groups()
        avro_bucket = storage.Client().get_bucket(avro_bucket_name)

        def arguments(key):
            return gcs_generate_avro_args(avro_bucket, avro_object_prefix, key)

        return arguments, True

    elif not any(is_gcs):
        if not args.references_path:
            raise ValueError(f"--references-path must be specified with local files")
        if gcs_pattern.match(args.references_path):
            raise ValueError(f"--references-path must refer to a local path")

        def arguments(key):
            return local_generate_avro_args(avro_path, key)

        return arguments, False

    else:
        raise ValueError("Arguments appear to be some unsavory mix of GCS and local paths, all or nothing please.")
