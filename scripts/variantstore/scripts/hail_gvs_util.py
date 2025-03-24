from google.cloud import storage
import re

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