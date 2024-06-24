"""

Convenience script to wrap `import_gvs.py` for creating a VDS from a collection of Avro files exported from GVS.

"""


from google.cloud import storage

import argparse
import os
import re

gcs_re = re.compile("^gs://(?P<bucket_name>[^/]+)/(?P<object_prefix>.*)$")


def create_vds(argsfn, vds_path, references_path, temp_path, use_classic_vqsr, intermediate_resume_point):
    import hail as hl
    import import_gvs
    from hail.utils.java import Env
    from hailtop.fs.router_fs import RouterFS

    hl.init(tmp_dir=f'{temp_path}/hail_tmp_general')
    hl._set_flags(use_new_shuffle='1')

    rg38 = hl.get_reference('GRCh38')
    rg38.add_sequence(f'{references_path}/Homo_sapiens_assembly38.fasta.gz',
                      f'{references_path}/Homo_sapiens_assembly38.fasta.fai')

    # A full description of the `import_gvs` function written by Hail for this process can be found in `import_gvs.py`:
    # https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/wdl/extract/import_gvs.py
    # Commented out parameters are ones where we are comfortable with the default, but want to make them easily
    # accessible to users.
    try:
        import_gvs.import_gvs(
            vets=argsfn('vets'),
            refs=argsfn('refs'),
            sample_mapping=argsfn('sample_mapping'),
            site_filtering_data=argsfn('site_filtering_data'),
            vqsr_filtering_data=argsfn('vqsr_filtering_data'),
            vqsr_tranche_data=argsfn('vqsr_tranche_data'),
            reference_genome=rg38,
            final_path=vds_path,
            tmp_dir=f'{temp_path}/hail_tmp_create_vds',
            # truth_sensitivity_snp_threshold: 'float' = 0.997,
            # truth_sensitivity_indel_threshold: 'float' = 0.990,
            # partitions_per_sample=0.35, # check with Hail about how to tune this for your large callset
            # intermediate_resume_point=0, # if your first run fails, and you want to use the intermediate files that already exist, check in with Hail to find out what stage to resume on
            # skip_final_merge=false, # if you want to create your VDS in two steps (because of mem issues) this can be skipped until the final run
            use_classic_vqsr=use_classic_vqsr,
            intermediate_resume_point=intermediate_resume_point
        )
    finally:
        local_hail_log_path = os.path.realpath(Env.hc()._log)
        fs = RouterFS()
        fs.copy(
            local_hail_log_path,
            f'{vds_path}.log'
        )


def gcs_generate_avro_args(bucket, blob_prefix, key):
    """
    Generate a list of the Avro arguments for the `import_gvs` invocation for the specified key. The datatype should
    match these parameters:

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
                        default="@AVRO_PATH@",
                        required=True)
    parser.add_argument('--vds-path', type=str, help='Path to which the VDS should be written',
                        default="@VDS_PATH@",
                        required=True)
    parser.add_argument('--temp-path', type=str, help='Path to temporary directory',
                        default="@TEMP_DIR@",
                        required=True)
    parser.add_argument('--references-path', type=str, help='Path to references, only required for local files',
                        required=False)
    parser.add_argument("--use-classic-vqsr", action="store_true",
                        help="If set, expect that the input GVS Avro files were generated using VQSR Classic")
    parser.add_argument('--intermediate-resume-point', type=int, required=False, default=0,
                        help='Intermediate VDS index at which to resume')

    args = parser.parse_args()

    # Remove trailing slashes if present.
    avro_path, temp_path, vds_path = [p if not p.endswith('/') else p[:-1] for p in
                                      [args.avro_path, args.temp_path, args.vds_path]]
    use_classic_vqsr =  args.use_classic_vqsr
    is_gcs = [gcs_re.match(p) for p in [avro_path, temp_path, vds_path]]
    is_not_gcs = [not g for g in is_gcs]

    if all(is_gcs):
        avro_bucket_name, avro_object_prefix = gcs_re.match(avro_path).groups()
        avro_bucket = storage.Client().get_bucket(avro_bucket_name)

        def arguments(key):
            return gcs_generate_avro_args(avro_bucket, avro_object_prefix, key)

        create_vds(arguments, vds_path, 'gs://hail-common/references', temp_path, use_classic_vqsr,
                   args.intermediate_resume_point)

    elif all(is_not_gcs):
        references_path = args.references_path
        if not references_path:
            raise ValueError(f"--references-path must be specified with local files")
        if gcs_re.match(references_path):
            raise ValueError(f"--references-path must refer to a local path")

        def arguments(key):
            return local_generate_avro_args(avro_path, key)

        create_vds(arguments, vds_path, references_path, temp_path, use_classic_vqsr,
                   args.intermediate_resume_point)
    else:
        raise ValueError("Arguments appear to be some unsavory mix of GCS and local paths, all or nothing please.")
