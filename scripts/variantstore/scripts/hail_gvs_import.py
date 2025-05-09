"""

Convenience script to wrap `import_gvs.py` for creating a VDS from a collection of Avro files exported from GVS.

"""

import argparse
import os

from hail_gvs_util import *


def create_vds(argsfn, vds_path, references_path, temp_path, intermediate_resume_point, skip_scoring):
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
            vets_filtering_data=argsfn('vets_filtering_data'),
            ploidy_data=argsfn('ploidy_data'),
            reference_genome=rg38,
            final_path=vds_path,
            tmp_dir=f'{temp_path}/hail_tmp_create_vds',
            # truth_sensitivity_snp_threshold: 'float' = 0.997,
            # truth_sensitivity_indel_threshold: 'float' = 0.990,
            # partitions_per_sample=0.35, # check with Hail about how to tune this for your large callset
            # intermediate_resume_point=0, # if your first run fails, and you want to use the intermediate files that already exist, check in with Hail to find out what stage to resume on
            # skip_final_merge=false, # if you want to create your VDS in two steps (because of mem issues) this can be skipped until the final run
            intermediate_resume_point=intermediate_resume_point,
            skip_scoring=skip_scoring
        )
    finally:
        local_hail_log_path = os.path.realpath(Env.hc()._log)
        fs = RouterFS()
        fs.copy(
            local_hail_log_path,
            f'{vds_path}.log'
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Create base Hail VDS from exported GVS Avro files.')
    parser.add_argument('--avro-path', type=str, help='Path at which exported GVS Avro files are found',
                        required=True)
    parser.add_argument('--vds-path', type=str, help='Path to which the VDS should be written',
                        required=True)
    parser.add_argument('--temp-path', type=str, help='Path to temporary directory',
                        required=True)
    parser.add_argument('--references-path', type=str, help='Path to references, only required for local files',
                        required=False)
    parser.add_argument('--intermediate-resume-point', type=int, required=False, default=0,
                        help='Intermediate VDS index at which to resume')
    parser.add_argument('--skip-scoring', type=bool, required=False, default=False, help='Whether to score the VDS or not')

    args = parser.parse_args()
    avro_path, temp_path, vds_path = remove_trailing_slashes(args.avro_path, args.temp_path, args.vds_path)

    arguments_fn, is_gcs = determine_arguments_function(args)

    references_path = 'gs://hail-common/references' if is_gcs else args.references_path
    create_vds(arguments_fn, vds_path, references_path, temp_path, args.intermediate_resume_point, args.skip_scoring)

