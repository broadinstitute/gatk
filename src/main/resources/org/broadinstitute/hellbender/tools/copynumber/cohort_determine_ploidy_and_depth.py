import os

# set theano flags
os.environ["THEANO_FLAGS"] = "device=cpu,floatX=float64,optimizer=fast_run,compute_test_value=ignore,openmp=true"

import argparse
import gcnvkernel
import shutil

parser = argparse.ArgumentParser(description="gCNV contig ploidy and read-depth determination tool",
                                 formatter_class=gcnvkernel.cli_commons.GCNVHelpFormatter)

# logging args
gcnvkernel.cli_commons.add_logging_args_to_argparse(parser)

# add tool-specific args
group = parser.add_argument_group(title="Required arguments")

group.add_argument("--interval_list",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Interval list of included genomic regions in the analysis (in .tsv format)")

group.add_argument("--contig_count_distribution_collection_files",
                   type=str,
                   required=True,
                   nargs='+',  # one or more
                   default=argparse.SUPPRESS,
                   help="List of per-contig count-distribution files for all samples (in .tsv format; must include sample name header)")

group.add_argument("--ploidy_state_priors_table",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Ploidy-state prior probabilities (in .tsv format)")

group.add_argument("--output_model_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Output path to write the ploidy model for future case-sample ploidy determination")

group.add_argument("--output_calls_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Output path to write posteriors")

# group.add_argument("--output_plots_path",
#                    type=str,
#                    required=True,
#                    default=argparse.SUPPRESS,
#                    help="Output path to write plots")

# optional arguments
gcnvkernel.PloidyModelConfig.expose_args(parser)

# override some inference parameters
gcnvkernel.HybridInferenceParameters.expose_args(parser)

if __name__ == "__main__":

    # parse arguments
    args = parser.parse_args()
    gcnvkernel.cli_commons.set_logging_config_from_args(args)

    # read ploidy-state prior map from file
    ploidy_state_priors_map = gcnvkernel.io_ploidy.get_ploidy_state_priors_map_from_tsv_file(
        args.ploidy_state_priors_table)

    # load interval list
    interval_list = gcnvkernel.io_intervals_and_counts.load_interval_list_tsv_file(args.interval_list)

    # load sample coverage metadata
    sample_metadata_collection: gcnvkernel.SampleMetadataCollection = gcnvkernel.SampleMetadataCollection()
    sample_names = gcnvkernel.io_metadata.read_sample_coverage_metadata(
        sample_metadata_collection, args.contig_count_distribution_collection_files)

    # generate interval list metadata
    intervals_metadata: gcnvkernel.IntervalListMetadata = gcnvkernel.IntervalListMetadata(interval_list)

    # inject ploidy-state priors map to the dictionary of parsed args
    args_dict = args.__dict__
    args_dict['ploidy_state_priors_map'] = ploidy_state_priors_map

    ploidy_config = gcnvkernel.PloidyModelConfig.from_args_dict(args_dict)
    ploidy_inference_params = gcnvkernel.HybridInferenceParameters.from_args_dict(args_dict)
    ploidy_workspace = gcnvkernel.PloidyWorkspace(ploidy_config, intervals_metadata, sample_names,
                                                  sample_metadata_collection)
    ploidy_task = gcnvkernel.CohortPloidyInferenceTask(ploidy_inference_params, ploidy_config, ploidy_workspace)

    # go!
    ploidy_task.engage()
    ploidy_task.disengage()

    # save model parameters
    gcnvkernel.io_ploidy.PloidyModelWriter(ploidy_config, ploidy_workspace,
                                           ploidy_task.continuous_model, ploidy_task.continuous_model_approx,
                                           args.output_model_path)()

    # sample sample-specific posteriors
    gcnvkernel.io_ploidy.SamplePloidyWriter(ploidy_config, ploidy_workspace,
                                            ploidy_task.continuous_model, ploidy_task.continuous_model_approx,
                                            args.output_calls_path)()

    # save a copy of interval list and ploidy priors as well
    shutil.copy(args.interval_list,
                os.path.join(args.output_model_path, gcnvkernel.io_consts.default_interval_list_filename))
    shutil.copy(args.ploidy_state_priors_table,
                os.path.join(args.output_model_path, gcnvkernel.io_consts.default_ploidy_state_prior_tsv_filename))
