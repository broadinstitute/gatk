import os

# set theano flags
os.environ["THEANO_FLAGS"] = "device=cpu,floatX=float64,optimizer=fast_run,compute_test_value=ignore,openmp=true"

import argparse
import gcnvkernel
import shutil

parser = argparse.ArgumentParser(description="gCNV contig ploidy and read depth determination tool",
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

group.add_argument("--sample_coverage_metadata",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Coverage metadata of all samples (in .tsv format)")

group.add_argument("--contig_ploidy_prior_table",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Contig ploidy prior probabilities (in .tsv format)")

group.add_argument("--output_model_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Output path to write the ploidy model for future case-sample ploidy determination use")

group.add_argument("--output_calls_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Output path to write posteriors")

# optional arguments
gcnvkernel.PloidyModelConfig.expose_args(parser)

# override some inference parameters
gcnvkernel.HybridInferenceParameters.expose_args(
    parser,
    override_default={
        "--learning_rate": 0.1,
        "--adamax_beta2": 0.999,
        "--log_emission_samples_per_round": 1000,
        "--log_emission_sampling_rounds": 50,
        "--log_emission_sampling_median_rel_error": 1e-3,
        "--max_advi_iter_first_epoch": 5000,
        "--max_advi_iter_subsequent_epochs": 1000,
        "--convergence_snr_averaging_window": 5000,
        "--convergence_snr_countdown_window": 100,
        "--num_thermal_advi_iters": 10000,
        "--max_calling_iters": 1,
        "--caller_update_convergence_threshold": 1e-3
    })

if __name__ == "__main__":

    # parse arguments
    args = parser.parse_args()
    gcnvkernel.cli_commons.set_logging_config_from_args(args)

    # read contig ploidy prior map from file
    contig_ploidy_prior_map = gcnvkernel.io_ploidy.get_contig_ploidy_prior_map_from_tsv_file(
        args.contig_ploidy_prior_table)

    # load interval list
    interval_list = gcnvkernel.io_intervals_and_counts.load_interval_list_tsv_file(args.interval_list)

    # load sample coverage metadata
    sample_metadata_collection: gcnvkernel.SampleMetadataCollection = gcnvkernel.SampleMetadataCollection()
    sample_names = gcnvkernel.io_metadata.read_sample_coverage_metadata(
        sample_metadata_collection, args.sample_coverage_metadata)

    # generate interval list metadata
    intervals_metadata: gcnvkernel.IntervalListMetadata = gcnvkernel.IntervalListMetadata(interval_list)

    # inject ploidy prior map to the dictionary of parsed args
    args_dict = args.__dict__
    args_dict['contig_ploidy_prior_map'] = contig_ploidy_prior_map

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
    shutil.copy(args.contig_ploidy_prior_table,
                os.path.join(args.output_model_path, gcnvkernel.io_consts.default_contig_ploidy_prior_tsv_filename))
