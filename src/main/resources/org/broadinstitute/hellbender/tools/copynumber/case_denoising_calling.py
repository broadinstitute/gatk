import os
import sys

# set theano flags
user_theano_flags = os.environ.get("THEANO_FLAGS")
default_theano_flags = "device=cpu,floatX=float64,optimizer=fast_run,compute_test_value=ignore," + \
                       "openmp=true,blas.ldflags=-lmkl_rt,openmp_elemwise_minsize=10"
theano_flags = default_theano_flags + ("" if user_theano_flags is None else "," + user_theano_flags)
os.environ["THEANO_FLAGS"] = theano_flags

import logging
import argparse
import gcnvkernel
import shutil
import json
from typing import Dict, Any

logger = logging.getLogger("case_denoising_calling")

parser = argparse.ArgumentParser(description="gCNV case calling tool based on a previously trained model",
                                 formatter_class=gcnvkernel.cli_commons.GCNVHelpFormatter)

# logging args
gcnvkernel.cli_commons.add_logging_args_to_argparse(parser)

# add tool-specific args
group = parser.add_argument_group(title="Required arguments")

group.add_argument("--input_model_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Path to denoising model parameters")

group.add_argument("--read_count_tsv_files",
                   type=str,
                   required=True,
                   nargs='+',  # one or more
                   default=argparse.SUPPRESS,
                   help="List of read count files in the cohort (in .tsv format; must include sample name header)")

group.add_argument("--ploidy_calls_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="The path to the results of ploidy determination tool")

group.add_argument("--output_calls_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Output path to write CNV calls")

group.add_argument("--output_opt_path",
                   type=str,
                   required=False,
                   default=argparse.SUPPRESS,
                   help="(advanced) Output path to write the latest optimizer state")

group.add_argument("--output_tracking_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Output path to write tracked parameters, ELBO, etc.")

group.add_argument("--input_calls_path",
                   type=str,
                   required=False,
                   default=argparse.SUPPRESS,
                   help="Path to previously obtained calls to take as starting point")

group.add_argument("--input_opt_path",
                   type=str,
                   required=False,
                   default=argparse.SUPPRESS,
                   help="(advanced) Path to saved optimizer state to take as the starting point")

# add denoising config args
# Note: we are hiding parameters that are either set by the model or are irrelevant to the case calling task
gcnvkernel.DenoisingModelConfig.expose_args(
    parser,
    hide={
        "--max_bias_factors",
        "--psi_t_scale",
        "--log_mean_bias_std",
        "--init_ard_rel_unexplained_variance",
        "--enable_bias_factors",
        "--enable_explicit_gc_bias_modeling",
        "--disable_bias_factors_in_active_class",
        "--num_gc_bins",
        "--gc_curve_sd",
    })

# add calling config args
# Note: we are hiding parameters that are either set by the model or are irrelevant to the case calling task
gcnvkernel.CopyNumberCallingConfig.expose_args(
    parser,
    hide={
        '--p_active',
        '--class_coherence_length'
    })

# override some inference parameters
gcnvkernel.HybridInferenceParameters.expose_args(parser)


def update_args_dict_from_saved_model(input_model_path: str,
                                      _args_dict: Dict[str, Any]):
    logging.info("Loading denoising model configuration from the provided model...")
    with open(os.path.join(input_model_path, "denoising_config.json"), 'r') as fp:
        loaded_denoising_config_dict = json.load(fp)

    # boolean flags
    _args_dict['enable_bias_factors'] = \
        loaded_denoising_config_dict['enable_bias_factors']
    _args_dict['enable_explicit_gc_bias_modeling'] = \
        loaded_denoising_config_dict['enable_explicit_gc_bias_modeling']
    _args_dict['disable_bias_factors_in_active_class'] = \
        loaded_denoising_config_dict['disable_bias_factors_in_active_class']

    # bias factor related
    _args_dict['max_bias_factors'] = \
        loaded_denoising_config_dict['max_bias_factors']

    # gc-related
    _args_dict['num_gc_bins'] = \
        loaded_denoising_config_dict['num_gc_bins']
    _args_dict['gc_curve_sd'] = \
        loaded_denoising_config_dict['gc_curve_sd']

    logging.info("- bias factors enabled: "
                 + repr(_args_dict['enable_bias_factors']))
    logging.info("- explicit GC bias modeling enabled: "
                 + repr(_args_dict['enable_explicit_gc_bias_modeling']))
    logging.info("- bias factors in active classes disabled: "
                 + repr(_args_dict['disable_bias_factors_in_active_class']))

    if _args_dict['enable_bias_factors']:
        logging.info("- maximum number of bias factors: "
                     + repr(_args_dict['max_bias_factors']))

    if _args_dict['enable_explicit_gc_bias_modeling']:
        logging.info("- number of GC curve knobs: "
                     + repr(_args_dict['num_gc_bins']))
        logging.info("- GC curve prior standard deviation: "
                     + repr(_args_dict['gc_curve_sd']))


if __name__ == "__main__":

    # parse arguments
    args = parser.parse_args()
    gcnvkernel.cli_commons.set_logging_config_from_args(args)

    logger.info("THEANO_FLAGS environment variable has been set to: {theano_flags}".format(theano_flags=theano_flags))

    # check gcnvkernel version in the input model path
    gcnvkernel.io_commons.check_gcnvkernel_version_from_path(args.input_model_path)

    # copy the intervals to the calls path
    # (we do this early to avoid inadvertent cleanup of temporary files)
    gcnvkernel.io_commons.assert_output_path_writable(args.output_calls_path)
    shutil.copy(os.path.join(args.input_model_path, gcnvkernel.io_consts.default_interval_list_filename),
                os.path.join(args.output_calls_path, gcnvkernel.io_consts.default_interval_list_filename))

    # load modeling interval list from the model
    logging.info("Loading modeling interval list from the provided model...")
    modeling_interval_list = gcnvkernel.io_intervals_and_counts.load_interval_list_tsv_file(
        os.path.join(args.input_model_path, gcnvkernel.io_consts.default_interval_list_filename))
    contigs_set = {target.contig for target in modeling_interval_list}
    logging.info("The model contains {0} intervals and {1} contig(s)".format(
        len(modeling_interval_list), len(contigs_set)))

    # load sample names, truncated counts, and interval list from the sample read counts table
    logging.info("Loading {0} read counts file(s)...".format(len(args.read_count_tsv_files)))
    sample_names, n_st = gcnvkernel.io_intervals_and_counts.load_counts_in_the_modeling_zone(
        args.read_count_tsv_files, modeling_interval_list)

    # load read depth and ploidy metadata
    sample_metadata_collection: gcnvkernel.SampleMetadataCollection = gcnvkernel.SampleMetadataCollection()
    gcnvkernel.io_metadata.update_sample_metadata_collection_from_ploidy_determination_calls(
        sample_metadata_collection, args.ploidy_calls_path)

    # setup the inference task
    args_dict = args.__dict__

    # read model configuration and update args dict
    update_args_dict_from_saved_model(args.input_model_path, args_dict)

    # instantiate config classes
    denoising_config = gcnvkernel.DenoisingModelConfig.from_args_dict(args_dict)
    calling_config = gcnvkernel.CopyNumberCallingConfig.from_args_dict(args_dict)
    inference_params = gcnvkernel.HybridInferenceParameters.from_args_dict(args_dict)

    # instantiate and initialize the workspace
    shared_workspace = gcnvkernel.DenoisingCallingWorkspace(
        denoising_config, calling_config, modeling_interval_list,
        n_st, sample_names, sample_metadata_collection)

    initial_params_supplier = gcnvkernel.DefaultDenoisingModelInitializer(
        denoising_config, calling_config, shared_workspace)

    task = gcnvkernel.CaseDenoisingCallingTask(
        denoising_config, calling_config, inference_params,
        shared_workspace, initial_params_supplier, args.input_model_path)

    if hasattr(args, 'input_calls_path'):
        logger.info("A call path was provided to use as starting point...")
        gcnvkernel.io_denoising_calling.SampleDenoisingAndCallingPosteriorsReader(
            shared_workspace, task.continuous_model, task.continuous_model_approx,
            args.input_calls_path)()

    if hasattr(args, 'input_opt_path'):
        logger.info("A saved optimizer state was provided to use as starting point...")
        task.fancy_opt.load(args.input_opt_path)

    try:
        # go!
        task.engage()
        task.disengage()
    except gcnvkernel.ConvergenceError as err:
        logger.info(err.message)
        # if inference diverged, pass an exit code to the Java side indicating that restart is needed
        sys.exit(gcnvkernel.io_consts.diverged_inference_exit_code)


    # save calls
    gcnvkernel.io_denoising_calling.SampleDenoisingAndCallingPosteriorsWriter(
        denoising_config, calling_config, shared_workspace, task.continuous_model, task.continuous_model_approx,
        args.output_calls_path)()

    # save optimizer state
    if hasattr(args, 'output_opt_path'):
        task.fancy_opt.save(args.output_opt_path)

    # save ELBO history
    if hasattr(args, 'output_tracking_path'):
        gcnvkernel.io_commons.assert_output_path_writable(args.output_tracking_path)

        elbo_hist_file = os.path.join(args.output_tracking_path, "elbo_history.tsv")
        task.save_elbo_history(elbo_hist_file)
