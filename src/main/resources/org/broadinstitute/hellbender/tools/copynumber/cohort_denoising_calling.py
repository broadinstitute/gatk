import os

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

logger = logging.getLogger("cohort_denoising_calling")

parser = argparse.ArgumentParser(description="gCNV cohort denoising and calling tool",
                                 formatter_class=gcnvkernel.cli_commons.GCNVHelpFormatter)

# logging args
gcnvkernel.cli_commons.add_logging_args_to_argparse(parser)

# add tool-specific args
group = parser.add_argument_group(title="Required arguments")

group.add_argument("--modeling_interval_list",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Full interval list, possibly including extra annotations (in .tsv format)")

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

group.add_argument("--output_model_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Output path to write model parameters")

group.add_argument("--output_calls_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Output path to write CNV calls")

group.add_argument("--output_tracking_path",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Output path to write tracked parameters, ELBO, etc.")

group.add_argument("--output_opt_path",
                   type=str,
                   required=False,
                   default=argparse.SUPPRESS,
                   help="(advanced) Output path to write the latest optimizer state")

group.add_argument("--input_model_path",
                   type=str,
                   required=False,
                   default=argparse.SUPPRESS,
                   help="(advanced) Path to previously obtained model parameters to take as starting point")

group.add_argument("--input_calls_path",
                   type=str,
                   required=False,
                   default=argparse.SUPPRESS,
                   help="(advanced) Path to previously obtained calls to take as starting point")

group.add_argument("--input_opt_path",
                   type=str,
                   required=False,
                   default=argparse.SUPPRESS,
                   help="(advanced) Path to saved optimizer state to take as the starting point")

# add denoising config args
gcnvkernel.DenoisingModelConfig.expose_args(parser)

# add calling config args
gcnvkernel.CopyNumberCallingConfig.expose_args(parser)

# override some inference parameters
gcnvkernel.HybridInferenceParameters.expose_args(
    parser)

if __name__ == "__main__":

    # parse arguments
    args = parser.parse_args()
    gcnvkernel.cli_commons.set_logging_config_from_args(args)

    logger.info("THEANO_FLAGS environment variable has been set to: {theano_flags}".format(theano_flags=theano_flags))

    # copy the intervals to the model and calls paths
    # (we do this early to avoid inadvertent cleanup of temporary files)
    gcnvkernel.io_commons.assert_output_path_writable(args.output_model_path)
    shutil.copy(args.modeling_interval_list,
                os.path.join(args.output_model_path, gcnvkernel.io_consts.default_interval_list_filename))
    gcnvkernel.io_commons.assert_output_path_writable(args.output_calls_path)
    shutil.copy(args.modeling_interval_list,
                os.path.join(args.output_calls_path, gcnvkernel.io_consts.default_interval_list_filename))

    # load modeling interval list
    modeling_interval_list = gcnvkernel.io_intervals_and_counts.load_interval_list_tsv_file(args.modeling_interval_list)

    # load sample names, truncated counts, and interval list from the sample read counts table
    logger.info("Loading {0} read counts file(s)...".format(len(args.read_count_tsv_files)))
    sample_names, n_st = gcnvkernel.io_intervals_and_counts.load_counts_in_the_modeling_zone(
        args.read_count_tsv_files, modeling_interval_list)

    # load read depth and ploidy metadata
    sample_metadata_collection: gcnvkernel.SampleMetadataCollection = gcnvkernel.SampleMetadataCollection()
    gcnvkernel.io_metadata.update_sample_metadata_collection_from_ploidy_determination_calls(
        sample_metadata_collection, args.ploidy_calls_path)

    # setup the inference task
    args_dict = args.__dict__

    # instantiate config classes for the main task
    main_calling_config = gcnvkernel.CopyNumberCallingConfig.from_args_dict(args_dict)
    main_denoising_config = gcnvkernel.DenoisingModelConfig.from_args_dict(args_dict)
    main_inference_params = gcnvkernel.HybridInferenceParameters.from_args_dict(args_dict)

    # instantiate config classes for the warm-up task
    warm_up_denoising_config = gcnvkernel.DenoisingModelConfig.from_args_dict(args_dict)
    warm_up_inference_params = gcnvkernel.HybridInferenceParameters.from_args_dict(args_dict)

    # update the contigs for the two-stage task
    warm_up_inference_params.min_training_epochs = 1
    warm_up_inference_params.max_training_epochs = 1
    warm_up_denoising_config.q_c_expectation_mode = 'marginalize'
    main_inference_params.max_advi_iter_first_epoch = 0  # the main task should start from sampling and calling

    # instantiate and initialize the workspace
    shared_workspace = gcnvkernel.DenoisingCallingWorkspace(
        main_denoising_config, main_calling_config, modeling_interval_list,
        n_st, sample_names, sample_metadata_collection)

    shared_workspace.initialize_copy_number_class_inference_vars()

    initial_params_supplier = gcnvkernel.DefaultDenoisingModelInitializer(
        main_denoising_config, main_calling_config, shared_workspace)

    # warm-up task
    warm_up_task = gcnvkernel.CohortDenoisingAndCallingWarmUpTask(
        warm_up_denoising_config, warm_up_inference_params,
        shared_workspace, initial_params_supplier)

    if hasattr(args, 'input_model_path'):
        logger.info("A model path was provided to use as starting point...")
        gcnvkernel.io_denoising_calling.DenoisingModelReader(
            warm_up_denoising_config, main_calling_config, shared_workspace, warm_up_task.continuous_model,
            warm_up_task.continuous_model_approx, args.input_model_path)()

    if hasattr(args, 'input_calls_path'):
        logger.info("A call path was provided to use as starting point...")
        gcnvkernel.io_denoising_calling.SampleDenoisingAndCallingPosteriorsReader(
            shared_workspace, warm_up_task.continuous_model, warm_up_task.continuous_model_approx,
            args.input_calls_path)()

    if hasattr(args, 'input_opt_path'):
        logger.info("A saved optimizer state was provided to use as starting point...")
        warm_up_task.fancy_opt.load(args.input_opt_path)


    try:
        # go!
        warm_up_task.engage()
        warm_up_task.disengage()
    except gcnvkernel.ConvergenceError as err:
        logger.info(err.message)
        # if inference diverged, pass an exit code to the Java side indicating that restart is needed
        sys.exit(gcnvkernel.io_consts.diverged_inference_exit_code)


    # main task
    main_task = gcnvkernel.CohortDenoisingAndCallingMainTask(
        main_denoising_config, main_calling_config, main_inference_params,
        shared_workspace, initial_params_supplier, warm_up_task)

    try:
        # go!
        main_task.engage()
        main_task.disengage()
    except gcnvkernel.ConvergenceError as err:
        logger.info(err.message)
        # if inference diverged, pass an exit code to the Java side indicating that restart is needed
        sys.exit(gcnvkernel.io_consts.diverged_inference_exit_code)


    # save model
    gcnvkernel.io_denoising_calling.DenoisingModelWriter(
        main_denoising_config, main_calling_config,
        shared_workspace, main_task.continuous_model, main_task.continuous_model_approx,
        args.output_model_path)()

    # save calls
    gcnvkernel.io_denoising_calling.SampleDenoisingAndCallingPosteriorsWriter(
        main_denoising_config, main_calling_config, shared_workspace, main_task.continuous_model,
        main_task.continuous_model_approx, args.output_calls_path)()

    # save optimizer state
    if hasattr(args, 'output_opt_path'):
        main_task.fancy_opt.save(args.output_opt_path)

    # save ELBO history
    if hasattr(args, 'output_tracking_path'):
        gcnvkernel.io_commons.assert_output_path_writable(args.output_tracking_path)

        warm_up_elbo_hist_file = os.path.join(args.output_tracking_path, "warm_up_elbo_history.tsv")
        warm_up_task.save_elbo_history(warm_up_elbo_hist_file)

        main_elbo_hist_file = os.path.join(args.output_tracking_path, "main_elbo_history.tsv")
        main_task.save_elbo_history(main_elbo_hist_file)
