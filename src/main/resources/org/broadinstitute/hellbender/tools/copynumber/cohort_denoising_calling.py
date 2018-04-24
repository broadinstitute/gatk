import os
import shutil

# set theano flags
os.environ["THEANO_FLAGS"] = "device=cpu,floatX=float64,optimizer=fast_run,compute_test_value=ignore,openmp=true"

import logging
import argparse
import gcnvkernel

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

    # instantiate config classes
    denoising_config = gcnvkernel.DenoisingModelConfig.from_args_dict(args_dict)
    calling_config = gcnvkernel.CopyNumberCallingConfig.from_args_dict(args_dict)
    inference_params = gcnvkernel.HybridInferenceParameters.from_args_dict(args_dict)

    # instantiate and initialize the workspace
    shared_workspace = gcnvkernel.DenoisingCallingWorkspace(
        denoising_config, calling_config, modeling_interval_list,
        n_st, sample_names, sample_metadata_collection)
    shared_workspace.initialize_copy_number_class_inference_vars()

    initial_params_supplier = gcnvkernel.DefaultDenoisingModelInitializer(
        denoising_config, calling_config, shared_workspace)

    task = gcnvkernel.CohortDenoisingAndCallingTask(
        denoising_config, calling_config, inference_params,
        shared_workspace, initial_params_supplier)

    if hasattr(args, 'input_model_path'):
        logger.info("A model path was provided to use as starting point...")
        gcnvkernel.io_denoising_calling.DenoisingModelReader(
            denoising_config, calling_config, shared_workspace, task.continuous_model,
            task.continuous_model_approx, args.input_model_path)()

    if hasattr(args, 'input_calls_path'):
        logger.info("A call path was provided to use as starting point...")
        gcnvkernel.io_denoising_calling.SampleDenoisingAndCallingPosteriorsReader(
            shared_workspace, task.continuous_model, task.continuous_model_approx,
            args.input_calls_path)()

    if hasattr(args, 'input_opt_path'):
        logger.info("A saved optimizer state was provided to use as starting point...")
        task.fancy_opt.load(args.input_opt_path)

    # go!
    task.engage()
    task.disengage()

    # save model
    gcnvkernel.io_denoising_calling.DenoisingModelWriter(
        denoising_config, calling_config,
        shared_workspace, task.continuous_model, task.continuous_model_approx,
        args.output_model_path)()

    # save a copy of targets in the model path
    shutil.copy(args.modeling_interval_list,
                os.path.join(args.output_model_path, gcnvkernel.io_consts.default_interval_list_filename))

    # save calls
    gcnvkernel.io_denoising_calling.SampleDenoisingAndCallingPosteriorsWriter(
        denoising_config, calling_config, shared_workspace, task.continuous_model, task.continuous_model_approx,
        args.output_calls_path)()

    # save a copy of targets in the calls path
    shutil.copy(args.modeling_interval_list,
                os.path.join(args.output_calls_path, gcnvkernel.io_consts.default_interval_list_filename))

    # save optimizer state
    if hasattr(args, 'output_opt_path'):
        task.fancy_opt.save(args.output_opt_path)
