import argparse
from modeled_segments_caller import LoadAndSampleCrAndAf
from modeled_segments_caller import ModeledSegmentsCaller



parser = argparse.ArgumentParser(description="gCNV contig ploidy and read depth determination tool")

# add tool-specific args
group = parser.add_argument_group(title="Required arguments")

group.add_argument("--input_file",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Input .seg file (which is an output of ModelSegments).")

group.add_argument("--output_image_dir",
                   type=str,
                   required=False,
                   default=argparse.SUPPRESS,
                   help="Output path to the image file showing the plots of the segments.")

group.add_argument("--output_calls_dir",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Output path to the file containing the list of called intervals.")

group.add_argument("--output_image_prefix",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Prefix of the output image filenames.")

group.add_argument("--output_calls_prefix",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Prefix of the output calls filenames.")

group.add_argument("--load_copy_ratio",
                   type=str,
                   required=False,
                   default="True",
                   help="Whether to load copy ratio data.")

group.add_argument("--load_allele_fraction",
                   type=str,
                   required=False,
                   default="True",
                   help="Whether to load allele fraction data.")

group.add_argument("--interactive",
                   type=str,
                   required=False,
                   default="True",
                   help="Whether auxiliary plots should be saved.")

group.add_argument("--interactive_output_del_ampl_image_suffix",
                   type=str,
                   required=False,
                   default="",
                   help="[Interactive mode only:] Suffix of plot showing deletions and insertions across the genome.")

group.add_argument("--interactive_output_scatter_plot_suffix",
                   type=str,
                   required=False,
                   default="",
                   help="[Interactive mode only:] Suffix of plot showing the normal and not normal segments on" +
                        " a scatter plot in copy ratio and allele fraction spaces.")

group.add_argument("--interactive_output_allele_fraction_plot_suffix",
                   type=str,
                   required=False,
                   default="",
                   help="[Interactive mode only:] Suffix of the plot showing allele fraction histograms of the" +
                        " data falling in the first and second segment of the copy ratio data.")

group.add_argument("--interactive_output_copy_ratio_suffix",
                   type=str,
                   required=False,
                   default="",
                   help="[Interactive mode only:] Suffix of the plot showing the histogram and the Gaussian fits to" +
                        " copy ratio data.")

group.add_argument("--interactive_output_copy_ratio_clustering_suffix",
                   type=str,
                   required=False,
                   default="",
                   help="[Interactive mode only:] Suffix of the plot showing the histogram of copy ratio data and" +
                        "the first two copy ratio clusters.")

group.add_argument("--log_filename_prefix",
                   type=str,
                   required=False,
                   default="",
                   help="Whether progress of the code should be saved in a log file.")


def str2bool(boolString):
    if boolString=="true" or boolString==True:
        return True
    return False


if __name__ == "__main__":

    # parse arguments
    args = parser.parse_args()

    # Run code
    data = LoadAndSampleCrAndAf(args.input_file,
                                load_CR=str2bool(args.load_copy_ratio),
                                load_AF=str2bool(args.load_allele_fraction),
                                log_filename_prefix=args.log_filename_prefix
                                )

    caller = ModeledSegmentsCaller(data, interactive=args.interactive,
                       output_image_dir=args.output_image_dir,
                       output_calls_dir=args.output_calls_dir,
                       output_image_prefix=args.output_image_prefix,
                       output_calls_prefix=args.output_calls_prefix,
                       interactive_output_del_ampl_image_suffix=args.interactive_output_del_ampl_image_suffix,
                       interactive_output_scatter_plot_suffix=args.interactive_output_scatter_plot_suffix,
                       interactive_output_allele_fraction_plot_suffix=args.interactive_output_allele_fraction_plot_suffix,
                       interactive_output_copy_ratio_suffix=args.interactive_output_copy_ratio_suffix,
                       interactive_output_copy_ratio_clustering_suffix=args.interactive_output_copy_ratio_clustering_suffix
                       )
