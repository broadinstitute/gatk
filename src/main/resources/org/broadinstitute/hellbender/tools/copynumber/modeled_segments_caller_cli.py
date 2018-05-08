import argparse
from modeled_segments_caller import LoadAndSampleCrAndAf
from modeled_segments_caller import ModeledSegmentsCaller



parser = argparse.ArgumentParser(description="gCNV contig ploidy and read depth determination tool")

# add tool-specific args
group = parser.add_argument_group(title="Required arguments")

group.add_argument("--input",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Input .seg file (which is an output of ModelSegments).")

group.add_argument("--output",
                   type=str,
                   required=False,
                   default=argparse.SUPPRESS,
                   help="Directory containing all output files.")

group.add_argument("--output_prefix",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Prefix of the output filenames.")

group.add_argument("--output_image_suffix",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Suffix of the output image filenames.")

group.add_argument("--output_calls_suffix",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Suffix of the output calls filenames.")

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

group.add_argument("--normal_minor_allele_fraction_threshold",
                   type=float,
                   required=False,
                   default=0.475,
                   help="If the allele fraction value of a peak fitted to the data is above this threshold and its "
                        "copy ratio value is within the appropriate region, then the peak is considered normal.")

group.add_argument("--copy_ratio_peak_min_weight",
                   type=float,
                   required=False,
                   default=0.03,
                   help="During the copy ratio clustering, peaks with weights smaller than this ratio are not "
                        "taken into account.")

group.add_argument("--min_fraction_of_points_in_normal_allele_fraction_region",
                   type=float,
                   required=False,
                   default=0.15,
                   help="The region of copy ratio values are is considered normal only if at least this "
                        "fraction of points are above the normalMinorAlleleFractionThreshold",)


def str2bool(boolString):
    if boolString=="true" or boolString==True:
        return True
    return False


if __name__ == "__main__":

    # parse arguments
    args = parser.parse_args()

    # Run code
    data = LoadAndSampleCrAndAf(args.input,
                                load_CR=str2bool(args.load_copy_ratio),
                                load_AF=str2bool(args.load_allele_fraction),
                                output_log_dir=args.output,
                                output_log_prefix=args.output_prefix
                                )

    caller = ModeledSegmentsCaller(data, interactive=args.interactive,
                       output_image_dir=args.output,
                       output_calls_dir=args.output,
                       output_image_prefix=args.output_prefix,
                       output_calls_prefix=args.output_prefix,
                       output_image_suffix=args.output_image_suffix,
                       output_calls_suffix=args.output_calls_suffix,
                       interactive_output_del_ampl_image_suffix=args.interactive_output_del_ampl_image_suffix,
                       interactive_output_scatter_plot_suffix=args.interactive_output_scatter_plot_suffix,
                       interactive_output_allele_fraction_plot_suffix=args.interactive_output_allele_fraction_plot_suffix,
                       interactive_output_copy_ratio_suffix=args.interactive_output_copy_ratio_suffix,
                       interactive_output_copy_ratio_clustering_suffix=args.interactive_output_copy_ratio_clustering_suffix,
                       normal_minor_allele_fraction_threshold=args.normal_minor_allele_fraction_threshold,
                       copy_ratio_peak_min_weight=args.copy_ratio_peak_min_weight,
                       min_fraction_of_points_in_normal_allele_fraction_region=args.min_fraction_of_points_in_normal_allele_fraction_region
                       )
