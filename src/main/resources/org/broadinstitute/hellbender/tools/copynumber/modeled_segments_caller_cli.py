import argparse
from modeled_segments_caller import LoadAndSampleCrAndAf
from modeled_segments_caller import ModeledSegmentsCaller

parser = argparse.ArgumentParser(description="CNV caller tool")

# add tool-specific args
group = parser.add_argument_group(title="Required arguments")

group.add_argument("--input",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Input modelFinal.seg file (which is an output of ModelSegments).")

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
                   help="Suffix of the output image filename.")

group.add_argument("--output_calls_suffix",
                   type=str,
                   required=True,
                   default=argparse.SUPPRESS,
                   help="Suffix of the output calls filename.")

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

group.add_argument("--weight_ratio_max",
                   type=float,
                   required=False,
                   default=10,
                   help="Upon loading the data from file, we set an upper cut-off to the segment weights, "
                        "given by weight_ratio_max * mean(weights). This is required as ModeledSegments might "
                        "over-estimate these weights.")

group.add_argument("--log",
                   type=str,
                   required=False,
                   default="True",
                   help="Whether to log the progress of the code.")

group.add_argument("--interactive",
                   type=str,
                   required=False,
                   default="True",
                   help="Whether auxiliary plots should be saved.")

group.add_argument("--interactive_output_calls_image_suffix",
                   type=str,
                   required=False,
                   default="",
                   help="[Interactive mode only:] Suffix of plot showing deletions and insertions across the genome.")

group.add_argument("--interactive_output_summary_plot_suffix",
                   type=str,
                   required=False,
                   default="",
                   help="[Interactive mode only:] Suffix of plot showing the histograms of copy ratio and " +
                        "allele fraction data as well as scatter plot of normal and not normal segments " +
                        "in copy ratio and allele fraction spaces.")

group.add_argument("--interactive_output_allele_fraction_plot_suffix",
                   type=str,
                   required=False,
                   default="",
                   help="[Interactive mode only:] Suffix of the plot showing allele fraction histograms of the " +
                        "segments whose copy ratio falls within the regions where we think normal segments should be.")

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
                        "the two segments that contain the clusters of the two lowest average copy ratio values.")

group.add_argument("--normal_minor_allele_fraction_threshold",
                   type=float,
                   required=False,
                   default=0.475,
                   help="If the allele fraction value of a peak fitted to the data is above this threshold and its "
                        "copy ratio value is within the appropriate region, then the peak is considered normal.")

group.add_argument("--copy_ratio_peak_min_relative_height",
                   type=float,
                   required=False,
                   default=0.02,
                   help="During the copy ratio clustering, peaks with weights smaller than this ratio are not "
                        "taken into account.")

group.add_argument("--copy_ratio_kernel_density_bandwidth",
                   type=float,
                   required=False,
                   default=0.015,
                   help="During the copy ratio clustering, we smoothen the data using a Gaussian kernel of "
                        "this bandwidth. If it is set to 0 or a negative value, its value is set automatically, "
                        "using the data.")

group.add_argument("--min_weight_first_cr_peak_cr_data_only",
                   type=float,
                   required=False,
                   default=0.40,
                   help="If only copy ratio data is taken into account, and we find more than one cluster in the "
                        "data, then the first peak is considered normal if its relative weight compared to the second"
                        "peak is above this threshold (or if the weight of the second peak is smaller than 5%.")

group.add_argument("--min_fraction_of_points_in_normal_allele_fraction_region",
                   type=float,
                   required=False,
                   default=0.05,
                   help="The region of copy ratio values are is considered normal only if at least this "
                        "fraction of points are above the normalMinorAlleleFractionThreshold",)

group.add_argument("--responsibility_threshold_normal",
                   type=float,
                   required=False,
                   default=0.5,
                   help="Segments are considered normal if the responsibility in the Gaussian mixture model "
                        "of them being normal exceeds this threshold (set to 0.5 by default).")

group.add_argument("--max_phred_score_normal",
                   type=float,
                   required=False,
                   default=100,
                   help="Upper cut-off for the PHRED scores on normal samples, to make sure that they are not off the "
                        "scale on the plots.")

group.add_argument("--n_inference_iterations",
                   type=int,
                   required=False,
                   default=120000,
                   help="Maximum number of iterations of the automatic differentiation variational inference "
                        "fitting a Gaussian mixture model to the segments in order to cluster them.")

group.add_argument("--inference_total_grad_norm_constraint",
                   type=float,
                   required=False,
                   default=50.,
                   help="In order to ensure convergence, this constant limits the value of the gradient in automatic"
                        "differentiation variational inference, fitting a Gaussian mixture model to the segments.")

group.add_argument("--n_extra_Gaussians_mixture_model",
                   type=int,
                   required=False,
                   default=12,
                   help="The number of Gaussians needed to classify segments is estimated using k-means clustering. "
                        "This parameters lets us add a few extra Gaussians that might cover some of the noise or "
                        "additional subcluster structure present in the data.")

group.add_argument("--max_n_peaks_in_copy_ratio",
                   type=int,
                   required=False,
                   default=10,
                   help="In order to discover the different copy number states, the clusters in the copy ratio data "
                        "is determined. This constant limits the maximum number of peaks that can be found in the data.")

group.add_argument("--gaussian_prior_standard_deviation",
                   type=int,
                   required=False,
                   default=0.0005,
                   help="The standard deviation of the prior of the mean and the standard deviation of the Gaussians "
                        "that we fit to the segment data in (copy ratio, allele fraction) space. Since the positions "
                        "of the peaks are estimated using k-means clustering, the prior around these points should be "
                        "tight.")

def str2bool(boolString):
    if boolString=="true" or boolString==True:
        return True
    return False


if __name__ == "__main__":

    # parse arguments
    args = parser.parse_args()

    # Load data from the input file and sample data points
    data = LoadAndSampleCrAndAf(args.input,
                                load_cr=str2bool(args.load_copy_ratio),
                                load_af=str2bool(args.load_allele_fraction),
                                output_log_dir=args.output,
                                output_log_prefix=args.output_prefix,
                                do_logging=args.log,
                                weight_ratio_max=args.weight_ratio_max
                                )

    # Run the caller
    caller = ModeledSegmentsCaller(data,
                                   interactive=args.interactive,
                                   output_image_dir=args.output,
                                   output_calls_dir=args.output,
                                   output_image_prefix=args.output_prefix,
                                   output_calls_prefix=args.output_prefix,
                                   output_image_suffix=args.output_image_suffix,
                                   output_calls_suffix=args.output_calls_suffix,
                                   interactive_output_calls_image_suffix=args.interactive_output_calls_image_suffix,
                                   interactive_output_summary_plot_suffix=args.interactive_output_summary_plot_suffix,
                                   interactive_output_allele_fraction_plot_suffix=args.interactive_output_allele_fraction_plot_suffix,
                                   interactive_output_copy_ratio_suffix=args.interactive_output_copy_ratio_suffix,
                                   interactive_output_copy_ratio_clustering_suffix=args.interactive_output_copy_ratio_clustering_suffix,
                                   normal_minor_allele_fraction_threshold=args.normal_minor_allele_fraction_threshold,
                                   copy_ratio_peak_min_relative_height=args.copy_ratio_peak_min_relative_height,
                                   copy_ratio_kernel_density_bandwidth=args.copy_ratio_kernel_density_bandwidth,
                                   min_weight_first_cr_peak_cr_data_only=args.min_weight_first_cr_peak_cr_data_only,
                                   min_fraction_of_points_in_normal_allele_fraction_region=args.min_fraction_of_points_in_normal_allele_fraction_region,
                                   responsibility_threshold_normal=args.responsibility_threshold_normal,
                                   max_phred_score_normal=args.max_phred_score_normal,
                                   n_inference_iterations=args.n_inference_iterations,
                                   inference_total_grad_norm_constraint=args.inference_total_grad_norm_constraint,
                                   n_extra_Gaussians_mixture_model=args.n_extra_Gaussians_mixture_model,
                                   max_n_peaks_in_copy_ratio=args.max_n_peaks_in_copy_ratio,
                                   gaussian_prior_standard_deviation=args.gaussian_prior_standard_deviation
                                   )
