package org.broadinstitute.hellbender.tools.copynumber;

import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.math3.special.Beta;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.*;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SampleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.MultidimensionalSegment;
import org.broadinstitute.hellbender.tools.copynumber.models.AlleleFractionModeller;
import org.broadinstitute.hellbender.tools.copynumber.models.AlleleFractionPrior;
import org.broadinstitute.hellbender.tools.copynumber.models.CopyRatioModeller;
import org.broadinstitute.hellbender.tools.copynumber.models.MultidimensionalModeller;
import org.broadinstitute.hellbender.tools.copynumber.segmentation.AlleleFractionKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.segmentation.CopyRatioKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.segmentation.MultidimensionalKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Model segmented copy ratio from denoised read counts and segmented minor-allele fraction from allelic counts.
 *
 * <p>
 *     Possible inputs are: 1) denoised copy ratios for the case sample, 2) allelic counts for the case sample,
 *     and 3) allelic counts for a matched-normal sample.  All available inputs will be used to to perform
 *     segmentation and model inference.
 * </p>
 *
 * <p>
 *     If allelic counts are available, the first step in the inference process is to genotype heterozygous sites,
 *     as the allelic counts at these sites will subsequently be modeled to infer segmented minor-allele fraction.
 *     We perform a relatively simple and naive genotyping based on the allele counts (i.e., pileups), which is
 *     controlled by a small number of parameters ({@code minimum-total-allele-count},
 *     {@code genotyping-homozygous-log-ratio-threshold}, and {@code genotyping-homozygous-log-ratio-threshold}).
 *     If the matched normal is available, its allelic counts will be used to genotype the sites, and
 *     we will simply assume these genotypes are the same in the case sample.  (This can be critical, for example,
 *     for determining sites with loss of heterozygosity in high purity case samples; such sites will be genotyped as
 *     homozygous if the matched-normal sample is not available.)
 * </p>
 *
 * <p>
 *     Next, we segment, if available, the denoised copy ratios and the alternate-allele fractions at the
 *     genotyped heterozygous sites.  This is done using kernel segmentation (see {@link KernelSegmenter});
 *     if both copy ratios and allele fractions are available, we use a combined kernel that is sensitive
 *     to changes that occur in only either of the two or in both.  Various segmentation parameters control
 *     the sensitivity of the segmentation and should be selected appropriately for each analysis.
 * </p>
 *
 * <p>
 *     After segmentation is complete, we run Markov-chain Monte Carlo (MCMC) to determine posteriors for
 *     segmented models for the log2 copy ratio and the minor-allele fraction; see {@link CopyRatioModeller}
 *     and {@link AlleleFractionModeller}, respectively.  After the first run of MCMC is complete,
 *     smoothing of the segmented posteriors is performed by merging adjacent segments whose posterior
 *     credible intervals sufficiently overlap according to specified segmentation-smoothing parameters.
 *     Then, additional rounds of segmentation smoothing (with intermediate MCMC optionally performed in between rounds)
 *     are performed until convergence, at which point a final round of MCMC is performed.
 * </p>
 *
 * <h3>Input</h3>
 *
 * <ul>
 *     <li>
 *         (Optional) Denoised-copy-ratios file (output of {@link DenoiseReadCounts}.
 *         If allelic counts are not provided, then this is required.
 *     </li>
 *     <li>
 *         (Optional) Allelic-counts file (output of {@link CollectAllelicCounts}.
 *         If denoised copy ratios are not provided, then this is required.
 *     </li>
 *     <li>
 *         (Optional) Matched-normal allelic-counts file (output of {@link CollectAllelicCounts}.
 *         This can only be provided if allelic counts for the case sample are also provided.
 *     </li>
 *     <li>
 *         Output prefix.
 *         This is used as a prefix for output filenames.
 *     </li>
 *     <li>
 *         Output directory.
 *         This must be a pre-existing directory.
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Modeled-segments files (.modelBegin.seg and .modelFinal.seg).
 *         These are TSVs with a SAM-style header containing a read-group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link ModeledSegmentCollection.ModeledSegmentTableColumn},
 *         and the corresponding entry rows.
 *         The initial result before segmentation smoothing is output to the .modelBegin.seg file
 *         and the final result after segmentation smoothing is output to the .modelFinal.seg file.
 *     </li>
 *     <li>
 *         Allele-fraction-model global-parameter files (.modelBegin.af.param and .modelFinal.af.param).
 *         These are TSVs with a SAM-style header containing a read-group sample name,
 *         a row specifying the column headers contained in {@link ParameterDecileCollection.ParameterTableColumn},
 *         and the corresponding entry rows.
 *         The initial result before segmentation smoothing is output to the .modelBegin.af.param file
 *         and the final result after segmentation smoothing is output to the .modelFinal.af.param file.
 *     </li>
 *     <li>
 *         Copy-ratio-model global-parameter files (.modelBegin.cr.param and .modelFinal.cr.param).
 *         These are TSVs with a SAM-style header containing a read-group sample name,
 *         a row specifying the column headers contained in {@link ParameterDecileCollection.ParameterTableColumn},
 *         and the corresponding entry rows.
 *         The initial result before segmentation smoothing is output to the .modelBegin.cr.param file
 *         and the final result after segmentation smoothing is output to the .modelFinal.cr.param file.
 *     </li>
 *     <li>
 *         Copy-ratio segment file (.cr.param).
 *         This is a TSV with a SAM-style header containing a read-group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link CopyRatioSegmentCollection.CopyRatioSegmentTableColumn},
 *         and the corresponding entry rows.
 *         It contains the segments from the .modelFinal.seg file converted to a format suitable for input to {@link CallCopyRatioSegments}.
 *     </li>
 *     <li>
 *         (Optional) Allelic-counts file containing the counts at sites genotyped as heterozygous in the case sample (.hets.tsv).
 *         This is a TSV with a SAM-style header containing a read-group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link AllelicCountCollection.AllelicCountTableColumn},
 *         and the corresponding entry rows.
 *         This is only output if normal allelic counts are provided as input.
 *     </li>
 *     <li>
 *         (Optional) Allelic-counts file containing the counts at sites genotyped as heterozygous in the matched-normal sample (.hets.normal.tsv).
 *         This is a TSV with a SAM-style header containing a read-group sample name, a sequence dictionary,
 *         a row specifying the column headers contained in {@link AllelicCountCollection.AllelicCountTableColumn},
 *         and the corresponding entry rows.
 *         This is only output if matched-normal allelic counts are provided as input.
 *     </li>
 * </ul>
 *
 * <h3>Examples</h3>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --denoised-copy-ratios tumor.denoisedCR.tsv \
 *          --allelic-counts tumor.allelicCounts.tsv \
 *          --normal-allelic-counts normal.allelicCounts.tsv \
 *          --output-prefix tumor \
 *          -O output_dir
 * </pre>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --allelic-counts tumor.allelicCounts.tsv \
 *          --normal-allelic-counts normal.allelicCounts.tsv \
 *          --output-prefix tumor \
 *          -O output_dir
 * </pre>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --denoised-copy-ratios normal.denoisedCR.tsv \
 *          --output-prefix normal \
 *          -O output_dir
 * </pre>
 *
 * <pre>
 *     gatk ModelSegments \
 *          --denoised-copy-ratios normal.denoisedCR.tsv \
 *          --allelic-counts normal.allelicCounts.tsv \
 *          --output-prefix normal \
 *          -O output_dir
 * </pre>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Model segmented copy ratio from denoised read counts and segmented minor-allele fraction from allelic counts.",
        oneLineSummary = "Model segmented copy ratio from denoised read counts and segmented minor-allele fraction from allelic counts.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class ModelSegments extends CommandLineProgram {
    //filename tags for output
    public static final String HET_ALLELIC_COUNTS_FILE_SUFFIX = ".hets.tsv";
    public static final String NORMAL_HET_ALLELIC_COUNTS_FILE_SUFFIX = ".hets.normal.tsv";
    public static final String SEGMENTS_FILE_SUFFIX = ".seg";
    public static final String BEGIN_FIT_FILE_TAG = ".modelBegin";
    public static final String FINAL_FIT_FILE_TAG = ".modelFinal";
    public static final String COPY_RATIO_MODEL_PARAMETER_FILE_SUFFIX = ".cr.param";
    public static final String ALLELE_FRACTION_MODEL_PARAMETER_FILE_SUFFIX = ".af.param";
    public static final String COPY_RATIO_SEGMENTS_FOR_CALLER_FILE = ".cr" + SEGMENTS_FILE_SUFFIX;

    //het genotyping argument names
    public static final String MINIMUM_TOTAL_ALLELE_COUNT_LONG_NAME = "minimum-total-allele-count";
    public static final String GENOTYPING_HOMOZYGOUS_LOG_RATIO_THRESHOLD_LONG_NAME = "genotyping-homozygous-log-ratio-threshold";
    public static final String GENOTYPING_BASE_ERROR_RATE_LONG_NAME = "genotyping-base-error-rate";

    //segmentation argument names
    public static final String MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_LONG_NAME = "maximum-number-of-segments-per-chromosome";
    public static final String KERNEL_VARIANCE_COPY_RATIO_LONG_NAME = "kernel-variance-copy-ratio";
    public static final String KERNEL_VARIANCE_ALLELE_FRACTION_LONG_NAME = "kernel-variance-allele-fraction";
    public static final String KERNEL_SCALING_ALLELE_FRACTION_LONG_NAME = "kernel-scaling-allele-fraction";
    public static final String KERNEL_APPROXIMATION_DIMENSION_LONG_NAME = "kernel-approximation-dimension";
    public static final String WINDOW_SIZE_LONG_NAME = "window-size";
    public static final String NUMBER_OF_CHANGEPOINTS_PENALTY_FACTOR_LONG_NAME = "number-of-changepoints-penalty-factor";

    //MCMC argument names
    public static final String MINOR_ALLELE_FRACTION_PRIOR_ALPHA_LONG_NAME = "minor-allele-fraction-prior-alpha";
    public static final String NUMBER_OF_SAMPLES_COPY_RATIO_LONG_NAME = "number-of-samples-copy-ratio";
    public static final String NUMBER_OF_BURN_IN_SAMPLES_COPY_RATIO_LONG_NAME = "number-of-burn-in-samples-copy-ratio";
    public static final String NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME = "number-of-samples-allele-fraction";
    public static final String NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME = "number-of-burn-in-samples-allele-fraction";

    //smoothing argument names
    public static final String SMOOTHING_CREDIBLE_INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME = "smoothing-credible-interval-threshold-copy-ratio";
    public static final String SMOOTHING_CREDIBLE_INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME = "smoothing-credible-interval-threshold-allele-fraction";
    public static final String MAXIMUM_NUMBER_OF_SMOOTHING_ITERATIONS_LONG_NAME = "maximum-number-of-smoothing-iterations";
    public static final String NUMBER_OF_SMOOTHING_ITERATIONS_PER_FIT_LONG_NAME = "number-of-smoothing-iterations-per-fit";

    @Argument(
            doc = "Input file containing denoised copy ratios (output of DenoiseReadCounts).",
            fullName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME,
            optional = true
    )
    private File inputDenoisedCopyRatiosFile = null;

    @Argument(
            doc = "Input file containing allelic counts (output of CollectAllelicCounts).",
            fullName = CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME,
            optional = true
    )
    private File inputAllelicCountsFile = null;

    @Argument(
            doc = "Input file containing allelic counts for a matched normal (output of CollectAllelicCounts).",
            fullName = CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME,
            optional = true
    )
    private File inputNormalAllelicCountsFile = null;

    @Argument(
            doc = "Prefix for output files.",
            fullName =  CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME
    )
    private String outputPrefix;

    @Argument(
            doc = "Output directory.",
            fullName =  StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputDir;

    @Argument(
            doc = "Minimum total count for filtering allelic counts, if available.",
            fullName = MINIMUM_TOTAL_ALLELE_COUNT_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int minTotalAlleleCount = 30;

    @Argument(
            doc = "Log-ratio threshold for genotyping and filtering homozygous allelic counts, if available.  " +
                    "Increasing this value will increase the number of sites assumed to be heterozygous for modeling.",
            fullName = GENOTYPING_HOMOZYGOUS_LOG_RATIO_THRESHOLD_LONG_NAME,
            optional = true
    )
    private double genotypingHomozygousLogRatioThreshold = -10.;

    @Argument(
            doc = "Maximum base-error rate for genotyping and filtering homozygous allelic counts, if available.  " +
                    "The likelihood for an allelic count to be generated from a homozygous site will be integrated " +
                    "from zero base-error rate up to this value.  Decreasing this value will increase " +
                    "the number of sites assumed to be heterozygous for modeling.",
            fullName = GENOTYPING_BASE_ERROR_RATE_LONG_NAME,
            optional = true
    )
    private double genotypingBaseErrorRate = 5E-2;

    @Argument(
            doc = "Maximum number of segments allowed per chromosome.",
            fullName = MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_LONG_NAME,
            minValue = 1,
            optional = true
    )
    private int maxNumSegmentsPerChromosome = 1000;

    @Argument(
            doc = "Variance of Gaussian kernel for copy-ratio segmentation, if performed.  If zero, a linear kernel will be used.",
            fullName = KERNEL_VARIANCE_COPY_RATIO_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double kernelVarianceCopyRatio = 0.;

    @Argument(
            doc = "Variance of Gaussian kernel for allele-fraction segmentation, if performed.  If zero, a linear kernel will be used.",
            fullName = KERNEL_VARIANCE_ALLELE_FRACTION_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double kernelVarianceAlleleFraction = 0.01;

    @Argument(
            doc = "Relative scaling S of the kernel K_AF for allele-fraction segmentation to the kernel K_CR for copy-ratio segmentation.  " +
                    "If multidimensional segmentation is performed, the total kernel used will be K_CR + S * K_AF.",
            fullName = KERNEL_SCALING_ALLELE_FRACTION_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double kernelScalingAlleleFraction = 1.0;

    @Argument(
            doc = "Dimension of the kernel approximation.  A subsample containing this number of data points " +
                    "will be used to construct the approximation for each chromosome.  " +
                    "If the total number of data points in a chromosome is greater " +
                    "than this number, then all data points in the chromosome will be used.  " +
                    "Time complexity scales quadratically and space complexity scales linearly with this parameter.",
            fullName = KERNEL_APPROXIMATION_DIMENSION_LONG_NAME,
            minValue = 1,
            optional = true
    )
    private int kernelApproximationDimension = 100;

    @Argument(
            doc = "Window sizes to use for calculating local changepoint costs.  " +
                    "For each window size, the cost for each data point to be a changepoint will be calculated " +
                    "assuming that the point demarcates two adjacent segments of that size.  " +
                    "Including small (large) window sizes will increase sensitivity to small (large) events.  " +
                    "Duplicate values will be ignored.",
            fullName = WINDOW_SIZE_LONG_NAME,
            minValue = 1,
            optional = true
    )
    private List<Integer> windowSizes = new ArrayList<>(Arrays.asList(8, 16, 32, 64, 128, 256));

    @Argument(
            doc = "Factor A for the penalty on the number of changepoints per chromosome for segmentation.  " +
                    "Adds a penalty of the form A * C * [1 + log (N / C)], " +
                    "where C is the number of changepoints in the chromosome, " +
                    "to the cost function for each chromosome.  " +
                    "Must be non-negative.",
            fullName = NUMBER_OF_CHANGEPOINTS_PENALTY_FACTOR_LONG_NAME,
            minValue = 0.,
            optional = true
    )
    private double numChangepointsPenaltyFactor = 1.;

    @Argument(
            doc = "Alpha hyperparameter for the 4-parameter beta-distribution prior on segment minor-allele fraction. " +
                    "The prior for the minor-allele fraction f in each segment is assumed to be Beta(alpha, 1, 0, 1/2). " +
                    "Increasing this hyperparameter will reduce the effect of reference bias at the expense of sensitivity.",
            fullName = MINOR_ALLELE_FRACTION_PRIOR_ALPHA_LONG_NAME,
            optional = true,
            minValue = 1
    )
    private double minorAlleleFractionPriorAlpha = 25.;

    @Argument(
            doc = "Total number of MCMC samples for copy-ratio model.",
            fullName = NUMBER_OF_SAMPLES_COPY_RATIO_LONG_NAME,
            optional = true,
            minValue = 1
    )
    private int numSamplesCopyRatio = 100;

    @Argument(
            doc = "Number of burn-in samples to discard for copy-ratio model.",
            fullName = NUMBER_OF_BURN_IN_SAMPLES_COPY_RATIO_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int numBurnInCopyRatio = 50;

    @Argument(
            doc = "Total number of MCMC samples for allele-fraction model.",
            fullName = NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME,
            optional = true,
            minValue = 1
    )
    private int numSamplesAlleleFraction = 100;

    @Argument(
            doc = "Number of burn-in samples to discard for allele-fraction model.",
            fullName = NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int numBurnInAlleleFraction = 50;

    @Argument(
            doc = "Number of 10% equal-tailed credible-interval widths to use for copy-ratio segmentation smoothing.",
            fullName = SMOOTHING_CREDIBLE_INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME,
            optional = true,
            minValue = 0.
    )
    private double smoothingCredibleIntervalThresholdCopyRatio = 2.;

    @Argument(
            doc = "Number of 10% equal-tailed credible-interval widths to use for allele-fraction segmentation smoothing.",
            fullName = SMOOTHING_CREDIBLE_INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME,
            optional = true,
            minValue = 0.
    )
    private double smoothingCredibleIntervalThresholdAlleleFraction = 2.;

    @Argument(
            doc = "Maximum number of iterations allowed for segmentation smoothing.",
            fullName = MAXIMUM_NUMBER_OF_SMOOTHING_ITERATIONS_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int maxNumSmoothingIterations = 25;

    @Argument(
            doc = "Number of segmentation-smoothing iterations per MCMC model refit. " +
                    "(Increasing this will decrease runtime, but the final number of segments may be higher. " +
                    "Setting this to 0 will completely disable model refitting between iterations.)",
            fullName = NUMBER_OF_SMOOTHING_ITERATIONS_PER_FIT_LONG_NAME,
            optional = true,
            minValue = 0
    )
    private int numSmoothingIterationsPerFit = 0;

    //initialize data variables, some of which may be optional
    private CopyRatioCollection denoisedCopyRatios = null;
    private AllelicCountCollection hetAllelicCounts = null;

    @Override
    protected Object doWork() {
        validateArguments();

        //perform one-dimensional or multidimensional segmentation as appropriate and write to file
        //(for use by CallCopyRatioSegments, if copy ratios are available)
        final MultidimensionalSegmentCollection multidimensionalSegments;
        if (inputDenoisedCopyRatiosFile != null && inputAllelicCountsFile == null) {
            readDenoisedCopyRatios();
            final CopyRatioSegmentCollection copyRatioSegments = performCopyRatioSegmentation();
            multidimensionalSegments = new MultidimensionalSegmentCollection(
                    copyRatioSegments.getMetadata(),
                    copyRatioSegments.getRecords().stream()
                            .map(s -> new MultidimensionalSegment(s.getInterval(), s.getNumPoints(), 0, s.getMeanLog2CopyRatio()))
                            .collect(Collectors.toList()));
            hetAllelicCounts = new AllelicCountCollection(denoisedCopyRatios.getMetadata(), Collections.emptyList()); //create an empty collection with the appropriate name
        } else if (inputDenoisedCopyRatiosFile == null && inputAllelicCountsFile != null) {
            readAndFilterAllelicCounts();
            final AlleleFractionSegmentCollection alleleFractionSegments = performAlleleFractionSegmentation();
            multidimensionalSegments = new MultidimensionalSegmentCollection(
                    alleleFractionSegments.getMetadata(),
                    alleleFractionSegments.getRecords().stream()
                            .map(s -> new MultidimensionalSegment(s.getInterval(), 0, s.getNumPoints(), Double.NaN))
                            .collect(Collectors.toList()));
            denoisedCopyRatios = new CopyRatioCollection(hetAllelicCounts.getMetadata(), Collections.emptyList());     //create an empty collection with the appropriate name
        } else {
            readDenoisedCopyRatios();
            readAndFilterAllelicCounts();
            multidimensionalSegments = new MultidimensionalKernelSegmenter(denoisedCopyRatios, hetAllelicCounts)
                    .findSegmentation(maxNumSegmentsPerChromosome,
                            kernelVarianceCopyRatio, kernelVarianceAlleleFraction, kernelScalingAlleleFraction, kernelApproximationDimension,
                            ImmutableSet.copyOf(windowSizes).asList(),
                            numChangepointsPenaltyFactor, numChangepointsPenaltyFactor);
        }

        logger.info("Modeling available denoised copy ratios and heterozygous allelic counts...");
        //initial MCMC model fitting performed by MultidimensionalModeller constructor
        final AlleleFractionPrior alleleFractionPrior = new AlleleFractionPrior(minorAlleleFractionPriorAlpha);
        final MultidimensionalModeller modeller = new MultidimensionalModeller(
                multidimensionalSegments, denoisedCopyRatios, hetAllelicCounts, alleleFractionPrior,
                numSamplesCopyRatio, numBurnInCopyRatio,
                numSamplesAlleleFraction, numBurnInAlleleFraction);

        //write initial segments and parameters to file
        writeModeledSegmentsAndParameterFiles(modeller, BEGIN_FIT_FILE_TAG);

        //segmentation smoothing
        modeller.smoothSegments(
                maxNumSmoothingIterations, numSmoothingIterationsPerFit,
                smoothingCredibleIntervalThresholdCopyRatio, smoothingCredibleIntervalThresholdAlleleFraction);

        //write final segments and parameters to file
        writeModeledSegmentsAndParameterFiles(modeller, FINAL_FIT_FILE_TAG);

        //write final segments for copy-ratio caller (TODO remove this and MEAN_LOG2_COPY_RATIO column when new caller is available)
        final OverlapDetector<CopyRatio> copyRatioMidpointOverlapDetector = denoisedCopyRatios.getMidpointOverlapDetector();
        final CopyRatioSegmentCollection copyRatioSegmentsFinal = new CopyRatioSegmentCollection(
                modeller.getModeledSegments().getMetadata(),
                modeller.getModeledSegments().getIntervals().stream()
                        .map(s -> new CopyRatioSegment(s, new ArrayList<>(copyRatioMidpointOverlapDetector.getOverlaps(s))))
                        .collect(Collectors.toList()));
        writeSegments(copyRatioSegmentsFinal, COPY_RATIO_SEGMENTS_FOR_CALLER_FILE);

        logger.info("SUCCESS: ModelSegments run complete.");

        return "SUCCESS";
    }

    private void validateArguments() {
        Utils.nonNull(outputPrefix);
        Utils.validateArg(!(inputDenoisedCopyRatiosFile == null && inputAllelicCountsFile == null),
                "Must provide at least a denoised-copy-ratios file or an allelic-counts file.");
        Utils.validateArg(!(inputAllelicCountsFile == null && inputNormalAllelicCountsFile != null),
                "Must provide an allelic-counts file for the case sample to run in matched-normal mode.");
        if (inputDenoisedCopyRatiosFile != null) {
            IOUtils.canReadFile(inputDenoisedCopyRatiosFile);
        }
        if (inputAllelicCountsFile != null) {
            IOUtils.canReadFile(inputAllelicCountsFile);
        }
        if (inputNormalAllelicCountsFile != null) {
            IOUtils.canReadFile(inputNormalAllelicCountsFile);
        }
        if (!new File(outputDir).exists()) {
            throw new UserException(String.format("Output directory %s does not exist.", outputDir));
        }
    }

    private void readDenoisedCopyRatios() {
        logger.info(String.format("Reading denoised-copy-ratios file (%s)...", inputDenoisedCopyRatiosFile));
        denoisedCopyRatios = new CopyRatioCollection(inputDenoisedCopyRatiosFile);
    }

    private CopyRatioSegmentCollection performCopyRatioSegmentation() {
        logger.info("Starting segmentation of denoised copy ratios...");
        final int maxNumChangepointsPerChromosome = maxNumSegmentsPerChromosome - 1;
        return new CopyRatioKernelSegmenter(denoisedCopyRatios)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceCopyRatio, kernelApproximationDimension,
                        ImmutableSet.copyOf(windowSizes).asList(),
                        numChangepointsPenaltyFactor, numChangepointsPenaltyFactor);
    }

    private void readAndFilterAllelicCounts() {
        //read in case sample
        logger.info(String.format("Reading allelic-counts file (%s)...", inputAllelicCountsFile));
        final AllelicCountCollection unfilteredAllelicCounts = new AllelicCountCollection(inputAllelicCountsFile);
        final SampleLocatableMetadata metadata = unfilteredAllelicCounts.getMetadata();

        //filter on total count in case sample
        logger.info(String.format("Filtering allelic counts with total count less than %d...", minTotalAlleleCount));
        AllelicCountCollection filteredAllelicCounts = new AllelicCountCollection(
                metadata,
                unfilteredAllelicCounts.getRecords().stream()
                        .filter(ac -> ac.getTotalReadCount() >= minTotalAlleleCount)
                        .collect(Collectors.toList()));
        logger.info(String.format("Retained %d / %d sites after filtering on total count...",
                filteredAllelicCounts.getRecords().size(), unfilteredAllelicCounts.getRecords().size()));

        //filter on overlap with copy-ratio intervals, if available
        if (denoisedCopyRatios != null) {
            logger.info("Filtering allelic-count sites not overlapping with copy-ratio intervals...");
            final OverlapDetector<CopyRatio> copyRatioOverlapDetector = denoisedCopyRatios.getOverlapDetector();
            filteredAllelicCounts = new AllelicCountCollection(
                    metadata,
                    filteredAllelicCounts.getRecords().stream()
                            .filter(copyRatioOverlapDetector::overlapsAny)
                            .collect(Collectors.toList()));
            logger.info(String.format("Retained %d / %d sites after filtering on overlap with copy-ratio intervals...",
                    filteredAllelicCounts.getRecords().size(), unfilteredAllelicCounts.getRecords().size()));
        }

        if (inputNormalAllelicCountsFile == null) {
            //filter on homozygosity in case sample
            logger.info("No matched normal was provided, not running in matched-normal mode...");
            logger.info("Performing binomial testing and filtering homozygous allelic counts...");
            hetAllelicCounts = new AllelicCountCollection(
                    metadata,
                    filteredAllelicCounts.getRecords().stream()
                            .filter(ac -> calculateHomozygousLogRatio(ac, genotypingBaseErrorRate) < genotypingHomozygousLogRatioThreshold)
                            .collect(Collectors.toList()));
            final File hetAllelicCountsFile = new File(outputDir, outputPrefix + HET_ALLELIC_COUNTS_FILE_SUFFIX);
            hetAllelicCounts.write(hetAllelicCountsFile);
            logger.info(String.format("Retained %d / %d sites after testing for heterozygosity...",
                    hetAllelicCounts.getRecords().size(), unfilteredAllelicCounts.getRecords().size()));
            logger.info(String.format("Heterozygous allelic counts written to %s.", hetAllelicCountsFile));
        } else {
            //read in matched normal
            logger.info("Matched normal was provided, running in matched-normal mode...");
            logger.info("Performing binomial testing and filtering homozygous allelic counts in matched normal...");
            final AllelicCountCollection unfilteredNormalAllelicCounts = new AllelicCountCollection(inputNormalAllelicCountsFile);
            if (!unfilteredNormalAllelicCounts.getIntervals().equals(unfilteredAllelicCounts.getIntervals())) {
                throw new UserException.BadInput("Allelic-count sites in case sample and matched normal do not match. " +
                        "Run CollectAllelicCounts using the same interval list of sites for both samples.");
            }
            final SampleLocatableMetadata normalMetadata = unfilteredNormalAllelicCounts.getMetadata();

            //filter on total count in matched normal
            logger.info(String.format("Filtering allelic counts in matched normal with total count less than %d...", minTotalAlleleCount));
            final AllelicCountCollection filteredNormalAllelicCounts = new AllelicCountCollection(
                    normalMetadata,
                    unfilteredNormalAllelicCounts.getRecords().stream()
                            .filter(ac -> ac.getTotalReadCount() >= minTotalAlleleCount)
                            .collect(Collectors.toList()));
            logger.info(String.format("Retained %d / %d sites in matched normal after filtering on total count...",
                    filteredNormalAllelicCounts.getRecords().size(), unfilteredNormalAllelicCounts.getRecords().size()));

            //filter on homozygosity in matched normal
            final AllelicCountCollection hetNormalAllelicCounts = new AllelicCountCollection(
                    normalMetadata,
                    filteredNormalAllelicCounts.getRecords().stream()
                            .filter(ac -> calculateHomozygousLogRatio(ac, genotypingBaseErrorRate) < genotypingHomozygousLogRatioThreshold)
                            .collect(Collectors.toList()));
            final File hetNormalAllelicCountsFile = new File(outputDir, outputPrefix + NORMAL_HET_ALLELIC_COUNTS_FILE_SUFFIX);
            hetNormalAllelicCounts.write(hetNormalAllelicCountsFile);
            logger.info(String.format("Retained %d / %d sites in matched normal after testing for heterozygosity...",
                    hetNormalAllelicCounts.getRecords().size(), unfilteredNormalAllelicCounts.getRecords().size()));
            logger.info(String.format("Heterozygous allelic counts for matched normal written to %s.", hetNormalAllelicCountsFile));

            //retrieve sites in case sample
            logger.info("Retrieving allelic counts at these sites in case sample...");
            final Set<SimpleInterval> hetNormalAllelicCountSites = new HashSet<>(hetNormalAllelicCounts.getIntervals());
            hetAllelicCounts = new AllelicCountCollection(
                    metadata,
                    filteredAllelicCounts.getRecords().stream()
                            .filter(ac -> hetNormalAllelicCountSites.contains(ac.getInterval()))
                            .collect(Collectors.toList()));
            final File hetAllelicCountsFile = new File(outputDir, outputPrefix + HET_ALLELIC_COUNTS_FILE_SUFFIX);
            hetAllelicCounts.write(hetAllelicCountsFile);
            logger.info(String.format("Allelic counts for case sample at heterozygous sites in matched normal written to %s.", hetAllelicCountsFile));
        }
    }

    private static double calculateHomozygousLogRatio(final AllelicCount allelicCount,
                                                      final double genotypingBaseErrorRate) {
        final int r = allelicCount.getRefReadCount();
        final int n = allelicCount.getTotalReadCount();
        final double betaAll = Beta.regularizedBeta(1, r + 1, n - r + 1);
        final double betaError = Beta.regularizedBeta(genotypingBaseErrorRate, r + 1, n - r + 1);
        final double betaOneMinusError = Beta.regularizedBeta(1 - genotypingBaseErrorRate, r + 1, n - r + 1);
        final double betaHom = betaError + betaAll - betaOneMinusError;
        final double betaHet = betaOneMinusError - betaError;
        return FastMath.log(betaHom) - FastMath.log(betaHet);
    }

    private AlleleFractionSegmentCollection performAlleleFractionSegmentation() {
        logger.info("Starting segmentation of heterozygous allelic counts...");
        final int maxNumChangepointsPerChromosome = maxNumSegmentsPerChromosome - 1;
        return new AlleleFractionKernelSegmenter(hetAllelicCounts)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceAlleleFraction, kernelApproximationDimension,
                        ImmutableSet.copyOf(windowSizes).asList(),
                        numChangepointsPenaltyFactor, numChangepointsPenaltyFactor);
    }

    private void writeModeledSegmentsAndParameterFiles(final MultidimensionalModeller modeller,
                                                       final String fileTag) {
        final ModeledSegmentCollection modeledSegments = modeller.getModeledSegments();
        writeSegments(modeledSegments, fileTag + SEGMENTS_FILE_SUFFIX);
        final File copyRatioParameterFile = new File(outputDir, outputPrefix + fileTag + COPY_RATIO_MODEL_PARAMETER_FILE_SUFFIX);
        final File alleleFractionParameterFile = new File(outputDir, outputPrefix + fileTag + ALLELE_FRACTION_MODEL_PARAMETER_FILE_SUFFIX);
        modeller.writeModelParameterFiles(copyRatioParameterFile, alleleFractionParameterFile);
    }

    private void writeSegments(final AbstractSampleLocatableCollection<?> segments,
                               final String fileSuffix) {
        final File segmentsFile = new File(outputDir, outputPrefix + fileSuffix);
        segments.write(segmentsFile);
        logger.info(String.format("Segments written to %s", segmentsFile));
    }
}