package org.broadinstitute.hellbender.tools.copynumber.denoising;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.CreateReadCountPanelOfNormals;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.MatrixSummaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.IntStream;

/**
 * Utility class for package-private methods for performing SVD-based denoising and related operations.
 *
 * These methods are specifically tailored for the SVD-denoising methods used in the GATK CNV pipeline
 * and are not intended for wide reuse.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class SVDDenoisingUtils {
    private static final Logger logger = LogManager.getLogger(SVDDenoisingUtils.class);

    private static final double EPSILON = 1E-9;
    private static final double INV_LN2 = MathUtils.INV_LOG_2;
    private static final double LN2_EPSILON = Math.log(EPSILON) * INV_LN2;

    private SVDDenoisingUtils() {}

    static final class PreprocessedStandardizedResult {
        final RealMatrix preprocessedStandardizedValues;
        final double[] panelIntervalFractionalMedians;
        final boolean[] filterSamples;
        final boolean[] filterIntervals;

        private PreprocessedStandardizedResult(final RealMatrix preprocessedStandardizedValues,
                                               final double[] panelIntervalFractionalMedians,
                                               final boolean[] filterSamples,
                                               final boolean[] filterIntervals) {
            this.preprocessedStandardizedValues = preprocessedStandardizedValues;
            this.panelIntervalFractionalMedians = panelIntervalFractionalMedians;
            this.filterSamples = filterSamples;
            this.filterIntervals = filterIntervals;
        }
    }

    /**
     * Preprocess (i.e., transform to fractional coverage, correct GC bias, filter, impute, and truncate)
     * and standardize read counts from a panel of normals.
     * All inputs are assumed to be valid.
     * The dimensions of {@code readCounts} should be samples x intervals.
     * To reduce memory footprint, {@code readCounts} is modified in place when possible.
     * Filtering is performed by using boolean arrays to keep track of intervals and samples
     * that have been filtered at any step and masking {@code readCounts} with them appropriately.
     * If {@code intervalGCContent} is null, GC-bias correction will not be performed.
     */
    static PreprocessedStandardizedResult preprocessAndStandardizePanel(final RealMatrix readCounts,
                                                                        final double[] intervalGCContent,
                                                                        final double minimumIntervalMedianPercentile,
                                                                        final double maximumZerosInSamplePercentage,
                                                                        final double maximumZerosInIntervalPercentage,
                                                                        final double extremeSampleMedianPercentile,
                                                                        final boolean doImputeZeros,
                                                                        final double extremeOutlierTruncationPercentile) {
        //preprocess (transform to fractional coverage, correct GC bias, filter, impute, truncate) and return copy of submatrix
        logger.info("Preprocessing read counts...");
        final PreprocessedStandardizedResult preprocessedStandardizedResult = preprocessPanel(readCounts, intervalGCContent,
                minimumIntervalMedianPercentile, maximumZerosInSamplePercentage, maximumZerosInIntervalPercentage,
                extremeSampleMedianPercentile, doImputeZeros, extremeOutlierTruncationPercentile);
        logger.info("Panel read counts preprocessed.");

        //standardize in place
        logger.info("Standardizing read counts...");
        divideBySampleMedianAndTransformToLog2(preprocessedStandardizedResult.preprocessedStandardizedValues);
        logger.info("Subtracting median of sample medians...");
        final double[] sampleLog2Medians = MatrixSummaryUtils.getRowMedians(preprocessedStandardizedResult.preprocessedStandardizedValues);
        final double medianOfSampleMedians = new Median().evaluate(sampleLog2Medians);
        preprocessedStandardizedResult.preprocessedStandardizedValues
                .walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int sampleIndex, int intervalIndex, double value) {
                return value - medianOfSampleMedians;
            }
        });
        logger.info("Panel read counts standardized.");

        return preprocessedStandardizedResult;
    }

    /**
     * Perform SVD-based denoising of integer read counts for a single sample using a panel of normals.
     * Only the eigensamples (which are sorted by singular value in decreasing order) specified by
     * {@code numEigensamples} are used to denoise.
     */
    static SVDDenoisedCopyRatioResult denoise(final SVDReadCountPanelOfNormals panelOfNormals,
                                              final SimpleCountCollection readCounts,
                                              final int numEigensamples) {
        Utils.nonNull(panelOfNormals);
        if (!CopyNumberArgumentValidationUtils.isSameDictionary(panelOfNormals.getSequenceDictionary(), readCounts.getMetadata().getSequenceDictionary())) {
            logger.warn("Sequence dictionaries in panel and case sample do not match.");
        }
        ParamUtils.isPositiveOrZero(numEigensamples, "Number of eigensamples to use for denoising must be non-negative.");
        Utils.validateArg(numEigensamples <= panelOfNormals.getNumEigensamples(),
                "Number of eigensamples to use for denoising is greater than the number available in the panel of normals.");

        logger.info("Validating sample intervals against original intervals used to build panel of normals...");
        Utils.validateArg(panelOfNormals.getOriginalIntervals().equals(readCounts.getIntervals()),
                "Sample intervals must be identical to the original intervals used to build the panel of normals.");

        logger.info("Preprocessing and standardizing sample read counts...");
        final RealMatrix standardizedCopyRatioValues = preprocessAndStandardizeSample(panelOfNormals, readCounts.getCounts());

        final RealMatrix denoisedCopyRatioValues;
        if (numEigensamples == 0 || panelOfNormals.getNumEigensamples() == 0) {
            logger.warn("A zero number of eigensamples was specified or no eigensamples were available to perform denoising; " +
                    "denoised copy ratios will be identical to the standardized copy ratios...");
            denoisedCopyRatioValues = standardizedCopyRatioValues;
        } else {
            logger.info(String.format("Using %d out of %d eigensamples to denoise...", numEigensamples, panelOfNormals.getNumEigensamples()));
            logger.info("Subtracting projection onto space spanned by eigensamples...");
            denoisedCopyRatioValues = subtractProjection(standardizedCopyRatioValues, panelOfNormals.getEigensampleVectors(), numEigensamples);
        }

        logger.info("Sample denoised.");

        //construct the result
        return new SVDDenoisedCopyRatioResult(
                readCounts.getMetadata(),
                panelOfNormals.getPanelIntervals(),
                standardizedCopyRatioValues,
                denoisedCopyRatioValues);
    }

    /**
     * Preprocess (i.e., transform to fractional coverage and correct GC bias)
     * and standardize read counts for samples when no panel of normals is available.
     * The original {@code readCounts} has dimensions 1 x intervals and is not modified.
     * If {@code intervalGCContent} is null, GC-bias correction will not be performed
     */
    public static RealMatrix preprocessAndStandardizeSample(final double[] readCounts,
                                                            final double[] intervalGCContent) {
        Utils.nonNull(readCounts);
        Utils.validateArg(intervalGCContent == null || readCounts.length == intervalGCContent.length,
                "Number of intervals for read counts must match those for GC-content annotations.");

        RealMatrix result = new Array2DRowRealMatrix(new double[][]{readCounts});

        //preprocess (transform to fractional coverage, correct GC bias) copy in place
        logger.info("Preprocessing read counts...");
        transformToFractionalCoverage(result);
        performOptionalGCBiasCorrection(result, intervalGCContent);
        logger.info("Sample read counts preprocessed.");

        //standardize copy in place
        logger.info("Standardizing read counts...");
        divideBySampleMedianAndTransformToLog2(result);
        logger.info("Subtracting sample median...");
        final double[] sampleLog2Medians = MatrixSummaryUtils.getRowMedians(result);
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int sampleIndex, int intervalIndex, double value) {
                return value - sampleLog2Medians[sampleIndex];
            }
        });
        logger.info("Sample read counts standardized.");

        return result;
    }

    /**
     * This method is purposely written as a single contiguous chunk of code rather than broken up into sub-methods.
     * This is to make efficient reuse of intermediate results without requiring them to be passed by reference.
     * Please do not refactor or extract code without a very good reason.
     */
    private static PreprocessedStandardizedResult preprocessPanel(final RealMatrix readCounts,
                                                                  final double[] intervalGCContent,
                                                                  final double minimumIntervalMedianPercentile,
                                                                  final double maximumZerosInSamplePercentage,
                                                                  final double maximumZerosInIntervalPercentage,
                                                                  final double extremeSampleMedianPercentile,
                                                                  final boolean doImputeZeros,
                                                                  final double extremeOutlierTruncationPercentile) {
        transformToFractionalCoverage(readCounts);
        performOptionalGCBiasCorrection(readCounts, intervalGCContent);

        final int numOriginalSamples = readCounts.getRowDimension();
        final int numOriginalIntervals = readCounts.getColumnDimension();

        final boolean[] filterSamples = new boolean[numOriginalSamples];
        final boolean[] filterIntervals = new boolean[numOriginalIntervals];

        //filter intervals by fractional median
        final double[] originalIntervalMedians = MatrixSummaryUtils.getColumnMedians(readCounts);
        if (minimumIntervalMedianPercentile == 0.) {
            logger.info(String.format("A value of 0 was provided for argument %s, so the corresponding filtering step will be skipped...",
                    CreateReadCountPanelOfNormals.MINIMUM_INTERVAL_MEDIAN_PERCENTILE_LONG_NAME));
        } else {
            logger.info(String.format("Filtering intervals with median (across samples) less than or equal to the %.2f percentile...", minimumIntervalMedianPercentile));
            //calculate percentile
            final double minimumIntervalMedianThreshold = new Percentile(minimumIntervalMedianPercentile).evaluate(originalIntervalMedians);
            //filter intervals
            IntStream.range(0, numOriginalIntervals)
                    .filter(intervalIndex -> originalIntervalMedians[intervalIndex] <= minimumIntervalMedianThreshold)
                    .forEach(intervalIndex -> filterIntervals[intervalIndex] = true);
            logger.info(String.format("After filtering, %d out of %d intervals remain...", countNumberPassingFilter(filterIntervals), numOriginalIntervals));
        }

        logger.info("Dividing by interval medians...");
        IntStream.range(0, numOriginalIntervals)
                .filter(intervalIndex -> !filterIntervals[intervalIndex])
                .forEach(intervalIndex -> IntStream.range(0, numOriginalSamples).filter(sampleIndex -> !filterSamples[sampleIndex]).forEach(sampleIndex -> {
                    final double value = readCounts.getEntry(sampleIndex, intervalIndex);
                    readCounts.setEntry(sampleIndex, intervalIndex,value / originalIntervalMedians[intervalIndex]);
                }));

        //filter samples by percentage of zero-coverage intervals not already filtered
        if (maximumZerosInSamplePercentage == 100.) {
            logger.info(String.format("A value of 100 was provided for argument %s, so the corresponding filtering step will be skipped...",
                    CreateReadCountPanelOfNormals.MAXIMUM_ZEROS_IN_SAMPLE_PERCENTAGE_LONG_NAME));
        } else {
            logger.info(String.format("Filtering samples with a fraction of zero-coverage intervals above %.2f percent...", maximumZerosInSamplePercentage));
            final int maxZerosInSample = calculateMaximumZerosCount(countNumberPassingFilter(filterIntervals), maximumZerosInSamplePercentage);
            IntStream.range(0, numOriginalSamples)
                    .filter(sampleIndex -> !filterSamples[sampleIndex])
                    .forEach(sampleIndex -> {
                        final int numZerosInSample = (int) IntStream.range(0, numOriginalIntervals)
                                .filter(intervalIndex -> !filterIntervals[intervalIndex] && readCounts.getEntry(sampleIndex, intervalIndex) == 0.)
                                .count();
                        if (numZerosInSample > maxZerosInSample) {
                            filterSamples[sampleIndex] = true;
                        }
                    });
            logger.info(String.format("After filtering, %d out of %d samples remain...", countNumberPassingFilter(filterSamples), numOriginalSamples));
        }

        //filter intervals by percentage of zero-coverage samples not already filtered
        if (maximumZerosInIntervalPercentage == 100.) {
            logger.info(String.format("A value of 100 was provided for argument %s, so the corresponding filtering step will be skipped...",
                    CreateReadCountPanelOfNormals.MAXIMUM_ZEROS_IN_INTERVAL_PERCENTAGE_LONG_NAME));
        } else {
            logger.info(String.format("Filtering intervals with a fraction of zero-coverage samples above %.2f percent...", maximumZerosInIntervalPercentage));
            final int maxZerosInInterval = calculateMaximumZerosCount(countNumberPassingFilter(filterSamples), maximumZerosInIntervalPercentage);
            IntStream.range(0, numOriginalIntervals)
                    .filter(intervalIndex -> !filterIntervals[intervalIndex])
                    .forEach(intervalIndex -> {
                        final int numZerosInInterval = (int) IntStream.range(0, numOriginalSamples)
                                .filter(sampleIndex -> !filterSamples[sampleIndex] && readCounts.getEntry(sampleIndex, intervalIndex) == 0.)
                                .count();
                        if (numZerosInInterval > maxZerosInInterval) {
                            filterIntervals[intervalIndex] = true;
                        }
                    });
            logger.info(String.format("After filtering, %d out of %d intervals remain...", countNumberPassingFilter(filterIntervals), numOriginalIntervals));
        }

        //filter samples with extreme medians
        if (extremeSampleMedianPercentile == 0.) {
            logger.info(String.format("A value of 0 was provided for argument %s, so the corresponding filtering step will be skipped...",
                    CreateReadCountPanelOfNormals.EXTREME_SAMPLE_MEDIAN_PERCENTILE_LONG_NAME));
        } else {
            logger.info(String.format("Filtering samples with a median (across intervals) below the %.2f percentile or above the %.2f percentile...",
                    extremeSampleMedianPercentile, 100. - extremeSampleMedianPercentile));
            //calculate the medians for all samples (which, although unnecessary, makes bookkeeping easier) across intervals not already filtered
            final double[] sampleMedians = IntStream.range(0, numOriginalSamples)
                    .mapToDouble(sampleIndex -> new Median().evaluate(IntStream.range(0, numOriginalIntervals)
                            .filter(intervalIndex -> !filterIntervals[intervalIndex])
                            .mapToDouble(intervalIndex -> readCounts.getEntry(sampleIndex, intervalIndex))
                            .toArray()))
                    .toArray();
            //calculate percentiles
            final double minimumSampleMedianThreshold = new Percentile(extremeSampleMedianPercentile).evaluate(sampleMedians);
            final double maximumSampleMedianThreshold = new Percentile(100. - extremeSampleMedianPercentile).evaluate(sampleMedians);
            //filter samples
            IntStream.range(0, numOriginalSamples)
                    .filter(sampleIndex -> sampleMedians[sampleIndex] < minimumSampleMedianThreshold || sampleMedians[sampleIndex] > maximumSampleMedianThreshold)
                    .forEach(sampleIndex -> filterSamples[sampleIndex] = true);
            logger.info(String.format("After filtering, %d out of %d samples remain...", countNumberPassingFilter(filterSamples), numOriginalSamples));
        }

        //construct the filtered results as a new matrix, which will be modified in place from this point on
        final int[] panelIntervalIndices = IntStream.range(0, numOriginalIntervals).filter(intervalIndex -> !filterIntervals[intervalIndex]).toArray();
        final int[] panelSampleIndices = IntStream.range(0, numOriginalSamples).filter(sampleIndex -> !filterSamples[sampleIndex]).toArray();
        final RealMatrix preprocessedReadCounts = readCounts.getSubMatrix(panelSampleIndices, panelIntervalIndices);
        final double[] panelIntervalFractionalMedians = IntStream.range(0, numOriginalIntervals)
                .filter(intervalIndex -> !filterIntervals[intervalIndex])
                .mapToDouble(intervalIndex -> originalIntervalMedians[intervalIndex]).toArray();

        //garbage collection to clean up readCounts
        logHeapUsage();
        logger.info("Performing garbage collection...");
        System.gc();
        logHeapUsage();

        //impute zeros as median of non-zero values in interval
        if (!doImputeZeros) {
            logger.info("Skipping imputation of zero-coverage values...");
        } else {
            final int numPanelIntervals = panelIntervalIndices.length;
            final double[] intervalNonZeroMedians = IntStream.range(0, numPanelIntervals)
                    .mapToObj(intervalIndex -> Arrays.stream(preprocessedReadCounts.getColumn(intervalIndex)).filter(value -> value > 0.).toArray())
                    .mapToDouble(nonZeroValues -> new Median().evaluate(nonZeroValues))
                    .toArray();
            final int[] numImputed = {0};  //needs to be effectively final to be used inside visitor
            preprocessedReadCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                @Override
                public double visit(int sampleIndex, int intervalIndex, double value) {
                    if (value == 0.) {
                        numImputed[0]++;
                        return intervalNonZeroMedians[intervalIndex];
                    }
                    return value;
                }
            });
            logger.info(String.format("%d zero-coverage values were imputed to the median of the non-zero values in the corresponding interval...",
                    numImputed[0]));
        }

        //truncate extreme values to the corresponding percentile
        if (extremeOutlierTruncationPercentile == 0.) {
            logger.info(String.format("A value of 0 was provided for argument %s, so the corresponding truncation step will be skipped...",
                    CreateReadCountPanelOfNormals.EXTREME_OUTLIER_TRUNCATION_PERCENTILE_LONG_NAME));
        } else if ((long) preprocessedReadCounts.getRowDimension() * preprocessedReadCounts.getColumnDimension() > Integer.MAX_VALUE) {
            logger.warn("The number of matrix elements exceeds Integer.MAX_VALUE, so outlier truncation will be skipped...");
        } else {
            final double[] values = Doubles.concat(preprocessedReadCounts.getData());
            final double minimumOutlierTruncationThreshold = new Percentile(extremeOutlierTruncationPercentile).evaluate(values);
            final double maximumOutlierTruncationThreshold = new Percentile(100. - extremeOutlierTruncationPercentile).evaluate(values);
            final int[] numTruncated = {0};  //needs to be effectively final to be used inside visitor
            preprocessedReadCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
                @Override
                public double visit(int sampleIndex, int intervalIndex, double value) {
                    if (value < minimumOutlierTruncationThreshold) {
                        numTruncated[0]++;
                        return minimumOutlierTruncationThreshold;
                    }
                    if (value > maximumOutlierTruncationThreshold) {
                        numTruncated[0]++;
                        return maximumOutlierTruncationThreshold;
                    }
                    return value;
                }
            });
            logger.info(String.format("%d values below the %.2f percentile or above the %.2f percentile were truncated to the corresponding value...",
                    numTruncated[0], extremeOutlierTruncationPercentile, 100. - extremeOutlierTruncationPercentile));
        }
        return new PreprocessedStandardizedResult(
                preprocessedReadCounts, panelIntervalFractionalMedians, filterSamples, filterIntervals);
    }

    private static void logHeapUsage() {
        final int mb = 1024 * 1024;
        final Runtime runtime = Runtime.getRuntime();
        logger.info("Heap utilization statistics [MB]:");
        logger.info("Used memory: " + (runtime.totalMemory() - runtime.freeMemory()) / mb);
        logger.info("Free memory: " + runtime.freeMemory() / mb);
        logger.info("Total memory: " + runtime.totalMemory() / mb);
        logger.info("Maximum memory: " + runtime.maxMemory() / mb);
    }

    /**
     * Preprocess (i.e., transform to fractional coverage, correct GC bias, subset, divide by fractional medians)
     * and standardize read counts for samples, using interval fractional medians from a panel of normals.
     * The original {@code readCounts} has dimensions 1 x intervals and is not modified.
     */
    private static RealMatrix preprocessAndStandardizeSample(final SVDReadCountPanelOfNormals panelOfNormals,
                                                             final double[] readCounts) {
        RealMatrix result = new Array2DRowRealMatrix(new double[][]{readCounts});

        //preprocess (transform to fractional coverage, correct GC bias, subset, divide by fractional medians) copy in place
        logger.info("Preprocessing read counts...");
        transformToFractionalCoverage(result);
        performOptionalGCBiasCorrection(result, panelOfNormals.getOriginalIntervalGCContent());

        logger.info("Subsetting sample intervals to post-filter panel intervals...");
        final Set<SimpleInterval> panelIntervals = new HashSet<>(panelOfNormals.getPanelIntervals());
        final int[] subsetIntervalIndices = IntStream.range(0, panelOfNormals.getOriginalIntervals().size())
                .filter(i -> panelIntervals.contains(panelOfNormals.getOriginalIntervals().get(i)))
                .toArray();
        result = result.getSubMatrix(new int[]{0}, subsetIntervalIndices);

        logger.info("Dividing by interval medians from the panel of normals...");
        final double[] intervalMedians = panelOfNormals.getPanelIntervalFractionalMedians();
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int sampleIndex, int intervalIndex, double value) {
                return value / intervalMedians[intervalIndex];
            }
        });
        logger.info("Sample read counts preprocessed.");

        //standardize copy in place
        logger.info("Standardizing read counts...");
        divideBySampleMedianAndTransformToLog2(result);
        logger.info("Subtracting sample median...");
        final double[] sampleLog2Medians = MatrixSummaryUtils.getRowMedians(result);
        result.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int sampleIndex, int intervalIndex, double value) {
                return value - sampleLog2Medians[sampleIndex];
            }
        });
        logger.info("Sample read counts standardized.");

        return result;
    }

    /**
     * Given standardized read counts specified by a row vector S (dimensions {@code 1 x M})
     * and all eigensample vectors U (dimensions {@code M x K}),
     * returns s - s U<sub>k</sub> U<sub>k</sub><sup>T</sup>,
     * where U<sub>k</sub> contains the first {@code numEigensamples}.
     */
    private static RealMatrix subtractProjection(final RealMatrix standardizedValues,
                                                 final double[][] eigensampleVectors,
                                                 final int numEigensamples) {
        if (numEigensamples == 0) {
            return standardizedValues.copy();
        }

        final int numIntervals = eigensampleVectors.length;
        final int numAllEigensamples = eigensampleVectors[0].length;

        logger.info("Distributing the standardized read counts...");

        logger.info("Composing eigensample matrix for the requested number of eigensamples and transposing them...");
        final RealMatrix eigensampleTruncatedMatrix = numEigensamples == numAllEigensamples
                ? new Array2DRowRealMatrix(eigensampleVectors, false)
                : new Array2DRowRealMatrix(eigensampleVectors, false).getSubMatrix(0, numIntervals - 1, 0, numEigensamples - 1);

        logger.info("Computing projection...");
        final RealMatrix projection = standardizedValues
                .multiply(eigensampleTruncatedMatrix)
                .multiply(eigensampleTruncatedMatrix.transpose());

        logger.info("Subtracting projection...");
        return standardizedValues.subtract(projection);
    }

    private static int countNumberPassingFilter(final boolean[] filter) {
        final int numPassingFilter = (int) IntStream.range(0, filter.length).filter(i -> !filter[i]).count();
        if (numPassingFilter == 0) {
            throw new UserException.BadInput("Filtering removed all samples or intervals.  Select less strict filtering criteria.");
        }
        return numPassingFilter;
    }

    private static void transformToFractionalCoverage(final RealMatrix matrix) {
        logger.info("Transforming read counts to fractional coverage...");
        final double[] sampleSums = MathUtils.rowSums(matrix);
        matrix.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int sampleIndex, int intervalIndex, double value) {
                return value / sampleSums[sampleIndex];
            }
        });
    }

    private static void performOptionalGCBiasCorrection(final RealMatrix matrix, 
                                                        final double[] intervalGCContent) {
        if (intervalGCContent != null) {
            logger.info("Performing GC-bias correction...");
            GCBiasCorrector.correctGCBias(matrix, intervalGCContent);
        }
    }

    private static void divideBySampleMedianAndTransformToLog2(final RealMatrix matrix) {
        logger.info("Dividing by sample medians and transforming to log2 space...");
        final double[] sampleMedians = MatrixSummaryUtils.getRowMedians(matrix);
        IntStream.range(0, sampleMedians.length).forEach(sampleIndex ->
                ParamUtils.isPositive(sampleMedians[sampleIndex],
                        sampleMedians.length == 1
                                ? "Sample does not have a non-negative sample median."
                                : String.format("Sample at index %s does not have a non-negative sample median.", sampleIndex)));
        matrix.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int sampleIndex, int intervalIndex, double value) {
                return safeLog2(value / sampleMedians[sampleIndex]);
            }
        });
    }

    private static int calculateMaximumZerosCount(final int numTotalCounts,
                                                  final double percentage) {
        return (int) Math.ceil(numTotalCounts * percentage / 100.0);
    }

    private static double safeLog2(final double x) {
        return x < EPSILON ? LN2_EPSILON : Math.log(x) * INV_LN2;
    }
}