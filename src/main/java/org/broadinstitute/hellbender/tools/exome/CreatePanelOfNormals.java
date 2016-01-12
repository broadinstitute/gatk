package org.broadinstitute.hellbender.tools.exome;

import com.google.api.client.repackaged.com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SparkToggleCommandLineProgram;
import org.broadinstitute.hellbender.utils.hdf5.HDF5File;
import org.broadinstitute.hellbender.utils.hdf5.HDF5PoN;
import org.broadinstitute.hellbender.utils.svd.SVD;
import org.broadinstitute.hellbender.utils.svd.SVDFactory;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

/**
 * Tool to create a panel of normals (PoN) given a collection of read-counts
 * for control samples.
 *
 * <p>
 * The input read-counts consists of a single file with counts for several
 * samples. This might be constructed from several single sample read counts using
 * the {@link CombineReadCounts} tool.
 * </p>
 * <p>
 * Accordingly the input format is exactly the same as the output format for
 * that tool, which is described in its documentation {@link CombineReadCounts here}.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Creates a Panel of Normals (PoN) given the proportional read counts for the samples that are part of the panel.  Supports Apache Spark for some operations.",
        oneLineSummary = "Creates a Panel of Normals.",
        programGroup = CopyNumberProgramGroup.class
)
public class CreatePanelOfNormals extends SparkToggleCommandLineProgram {

    static final long serialVersionUID = 42123132L;

    private static final double INV_LN_2 = 1.0 / Math.log(2);

    public static final double EPSILON = 10e-10;

    public static final double DEFAULT_TARGET_FACTOR_THRESHOLD_PERCENTILE = 25.0;

    public static final double DEFAULT_COLUMN_OUTLIER_DROP_THRESHOLD_PERCENTILE = 2.5;

    // 0.1% by default as copied from python code... perhaps a bug there as most of these constants in that code
    // are then multiplied by 100 or subtracted from a 100, check out PanelCleaner.py in ReCapSeg repo.
    public static final double DEFAULT_OUTLIER_TRUNCATE_PERCENTILE_THRESHOLD = 0.1;

    public static final double DEFAULT_MAXIMUM_PERCENT_ZEROS_IN_COLUMN = 2.0;

    public static final double DEFAULT_MAXIMUM_PERCENT_ZEROS_IN_TARGET = 5.0;

    public static final String INFER_NUMBER_OF_EIGEN_SAMPLES = "auto";

    public static final String DEFAULT_NUMBER_OF_EIGEN_SAMPLES = INFER_NUMBER_OF_EIGEN_SAMPLES;

    public static final String NUMBER_OF_EIGEN_SAMPLES_DOCUMENTATION =
            "Number of eigen samples to use for the reduced PoN. " +
            "By default it will infer the appropriate number of eigen samples (value " + INFER_NUMBER_OF_EIGEN_SAMPLES + ")";

    public static final String COLUMN_EXTREME_THRESHOLD_PERCENTILE_DOCUMENTATION =
            "Percentile for the two-tailed extreme median column coverage filter. " +
            "Columns that have a median count in the that bottom or top percentile are excluded from the PoN. " +
            "Any value in the range (0, 50) (default is " + DEFAULT_COLUMN_OUTLIER_DROP_THRESHOLD_PERCENTILE + ")";

    public static final String TARGET_FACTOR_THRESHOLD_PERCENTILE_DOCUMENTATION =
            "Percentile to determine the minimum target factor for any target to be considered part of the PoN. " +
            "Any value in the range (0, 100) (default is " + DEFAULT_TARGET_FACTOR_THRESHOLD_PERCENTILE + ")";

    public static final String MAXIMUM_PERCENT_ZEROS_IN_COLUMN_DOCUMENTATION =
            "Maximum percentage of 0 counts in a column to drop it from the panel. " +
            "Any value in the range (0, 100) (default is " + DEFAULT_MAXIMUM_PERCENT_ZEROS_IN_COLUMN + ")";

    public static final String MAXIMUM_PERCENT_ZEROS_IN_TARGET_DOCUMENTATION =
            "Maximum percentage of 0 counts for a target to drop it from the panel. " +
                    "Any value in the range (0, 100) (default is " + DEFAULT_MAXIMUM_PERCENT_ZEROS_IN_TARGET + ")";

    public static final String COUNT_TRUNCATE_PERCENTILE_DOCUMENTATION =
            "Percentiles to obtain the maximum and minimum value for any count in the panel. " +
            "Any value outside the resulting range would be set will be truncated to that minimum or maximum. " +
            "Valid values are within the range (0, 50) (default is " + DEFAULT_OUTLIER_TRUNCATE_PERCENTILE_THRESHOLD + ")";


    public static final String TARGET_FACTOR_THRESHOLD_PERCENTILE_SHORT_NAME = "minTFPcTh";
    public static final String TARGET_FACTOR_THRESHOLD_PERCENTILE_FULL_NAME = "minimumTargetFactorPercentileThreshold";
    public static final String MAXIMUM_PERCENT_ZEROS_IN_COLUMN_SHORT_NAME = "maxCol0sPc";
    public static final String MAXIMUM_PERCENT_ZEROS_IN_COLUMN_FULL_NAME = "maximumColumnZerosPercentage";
    public static final String MAXIMUM_PERCENT_ZEROS_IN_TARGET_SHORT_NAME = "maxTrg0sPc";
    public static final String MAXIMUM_PERCENT_ZEROS_IN_TARGET_FULL_NAME = "maximumTargetZerosPercentage";
    public static final String COLUMN_EXTREME_THRESHOLD_PERCENTILE_SHORT_NAME = "extremeColMedPrTh";
    public static final String COLUMN_EXTREME_THRESHOLD_PERCENTILE_FULL_NAME = "extremeColumnMedianCountPercentileThreshold";
    public static final String COUNT_TRUNCATE_PERCENTILE_SHORT_NAME = "truncPcThr";
    public static final String COUNT_TRUNCATE_PERCENTILE_FULL_NAME = "truncatePercentileThreshold";
    public static final String NUMBER_OF_EIGEN_SAMPLES_SHORT_NAME = "numEigen";
    public static final String NUMBER_OF_EIGEN_SAMPLES_FULL_NAME = "numberOfEigenSamples";
    public static final String DRY_RUN_SHORT_NAME = "dryRun";
    public static final String DRY_RUN_FULL_NAME = DRY_RUN_SHORT_NAME;

    public static final double JOLLIFES_RULE_MEAN_FACTOR = 0.7;

    @Argument(
            doc = TARGET_FACTOR_THRESHOLD_PERCENTILE_DOCUMENTATION,
            shortName = TARGET_FACTOR_THRESHOLD_PERCENTILE_SHORT_NAME,
            fullName  = TARGET_FACTOR_THRESHOLD_PERCENTILE_FULL_NAME,
            optional  = true
    )
    protected double targetFactorThreshold = DEFAULT_TARGET_FACTOR_THRESHOLD_PERCENTILE;

    @Argument(
            doc = MAXIMUM_PERCENT_ZEROS_IN_COLUMN_DOCUMENTATION,
            shortName = MAXIMUM_PERCENT_ZEROS_IN_COLUMN_SHORT_NAME,
            fullName = MAXIMUM_PERCENT_ZEROS_IN_COLUMN_FULL_NAME,
            optional = true
    )
    protected double maximumPercentZerosInColumn = DEFAULT_MAXIMUM_PERCENT_ZEROS_IN_COLUMN;

    @Argument(
            doc = MAXIMUM_PERCENT_ZEROS_IN_TARGET_DOCUMENTATION,
            shortName = MAXIMUM_PERCENT_ZEROS_IN_TARGET_SHORT_NAME,
            fullName = MAXIMUM_PERCENT_ZEROS_IN_TARGET_FULL_NAME,
            optional = true
    )
    protected double maximumPercentZerosInTarget = DEFAULT_MAXIMUM_PERCENT_ZEROS_IN_TARGET;

    @Argument(
            doc = COLUMN_EXTREME_THRESHOLD_PERCENTILE_DOCUMENTATION,
            shortName = COLUMN_EXTREME_THRESHOLD_PERCENTILE_SHORT_NAME,
            fullName = COLUMN_EXTREME_THRESHOLD_PERCENTILE_FULL_NAME,
            optional = false
    )
    protected double columnExtremeThresholdPercentile = DEFAULT_COLUMN_OUTLIER_DROP_THRESHOLD_PERCENTILE;

    @Argument(
            doc = COUNT_TRUNCATE_PERCENTILE_DOCUMENTATION,
            shortName = COUNT_TRUNCATE_PERCENTILE_SHORT_NAME,
            fullName = COUNT_TRUNCATE_PERCENTILE_FULL_NAME,
            optional = false
    )
    protected double outlierTruncatePercentileThresh = DEFAULT_OUTLIER_TRUNCATE_PERCENTILE_THRESHOLD;

    @Argument(
            doc = NUMBER_OF_EIGEN_SAMPLES_DOCUMENTATION,
            shortName = NUMBER_OF_EIGEN_SAMPLES_SHORT_NAME,
            fullName = NUMBER_OF_EIGEN_SAMPLES_FULL_NAME,
            optional = true
    )
    protected String numberOfEigenSamples = DEFAULT_NUMBER_OF_EIGEN_SAMPLES;

    @Argument(
            doc = "Input proportional read counts for samples in the panel of normals.",
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            fullName  = StandardArgumentDefinitions.INPUT_LONG_NAME,
            optional  = false
    )
    protected File inputFile = null;

    @ArgumentCollection
    protected TargetArgumentCollection targetArguments = new TargetArgumentCollection(() -> inputFile);

    @Argument(
            doc = "Output HDF5 file name.",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            optional = false
    )
    protected File outFile = null;


    // This option is useful to test performance when the HDF5 lib is not
    // present for whatever reason.
    @Argument(
            doc = "Dry-run, skip the creation of the HDF5 output file.",
            shortName = DRY_RUN_SHORT_NAME,
            fullName  = DRY_RUN_FULL_NAME,
            optional  = true
    )
    protected boolean dryRun = false;


    protected Object createPoN(final JavaSparkContext ctx) {

        final OptionalInt numberOfEigenSamples = calculatePreferredNumberOfEigenSamples();

        final double targetFactorPercentileThreshold = calculateTargetFactorsPercentileThreshold();
        final double extremeColumnMedianCountPercentileThreshold = calculateExtremeColumnMedianCountsPercentileThreshold();
        final double countTruncatePercentile = checkCountTruncatePercentile();

        // Remove low coverage targets:
        final TargetCollection<Target> targets = targetArguments.readTargetCollection(true);
        final Pair<ReadCountCollection, double[]> inputSubsetByUsableTargets = subsetReadCountsToUsableTargets(readReadCountsFromFile(inputFile, targets), targetFactorPercentileThreshold, logger);
        ReadCountCollection readCounts = inputSubsetByUsableTargets.getLeft();
        final double[] targetFactors = inputSubsetByUsableTargets.getRight();

        // Normalize read-counts by removing the target factor component:
        normalizeReadCountsByTargetFactors(readCounts, targetFactors);

        // Create and write the first part of the PoN HDF5 file.
        if (!dryRun) {
            writeTargetFactorNormalizeReadCountsAndTargetFactors(outFile, readCounts, targetFactors);
        }

        // Remove column and targets with too many zeros and targets with
        // extreme median coverage:
        final int maximumColumnZeros = calculateMaximumColumnZerosPercentage(readCounts.targets().size());
        final int maximumTargetZeros = calculateMaximumTargetZerosPercentage(readCounts.columnNames().size());
        readCounts = removeColumnsWithTooManyZeros(readCounts, maximumColumnZeros, logger);
        readCounts = removeTargetsWithTooManyZeros(readCounts, maximumTargetZeros, logger);
        readCounts = removeColumnsWithExtremeMedianCounts(readCounts, extremeColumnMedianCountPercentileThreshold, logger);

        // Impute zero counts to be same as the median for the same target.
        // This happens in-place.
        imputeZerosCounts(readCounts, logger);

        // Truncate extreme count values.
        // This happens in-place.
        truncateExtremeCounts(readCounts, countTruncatePercentile, logger);

        // Normalized by the median and log scale the read counts.
        normalizeAndLogReadCounts(readCounts, logger);
        subtractBGSCenter(readCounts, logger);

        final ReductionResult reduction = calculateReducedPanelAndPInverses(readCounts, numberOfEigenSamples, logger, ctx);

        // Write the log-normals and reduced panel matrices.
        if (!dryRun) {
            writeLogNormalsReducedPanel(outFile, readCounts, reduction);
        }
        return "SUCCESS";
    }

    @Override
    protected void runPipeline(JavaSparkContext ctx) {
        createPoN(ctx);
    }

    // Writes normalized and reduced panel matrices to the output HDF5 file.
    private void writeLogNormalsReducedPanel(final File outFile, final ReadCountCollection readCounts, final ReductionResult reduction) {
        try (final HDF5File file = new HDF5File(outFile, HDF5File.OpenMode.READ_WRITE)) {
            final HDF5PoN pon = new HDF5PoN(file);
            pon.setPanelSampleNames(readCounts.columnNames());
            pon.setPanelTargetNames(readCounts.targets().stream()
                    .map(Target::getName)
                    .collect(Collectors.toList()));
            pon.setLogNormalCounts(readCounts.counts());
            pon.setLogNormalPInverseCounts(reduction.pseudoInverse);
            pon.setReducedPanelCounts(reduction.reduced);
            pon.setReducedPanelPInverseCounts(reduction.reducedInverse);
            pon.setVersion(4.0);
        }
    }

    private static void writeTargetFactorNormalizeReadCountsAndTargetFactors(final File outFile, final ReadCountCollection readCounts, final double[] targetFactors) {
        final List<String> sampleNames = readCounts.columnNames();
        final List<String> targetNames = readCounts.targets().stream()
                .map(Target::getName)
                .collect(Collectors.toList());
        try (final HDF5File file = new HDF5File(outFile, HDF5File.OpenMode.CREATE)) {
            final HDF5PoN pon = new HDF5PoN(file);
            pon.setSampleNames(sampleNames);
            pon.setTargetNames(targetNames);
            pon.setTargetFactors(new Array2DRowRealMatrix(targetFactors));
            pon.setNormalCounts(readCounts.counts());
        }
    }

    /**
     * Composes the preferred number of eigen values optional given the user input.
     *
     * @return an empty optional if the user elected to use the automatic/inferred value, otherwise
     * a strictly positive integer.
     */
    private OptionalInt calculatePreferredNumberOfEigenSamples() {
        if (numberOfEigenSamples.equalsIgnoreCase(INFER_NUMBER_OF_EIGEN_SAMPLES)) {
            return OptionalInt.empty();
        } else {
            try {
                final int result = Integer.parseInt(numberOfEigenSamples);
                if (result <= 0) {
                    throw new UserException.BadArgumentValue(NUMBER_OF_EIGEN_SAMPLES_FULL_NAME, "0 or negative values are not allowed: " + numberOfEigenSamples);
                } else {
                    return OptionalInt.of(result);
                }
            } catch (final NumberFormatException ex) {
                throw new UserException.BadArgumentValue(NUMBER_OF_EIGEN_SAMPLES_FULL_NAME,
                        "it must be either '" + INFER_NUMBER_OF_EIGEN_SAMPLES + "' or an integer value");
            }
        }
    }

    /**
     * Structure to hold the result of the SVD/reduction steps.
     */
    @VisibleForTesting
    static final class ReductionResult {
        final RealMatrix pseudoInverse;
        final RealMatrix reduced;
        final RealMatrix reducedInverse;
        final double[] allSingularValues;

        public ReductionResult(final RealMatrix pseudoInverse, final RealMatrix reduced, final RealMatrix reduceInverse, final double[] allSingularValues) {
            this.pseudoInverse = pseudoInverse;
            this.reduced = reduced;
            this.reducedInverse = reduceInverse;
            this.allSingularValues = allSingularValues;
        }
    }

    /**
     * SVD and Pseudo inverse calculation.
     *
     * @param logNormals the input counts for the SVD and reduction steps, fully normalized and already logged.
     * @param requestedNumberOfEigenSamples user requested number of eigen samples for the reduced panel.
     * @return never {@code null}.
     */
    @VisibleForTesting
    static ReductionResult calculateReducedPanelAndPInverses(final ReadCountCollection logNormals, final OptionalInt requestedNumberOfEigenSamples, final Logger logger, final JavaSparkContext ctx) {

        if (ctx == null) {
            logger.warn("No Spark context provided, not going to use Spark...");
        }

        final RealMatrix logNormalCounts = logNormals.counts();
        final int numberOfCountColumns = logNormalCounts.getColumnDimension();

        logger.info("Starting the SVD decomposition of the log-normal counts ...");
        final long svdStartTime = System.currentTimeMillis();
        final SVD logNormalsSVD = SVDFactory.createSVD(logNormals.counts(), ctx);
        final long svdEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished the SVD decomposition of the log-normal counts. Elapse of %d seconds", (svdEndTime - svdStartTime) / 1000));

        final int numberOfEigenSamples = determineNumberOfEigenSamples(requestedNumberOfEigenSamples, numberOfCountColumns, logNormalsSVD, logger);
        logger.info(String.format("Including %d eigen samples in the reduced PoN", numberOfEigenSamples));
        final RealMatrix logNormalsV = logNormalsSVD.getV();
        final RealMatrix logNormalsEigenV = logNormalsV.getSubMatrix(0, numberOfCountColumns - 1, 0, numberOfEigenSamples - 1);

        final RealMatrix reducedCounts = logNormalCounts.multiply(logNormalsEigenV);

        // The Pseudo inverse comes nearly for free from having run the SVD decomposition.
        final RealMatrix logNormalsPseudoInverse = logNormalsSVD.getPinv();

        logger.info("Calculating the reduced PoN inverse matrix...");
        final long riStartTime = System.currentTimeMillis();
        final RealMatrix reducedCountsInverse = SVDFactory.createSVD(reducedCounts, ctx).getPinv();
        final long riEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished calculating the reduced PoN inverse matrix. Elapse of %d seconds", (riEndTime - riStartTime) / 1000));
        return new ReductionResult(logNormalsPseudoInverse, reducedCounts, reducedCountsInverse, logNormalsSVD.getSingularValues());
    }

    /**
     * Calculate the number of eigen sample given the user preferences and the
     * result of the log-normals SVD.
     *
     * @param requestedNumberOfEigenSamples the user requested eigen values (empty means the user didn't specify any in particular).
     * @param numberOfCountColumns number of count columns in the original input.
     * @param logNormalsSVD SVD results on the log-normal counts.
     * @return always greater than 0.
     */
    @VisibleForTesting
    static int determineNumberOfEigenSamples(final OptionalInt requestedNumberOfEigenSamples, final int numberOfCountColumns, final SVD logNormalsSVD, final Logger logger) {
        final int numberOfEigenSamples;
        if (requestedNumberOfEigenSamples.isPresent()) {
            if (requestedNumberOfEigenSamples.getAsInt() > numberOfCountColumns) {
                logger.warn(String.format("The number of requested eigen samples (%d) is larger than the available number of read count columns after filtering (%d), thus we will have to use the latter.", requestedNumberOfEigenSamples.getAsInt(), numberOfCountColumns));
            }
            numberOfEigenSamples = Math.min(requestedNumberOfEigenSamples.getAsInt(), numberOfCountColumns);
        } else {
            final double[] singularValues = logNormalsSVD.getSingularValues();
            final double mean = MathUtils.mean(singularValues, 0, singularValues.length);
            final double eigenValueCutoff = mean * JOLLIFES_RULE_MEAN_FACTOR; // Jollife's less strict version of Kaiser' Rule.
            numberOfEigenSamples = (int) DoubleStream.of(singularValues).filter(v -> v > eigenValueCutoff).count();
            logger.info(String.format("Jollife's rule produced %d eigen samples out of %d possibles for the reduced PoN", numberOfEigenSamples, numberOfCountColumns));
        }
        return numberOfEigenSamples;
    }

    /**
     * Calculates the median of column medians and subtract it from all counts.
     * @param readCounts the input counts to center.
     * @return the median of medians that has been subtracted from all counts.
     */
    @VisibleForTesting
    static double subtractBGSCenter(final ReadCountCollection readCounts, final Logger logger) {
        final RealMatrix counts = readCounts.counts();
        final int columnCount = counts.getColumnDimension();
        final double[] columnMedians = new double[columnCount];
        final double[][] data = counts.getData();
        final Median medianCalculator = new Median();
        final double[] columnData = new double[data.length];
        for (int i = 0; i < columnCount; i++) {
            for (int j = 0; j < data.length; j++) {
                columnData[j] = data[j][i];
            }
            columnMedians[i] = medianCalculator.evaluate(columnData);
        }
        final double medianOfMedians = medianCalculator.evaluate(columnMedians);
        counts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) {
                return value - medianOfMedians;
            }
        });
        logger.info(String.format("Counts centered around the BGS center %.2f", medianOfMedians));
        return medianOfMedians;
    }

    /**
     * Final pre-panel normalization that consists on dividing all counts by the median of
     * its column and log it with base 2.
     * <p>
     *     The normalization occurs in-place.
     * </p>
     *
     * @param readCounts the input counts to normalize.
     */
    @VisibleForTesting
    static void normalizeAndLogReadCounts(final ReadCountCollection readCounts, final Logger logger) {
        final RealMatrix counts = readCounts.counts();
        final int columnCount = counts.getColumnDimension();
        final Median medianCalculator = new Median();
        for (int i = 0; i < columnCount; i++) {
            final double[] colCounts = counts.getColumn(i);
            final double median = medianCalculator.evaluate(colCounts);
            final double invMedian = 1.0 / median;
            for (int j = 0; j < colCounts.length; j++) {
                colCounts[j] = Math.log(Math.max(EPSILON, colCounts[j] * invMedian)) * INV_LN_2;
            }
            counts.setColumn(i, colCounts);
        }
        logger.info("Counts normalized by the column median and changed into bits.");
    }

    /**
     * Checks that the user input {@link #outlierTruncatePercentileThresh} is correct.
     * @return the user input {@link #outlierTruncatePercentileThresh}.
     */
    private double checkCountTruncatePercentile() {
        if (outlierTruncatePercentileThresh < 0 || outlierTruncatePercentileThresh > 50.0 || Double.isNaN(outlierTruncatePercentileThresh)) {
            throw new UserException.BadArgumentValue(COUNT_TRUNCATE_PERCENTILE_FULL_NAME, "the value must be in the [0, 50.0) range");
        }
        return outlierTruncatePercentileThresh;
    }

    /**
     * Impute zero counts to the median of non-zero values in the enclosing target row.
     *
     * <p>The imputation is done in-place, thus the input matrix is well be modified as a result of this call.</p>
     *
     * @param readCounts the input and output read-count matrix.
     */
    @VisibleForTesting
    static void imputeZerosCounts(final ReadCountCollection readCounts, final Logger logger) {

        final RealMatrix counts = readCounts.counts();
        final int targetCount = counts.getRowDimension();
        final int columnCount = counts.getColumnDimension();

        // Declared once outside loop to reuse across rows:
        final int[] zeroIndexes = new int[columnCount];
        final double[] nonZeroValues = new double[columnCount];
        final Median medianCalculator = new Median();
        int totalCounts = 0;
        int totalZeroCounts = 0;

        // for each row aka target:
        for (int i = 0; i < targetCount; i++) {
            // We classify columns as having a zero or not.
            // we accumulate nonZero values in nonZeroValues to calculate the median
            // where as we accumulate zero indexes that need to be updated later.
            final double[] rowCounts = counts.getRow(i);
            int zeroCount = 0;
            int nonZeroCount = 0;
            for (int j = 0; j <  rowCounts.length; j++) {
                if (rowCounts[j] != 0.0) {
                    nonZeroValues[nonZeroCount++] = rowCounts[j];
                } else {
                    zeroIndexes[zeroCount++] = j;
                }
            }
            totalCounts += rowCounts.length;
            if (zeroCount > 0) { // We only have to do anything if if there is any zeros in the row.
                final double nonZeroMedian = medianCalculator.evaluate(nonZeroValues, 0, nonZeroCount);
                for (int j = 0; j < zeroCount; j++) {
                    counts.setEntry(i, zeroIndexes[j], nonZeroMedian);
                }
                totalZeroCounts += zeroCount;
            }
        }
        if (totalZeroCounts > 0) {
            logger.info(String.format("Some 0.0 counts, %d of %d (%.2f%%), where imputed to their enclosing target non-zero median", totalZeroCounts, totalZeroCounts, 100.0 * (totalZeroCounts / totalCounts)));
        } else {
            logger.info("No count is 0.0 thus no count needed to be imputed.");
        }
    }
    /**
     * Truncates the extreme count values in the input read-count collection.
     * Values are forced to be bound by the percentile indicated with the input {@code percentile} which must be
     * in the range [0 .. 50.0]. Values under that percentile and the complementary (1 - percentile) are set to the
     * corresponding threshold value.
     *
     * <p>The imputation is done in-place, thus the input matrix is well be modified as a result of this call.</p>
     *
     * @param readCounts the input and output read-count matrix.
     */
    @VisibleForTesting
    static void truncateExtremeCounts(final ReadCountCollection readCounts, final double percentile, final Logger logger) {

        final RealMatrix counts = readCounts.counts();
        final int targetCount = counts.getRowDimension();
        final int columnCount = counts.getColumnDimension();
        final double[] values = new double[targetCount * columnCount];
        for (int i = 0; i < targetCount; i++) {
            final double[] rowCounts = counts.getRow(i);
            System.arraycopy(rowCounts, 0, values, i * columnCount, columnCount);
        }
        final Percentile bottomPercentileEvaluator = new Percentile(percentile);
        final Percentile topPercentileEvaluator = new Percentile(100.0 - percentile);
        final double bottomPercentileThreshold = bottomPercentileEvaluator.evaluate(values);
        final double topPercentileThreshold = topPercentileEvaluator.evaluate(values);
        long totalCounts = 0;
        long bottomTruncatedCounts = 0;
        long topTruncatedCounts = 0;

        for (int i = 0; i < targetCount; i++) {
            final double[] rowCounts = counts.getRow(i);
            for (int j = 0; j < columnCount; j++) {
                final double count = rowCounts[j];
                totalCounts++;
                if (count < bottomPercentileThreshold) {
                    counts.setEntry(i, j, bottomPercentileThreshold);
                    bottomTruncatedCounts++;
                } else if (count > topPercentileThreshold) {
                    counts.setEntry(i, j, topPercentileThreshold);
                    topTruncatedCounts++;
                }
            }
        }
        if (topTruncatedCounts == 0 && bottomTruncatedCounts == 0) {
            logger.info(String.format("None of the %d counts where truncated as they all fall in the non-extreme range [%.2f, %.2f]", totalCounts, bottomPercentileThreshold, topPercentileThreshold));
        } else {
            final double truncatedPercentage = ((double)(topTruncatedCounts + bottomTruncatedCounts) / totalCounts) * 100;
            logger.info(String.format("Some counts, %d out of %d (%.2f%%), where truncated as they fall out of the non-extreme range [%.2f, %.2f]",
                    topTruncatedCounts + bottomTruncatedCounts, totalCounts, truncatedPercentage, bottomPercentileThreshold, topPercentileThreshold));
        }
    }

    /**
     * Creates a new read-count collection that is a copy of the input but
     * dropping columns with extreme median coverage.
     *
     * @param readCounts the input read-counts.
     * @param extremeColumnMedianCountPercentileThreshold bottom percentile to use as a exclusion threshold.
     * @return never {@code null}. It might be a reference to the input read-counts if
     * there is no columns to be dropped.
     */
    @VisibleForTesting
    static ReadCountCollection removeColumnsWithExtremeMedianCounts(final ReadCountCollection readCounts, final double extremeColumnMedianCountPercentileThreshold, final Logger logger) {
        final List<String> columnNames = readCounts.columnNames();
        final double[] columnMedians = new double[columnNames.size()];
        final RealMatrix counts = readCounts.counts();
        final Median medianCalculator = new Median();
        for (int i = 0; i < columnMedians.length; i++) {
            final double[] columnCounts = counts.getColumn(i);
            columnMedians[i] = medianCalculator.evaluate(columnCounts);
        }

        // Calculate thresholds:
        final double bottomExtremeThreshold = new Percentile(extremeColumnMedianCountPercentileThreshold).evaluate(columnMedians);
        final double topExtremeThreshold = new Percentile(100 - extremeColumnMedianCountPercentileThreshold).evaluate(columnMedians);

        // Determine kept and dropped column sets.
        final Set<String> columnsToKeep = new HashSet<>(readCounts.columnNames().size());
        final Map<String, Double> columnsToDrop = new LinkedHashMap<>((int) (extremeColumnMedianCountPercentileThreshold * 4.0 / 100.0));
        for (int i = 0; i < columnMedians.length; i++) {
            if (columnMedians[i] >= bottomExtremeThreshold && columnMedians[i] <= topExtremeThreshold) {
                columnsToKeep.add(columnNames.get(i));
            } else {
                columnsToDrop.put(columnNames.get(i), columnMedians[i]);
            }
        }

        // log and drop columns if it applies
        if (columnsToKeep.isEmpty()) {
            throw new UserException.BadInput("No column count left after applying the extreme counts outlier filter");
        } else if (columnsToKeep.size() == columnNames.size()) {
            logger.info(String.format("No column dropped due to extreme counts outside [%.10f, %.10f]", bottomExtremeThreshold, topExtremeThreshold));
            return readCounts;
        } else {
            logger.info(String.format("Some columns dropped, %d of %d, as they are classified as having extreme median counts across targets outside [%.10f, %.10f]: e.g. %s",
                    columnNames.size() - columnsToKeep.size(), columnNames.size(), bottomExtremeThreshold, topExtremeThreshold,
                    columnsToDrop.entrySet().stream().limit(10).map(kv -> kv.getKey() + " (" + kv.getValue() + ")").collect(Collectors.joining(", "))));
            return readCounts.subsetColumns(columnsToKeep);
        }

    }

    /**
     * Remove targets that have too many counts equal to 0.
     * <p>
     *     It will return a copy of the input read-count collection with such targets dropped.
     * </p>
     *
     * @param readCounts the input read counts.
     * @param maximumTargetZeros maximum number of counts equal to 0. per target tolerated.
     * @return never {@code null}. It might be a reference to the input read-counts if there is
     *   is no target to be dropped.
     */
    @VisibleForTesting
    static ReadCountCollection removeTargetsWithTooManyZeros(final ReadCountCollection readCounts, final int maximumTargetZeros, final Logger logger) {
        final Set<Target> targetsToKeep = new HashSet<>(readCounts.targets().size());
        final RealMatrix counts = readCounts.counts();
        for (int i = 0; i < counts.getRowDimension(); i++) {
            final double[] rowCounts = counts.getRow(i);
            int zeroCounts = 0;
            for (final double rowCount : rowCounts) {
                if (rowCount == 0.0) {
                    zeroCounts++;
                }
            }
            if (zeroCounts <= maximumTargetZeros) {
                targetsToKeep.add(readCounts.targets().get(i));
            }
        }

        final int targetsToDropCount = readCounts.targets().size() - targetsToKeep.size();
        if (targetsToDropCount == 0) {
            logger.info(
                    String.format("No target drop due to too many zero counts; there is no target with a large number of columns with zero counts (<= %d of %d) ", maximumTargetZeros, readCounts.targets().size()));
            return readCounts;
        } else if (targetsToDropCount == readCounts.targets().size()) {
            throw new UserException.BadInput("the number of zeros per target in the input is too large resulting in all targets being dropped");
        } else {
            logger.info(
                    String.format("Some targets dropped, %d of %d, as they have too many zeros (> %d of %d)."
                            , targetsToDropCount, readCounts.columnNames().size(), maximumTargetZeros, readCounts.targets().size()));
            return readCounts.subsetTargets(targetsToKeep);
        }
    }

    /**
     * Remove columns that have too many counts equal to 0.
     * <p>
     *     It will return a copy of the input read-count collection with such columns dropped.
     * </p>
     *
     * @param readCounts the input read counts.
     * @param maximumColumnZeros maximum number of counts equal to 0. per column tolerated.
     * @return never {@code null}. It might be a reference to the input read-counts if there is
     *   is no column to be dropped.
     */
    @VisibleForTesting
    static ReadCountCollection removeColumnsWithTooManyZeros(final ReadCountCollection readCounts, final int maximumColumnZeros, final Logger logger) {
        final Set<String> columnsToKeep = new HashSet<>(readCounts.columnNames().size());
        final RealMatrix counts = readCounts.counts();
        for (int i = 0; i < counts.getColumnDimension(); i++) {
            final double[] columnCounts = counts.getColumn(i);
            int zeroCounts = 0;
            for (final double columnCount : columnCounts) {
                if (columnCount == 0.0) {
                    zeroCounts++;
                }
            }
            if (zeroCounts <= maximumColumnZeros) {
                columnsToKeep.add(readCounts.columnNames().get(i));
            }
        }

        final int columnsToDropCount = readCounts.columnNames().size() - columnsToKeep.size();
        if (columnsToDropCount == 0) {
            logger.info(
                    String.format("No count column dropped due to zero counts; there is no column with a large number of targets with zero counts (<= %d of %d) ", maximumColumnZeros, readCounts.targets().size()));
            return readCounts;
        } else if (columnsToDropCount == readCounts.columnNames().size()) {
            throw new UserException.BadInput("The number of zeros per count column is too large resulting in all count columns to be dropped");
        } else {
            logger.info(
                    String.format("Some counts columns dropped, %d of %d, as they have too many targets with zeros (> %d of %d)."
                            , columnsToDropCount, readCounts.columnNames().size(), maximumColumnZeros, readCounts.targets().size()));
            return readCounts.subsetColumns(columnsToKeep);
        }
    }

    /**
     * Calculates the absolute number of zero tolerated given the user requested percentage
     * and the number of counts involved.
     * @param totalCounts number of counts to check.
     * @param parameterValue a value from 0 to 100, otherwise a exception will be thrown.
     * @param exceptionFactory return exception to throw in case the value provided is invalid.
     * @return a value between 0 and {@code totalCounts}.
     */
    private int calculateMaximumZerosPercentage(final int totalCounts, final double parameterValue, final Supplier<UserException.BadArgumentValue> exceptionFactory) {
        if (parameterValue < 0 || parameterValue > 100 || Double.isNaN(parameterValue)) {
            throw exceptionFactory.get();
        } else {
            return (int) Math.ceil(totalCounts * parameterValue / 100.0);
        }
    }

    /**
     * Calculate the maximum number of target counts that can be zeros.
     * @param columnCount number of columns.
     * @return a number between 0 and {@code columnCount}.
     */
    private int calculateMaximumTargetZerosPercentage(final int columnCount) {
        return calculateMaximumZerosPercentage(columnCount, maximumPercentZerosInTarget,
                () -> new UserException.BadArgumentValue(MAXIMUM_PERCENT_ZEROS_IN_TARGET_FULL_NAME,
                    String.format("it must have a value in (0, 100) excluding limits: %.2f ", maximumPercentZerosInTarget)));
    }

    /**
     * Calculate the maximum number of column counts that can be zeros.
     * @param targetCount number of targets.
     * @return a number between 0 and {@code targetCount}.
     */
    private int calculateMaximumColumnZerosPercentage(final int targetCount) {
        return calculateMaximumZerosPercentage(targetCount, maximumPercentZerosInColumn,
                () -> new UserException.BadArgumentValue(MAXIMUM_PERCENT_ZEROS_IN_COLUMN_FULL_NAME,
                        String.format("it must have a value in (0, 100) excluding limits: %.2f ", maximumPercentZerosInColumn)));
    }

    /**
     * Normalizes read-counts values by the target factors in-place.
     * <p>
     *    Thi basically consists in dividing each count by the corresponding target factor
     *    (which in turn is typically the median across columns although this is not assured by this method.).
     * </p>
     * @param readCounts read-counts to normalize.
     * @param targetFactors the target factors
     */
    @VisibleForTesting
    static void normalizeReadCountsByTargetFactors(final ReadCountCollection readCounts, final double[] targetFactors) {
        final RealMatrix counts = readCounts.counts();
        final int targetCount = counts.getRowDimension();
        final int columnCount = counts.getColumnDimension();
        for (int i = 0; i < targetCount; i++) {
            final double factorInverse = 1.0 / targetFactors[i];

            for (int j = 0; j < columnCount; j++) {
                counts.setEntry(i, j, counts.getEntry(i, j) * factorInverse);
            }
        }
    }

    private double calculateExtremeColumnMedianCountsPercentileThreshold() {
        if (columnExtremeThresholdPercentile < 0 || columnExtremeThresholdPercentile > 50 || Double.isNaN(columnExtremeThresholdPercentile)) {
            throw new UserException.BadArgumentValue(COLUMN_EXTREME_THRESHOLD_PERCENTILE_FULL_NAME, "the value must be in the range [0, 50]");
        }
        return columnExtremeThresholdPercentile;
    }

    private double calculateTargetFactorsPercentileThreshold() {
        if (targetFactorThreshold < 0 || targetFactorThreshold > 100 || Double.isNaN(targetFactorThreshold)) {
            throw new UserException.BadArgumentValue(TARGET_FACTOR_THRESHOLD_PERCENTILE_FULL_NAME, "the value must be in the range [0, 100]");
        }
        return targetFactorThreshold;
    }

    /**
     * Calculate the target factors:
     *
     * <p>These are the median coverage per target</p>.
     * @param readCounts the input read counts.
     * @return never {@code null}, with as many elements as targets in {@code readCounts}.
     */
    private static double[] calculateTargetFactors(final ReadCountCollection readCounts) {
        final double[] result = new double[readCounts.targets().size()];
        final RealMatrix counts = readCounts.counts();
        final Median medianCalculator = new Median();
        for (int i = 0; i < result.length; i++) {
            result[i] = medianCalculator.evaluate(counts.getRow(i));
        }
        return result;
    }

    /**
     * Subsets targets in the input count to the usable ones based on the percentile threshold indicated
     * by the user.
     *
     * <p>
     *     It returns a pair of object, where the left one is the updated read-counts with only the usable
     *     targets, and the right one is the corresponding target factors.
     * </p>
     *
     * @param readCounts the input read-counts.
     * @param targetFactorPercentileThreshold the minimum median count percentile under which targets are not considered useful.
     * @return never {@code null}.
     */
    @VisibleForTesting
    static Pair<ReadCountCollection, double[]> subsetReadCountsToUsableTargets(final ReadCountCollection readCounts,
                                                                final double targetFactorPercentileThreshold, final Logger logger) {
        final double[] targetFactors = calculateTargetFactors(readCounts);
        final Percentile percentile = new Percentile(targetFactorPercentileThreshold);
        final double threshold = percentile.evaluate(targetFactors);
        final List<Target> targetByIndex = readCounts.targets();
        final Set<Target> result = IntStream.range(0, targetFactors.length).filter(i -> targetFactors[i] >= threshold)
                .mapToObj(targetByIndex::get)
                .collect(Collectors.toCollection(LinkedHashSet::new));
        if (result.size() == targetByIndex.size()) {
            logger.info(String.format("All %d targets are kept", targetByIndex.size()));
            return new ImmutablePair<>(readCounts, targetFactors);
        } else {
            final int discardedCount = targetFactors.length - result.size();
            logger.info(String.format("Discarded %d target(s) out of %d with factors below %.2g (%.2f percentile)", discardedCount, targetFactors.length, threshold, targetFactorPercentileThreshold  ));
            final double[] targetFactorSubset = DoubleStream.of(targetFactors).filter(i -> i >= threshold).toArray();
            return new ImmutablePair<>(readCounts.subsetTargets(result), targetFactorSubset);
        }
    }

    /**
     * Reads the read-counts from a file.
     * <p>
     *     The read-counts are subset
     * </p>
     *
     * @param inputFile the input file containing the read-counts.
     * @param targets the target collection to subset to.
     * @return never {@code null}.
     */
    private static ReadCountCollection readReadCountsFromFile(final File inputFile, final TargetCollection<Target> targets) {
        try {
            return ReadCountCollectionUtils.parse(inputFile, targets, true);
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile);
        }
    }
}
