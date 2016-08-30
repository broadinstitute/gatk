package org.broadinstitute.hellbender.tools.pon.coverage.pca;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.google.common.primitives.Doubles;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hdf5.HDF5File;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.MatrixSummaryUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.svd.SVD;
import org.broadinstitute.hellbender.utils.svd.SVDFactory;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

/**
 * This class contains utility methods for creating and writing the PCA coverage panel of normals to an {@link HDF5PCACoveragePoN}.
 */
public final class HDF5PCACoveragePoNCreationUtils {
    public static final double JOLLIFES_RULE_MEAN_FACTOR = 0.7;
    public static final double EPSILON = 1E-9;
    private static final double INV_LN_2 = 1.0 / Math.log(2);

    private static final Logger logger = LogManager.getLogger(HDF5PCACoveragePoNCreationUtils.class);

    private HDF5PCACoveragePoNCreationUtils() {}

    /**
     * Create an HDF5 coverage PoN file with the output from {@link CombineReadCounts}.
     *
     * When using this method in testing, see {@link CreatePanelOfNormals} for good default values.
     * @param ctx If no spark context is available, specify {@code null}
     * @param outputHDF5Filename  final filename for the output PoN HDF5 file
     * @param openMode            desired {@link HDF5File.OpenMode} (if {@code HDF5File.OpenMode.READ_ONLY}, an exception will be thrown if not a dry run)
     * @param inputPCovFile  pcov file from {@link CombineReadCounts} with samples as columns.  Never {@code null}
     * @param initialTargets  the set of targets that was used to generate the {@code inputPCovFile}.
     * @param sampleNameBlacklist  sample names in {@code inputPCovFile} that should be ignored.  Names must match exactly.  Never {@code null}.  Entires in this parameter that are not in {@code inputPCovFile} are ignored.
     * @param targetFactorPercentileThreshold  the percentile of extreme target factor values to cut.
     * @param maximumPercentageZeroColumns  the maximum percentage of zero values in a sample (across targets) before the sample is filtered.
     * @param maximumPercentageZeroTargets  the maximum percentage of zero values in a target (across samples) before the target is filtered.
     * @param extremeColumnMedianCountPercentileThreshold Percentile to cut columns after they are ranked by median value
     * @param countTruncatePercentile percentile (on either end) to truncate extreme values (i.e. extreme high or low values will be set to the high or low count truncation value, which is based on this percentlie)
     * @param numberOfEigensamples  desired number of eigensamples in the final PoN.  If missing, will apply a rule to determine.
     * @param isDryRun  whether this is a dry run.  If you wish to do diagnostics and skip the actual writing of a PoN file, set this to true.
     */
    public static void create(final JavaSparkContext ctx,
                              final File outputHDF5Filename,
                              final HDF5File.OpenMode openMode,
                              final File inputPCovFile,
                              final TargetCollection<Target> initialTargets,
                              final List<String> sampleNameBlacklist,
                              final double targetFactorPercentileThreshold,
                              final double maximumPercentageZeroColumns,
                              final double maximumPercentageZeroTargets,
                              final double extremeColumnMedianCountPercentileThreshold,
                              final double countTruncatePercentile,
                              final OptionalInt numberOfEigensamples,
                              final boolean isDryRun) {
        Utils.nonNull(outputHDF5Filename);
        Utils.regularReadableUserFile(inputPCovFile);
        Utils.nonNull(initialTargets, "Target collection cannot be null.");
        Utils.nonNull(sampleNameBlacklist, "Blacklist sample list cannot be null.  Use empty list if no blacklisting is desired.");
        ParamUtils.inRange(targetFactorPercentileThreshold, 0, 100, "Target factor percentile threshold must be in range [0, 100].");
        ParamUtils.inRange(maximumPercentageZeroColumns, 0, 100, "Maximum percentage of zero-columns must be in range [0, 100].");
        ParamUtils.inRange(maximumPercentageZeroTargets, 0, 100, "Maximum percentage of zero-targets must be in range [0, 100].");
        ParamUtils.inRange(extremeColumnMedianCountPercentileThreshold, 0, 50, "Extreme column median percentile threshold must be in range [0, 50].");
        ParamUtils.inRange(countTruncatePercentile, 0, 50, "Count truncation threshold percentile threshold must be in range [0, 50].");
        Utils.nonNull(numberOfEigensamples, "Number of eigensamples cannot be null.");

        final ReadCountCollection inputPCov = readReadCountsFromFile(inputPCovFile, initialTargets);    //ignore missing targets

        final Pair<ReadCountCollection, double[]> inputSubsetByUsableTargets = subsetReadCountsToUsableTargets(inputPCov, targetFactorPercentileThreshold, logger);

        // Remove samples listed on the blacklist and normalize read-counts by removing the target factor component
        ReadCountCollection normalizedCounts = inputSubsetByUsableTargets.getLeft();
        normalizedCounts = normalizedCounts.subsetColumns(Sets.difference(new HashSet<>(normalizedCounts.columnNames()), new HashSet<>(sampleNameBlacklist)) );
        final double[] targetFactors = inputSubsetByUsableTargets.getRight();
        PCATangentNormalizationUtils.factorNormalize(normalizedCounts.counts(), targetFactors);

        // Impute zeros as median values for columns and targets with too many zeros, remove targets with extreme medians, and truncate extreme counts
        final ReadCountCollection logNormalizedCounts = cleanNormalizedCounts(normalizedCounts, logger, maximumPercentageZeroColumns, maximumPercentageZeroTargets, extremeColumnMedianCountPercentileThreshold, countTruncatePercentile);

        // Normalize by the median and log_2 scale the read counts.
        normalizeAndLogReadCounts(logNormalizedCounts, logger);

        // Subtract the median of medians
        subtractMedianOfMedians(logNormalizedCounts, logger);

        // Perform the SVD and calculate the pseudoinverse
        final ReductionResult reduction = calculateReducedPanelAndPInverses(logNormalizedCounts, numberOfEigensamples, logger, ctx);

        // Calculate the target variances
        final List<String> panelTargetNames = logNormalizedCounts.targets().stream().map(Target::getName).collect(Collectors.toList());
        final double[] targetVariances = calculateTargetVariances(normalizedCounts, panelTargetNames, reduction, ctx);

        // Write the PoN to HDF5 file
        if (!isDryRun) {
            HDF5PCACoveragePoN.write(outputHDF5Filename, openMode, initialTargets.targets(), normalizedCounts, logNormalizedCounts, targetFactors, targetVariances, reduction);
        }
    }

    /**
     * Creates a new PoN file using values from a given PoN file, except for reduction.  Does the reduction step a
     * second time and populates a new PoN file.
     *
     * Please note that this code will redo the reduction even if the number of eigensamples is the same in both the input and desired output.
     *
     * @param ctx  {@code null} is okay if not using Spark
     * @param newNumberOfEigensamples  desired number of eigensamples in the new PoN reduction.
     * @param inputHDF5Filename  input PoN file
     * @param outputHDF5Filename  output PoN file, which will be the same as the input, except for fields regarding the PoN reduction
     * @param openMode            desired {@link HDF5File.OpenMode} (if {@code HDF5File.OpenMode.READ_ONLY}, an exception will be thrown)
     */
    public static void redoReduction(final JavaSparkContext ctx, final OptionalInt newNumberOfEigensamples, final File inputHDF5Filename, final File outputHDF5Filename, final HDF5File.OpenMode openMode) {
        Utils.nonNull(newNumberOfEigensamples);
        Utils.regularReadableUserFile(inputHDF5Filename);
        Utils.nonNull(outputHDF5Filename);
        if (inputHDF5Filename.getAbsolutePath().equals(outputHDF5Filename.getAbsolutePath())) {
            throw new UserException.CouldNotCreateOutputFile(outputHDF5Filename, "Cannot create a new PoN overwriting an old one.");
        }

        try (final HDF5File ponReader = new HDF5File(inputHDF5Filename, HDF5File.OpenMode.READ_ONLY)) {
            final PCACoveragePoN inputPoN = new HDF5PCACoveragePoN(ponReader);
            final ReadCountCollection normalizedCounts = new ReadCountCollection(inputPoN.getTargets(), inputPoN.getSampleNames(), inputPoN.getNormalizedCounts());
            final ReadCountCollection logNormalizedCounts = new ReadCountCollection(inputPoN.getPanelTargets(), inputPoN.getPanelSampleNames(), inputPoN.getLogNormalizedCounts());
            final ReductionResult newReduction = calculateReducedPanelAndPInverses(logNormalizedCounts, newNumberOfEigensamples, logger, ctx);
            final List<String> panelTargetNames = logNormalizedCounts.targets().stream().map(Target::getName).collect(Collectors.toList());
            final double[] targetVariances = calculateTargetVariances(normalizedCounts, panelTargetNames, newReduction, ctx);

            HDF5PCACoveragePoN.write(outputHDF5Filename, openMode, inputPoN.getRawTargets(), normalizedCounts, logNormalizedCounts, inputPoN.getTargetFactors(), targetVariances, newReduction);
        }
    }

    /*===============================================================================================================*
     * PRIVATE METHODS (SOME VISIBLE FOR TESTING)                                                                    *
     * These methods perform all of the steps needed to calculate the fields of the coverage panel of normals.       *
     *===============================================================================================================*/

    private static ReadCountCollection readReadCountsFromFile(final File inputFile, final TargetCollection<Target> targets) {
        try {
            return ReadCountCollectionUtils.parse(inputFile, targets, true);    //ignore missing targets = true
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile);
        }
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
        final double threshold = new Percentile(targetFactorPercentileThreshold).evaluate(targetFactors);
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
     * Performs
     * {@link ReadCountCollectionUtils#removeColumnsWithTooManyZeros},
     * {@link ReadCountCollectionUtils#removeTargetsWithTooManyZeros},
     * {@link ReadCountCollectionUtils#removeColumnsWithExtremeMedianCounts},
     * {@link ReadCountCollectionUtils#imputeZeroCountsAsTargetMedians}, and
     * {@link ReadCountCollectionUtils#truncateExtremeCounts}.
     */
    private static ReadCountCollection cleanNormalizedCounts(final ReadCountCollection readCounts, final Logger logger,
                                                             double maximumPercentageZeroColumns,
                                                             double maximumPercentageZeroTargets,
                                                             double extremeColumnMedianCountPercentileThreshold,
                                                             double countTruncatePercentile) {
        // Remove column and targets with too many zeros and targets with extreme median coverage:
        final int maximumColumnZerosCount = calculateMaximumZerosCount(readCounts.targets().size(), maximumPercentageZeroColumns);
        final int maximumTargetZerosCount = calculateMaximumZerosCount(readCounts.columnNames().size(), maximumPercentageZeroTargets);
        ReadCountCollection cleanedCounts = ReadCountCollectionUtils.removeColumnsWithTooManyZeros(readCounts, maximumColumnZerosCount, logger);
        cleanedCounts = ReadCountCollectionUtils.removeTargetsWithTooManyZeros(cleanedCounts, maximumTargetZerosCount, logger);
        cleanedCounts = ReadCountCollectionUtils.removeColumnsWithExtremeMedianCounts(cleanedCounts, extremeColumnMedianCountPercentileThreshold, logger);

        // Impute zero counts to be same as the median for the same target.
        // This happens in-place.
        ReadCountCollectionUtils.imputeZeroCountsAsTargetMedians(cleanedCounts, logger);

        // Truncate extreme count values.
        // This happens in-place.
        ReadCountCollectionUtils.truncateExtremeCounts(cleanedCounts, countTruncatePercentile, logger);
        return cleanedCounts;
    }

    /**
     * Final pre-panel normalization that consists of dividing all counts by the median of
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
        final Median medianCalculator = new Median();

        final double[] medians = IntStream.range(0, counts.getColumnDimension()).mapToDouble(col -> medianCalculator.evaluate(counts.getColumn(col))).toArray();
        counts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) {
                return Math.log(Math.max(EPSILON, value / medians[column])) * INV_LN_2;
            }
        });

        logger.info("Counts normalized by the column median and log2'd.");
    }

    /**
     * Calculates the median of column medians and subtract it from all counts.
     * @param readCounts the input counts to center.
     * @return the median of medians that has been subtracted from all counts.
     */
    @VisibleForTesting
    static double subtractMedianOfMedians(final ReadCountCollection readCounts, final Logger logger) {
        final RealMatrix counts = readCounts.counts();
        final Median medianCalculator = new Median();
        final double[] columnMedians = MatrixSummaryUtils.getColumnMedians(counts);

        final double medianOfMedians = medianCalculator.evaluate(columnMedians);
        counts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) {
                return value - medianOfMedians;
            }
        });
        logger.info(String.format("Counts centered around the median of medians %.2f", medianOfMedians));
        return medianOfMedians;
    }

    /**
     * SVD and Pseudo inverse calculation.
     *
     * @param logNormalized the input counts for the SVD and reduction steps, fully normalized and already logged.
     * @param requestedNumberOfEigensamples user requested number of eigensamples for the reduced panel.
     * @return never {@code null}.
     */
    @VisibleForTesting
    static ReductionResult calculateReducedPanelAndPInverses(final ReadCountCollection logNormalized,
                                                             final OptionalInt requestedNumberOfEigensamples,
                                                             final Logger logger,
                                                             final JavaSparkContext ctx) {

        if (ctx == null) {
            logger.warn("No Spark context provided, not going to use Spark...");
        }

        final RealMatrix logNormalizedCounts = logNormalized.counts();
        final int numberOfCountColumns = logNormalizedCounts.getColumnDimension();

        logger.info("Starting the SVD decomposition of the log-normalized counts ...");
        final long svdStartTime = System.currentTimeMillis();
        final SVD logNormalizedSVD = SVDFactory.createSVD(logNormalized.counts(), ctx);
        final long svdEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished the SVD decomposition of the log-normal counts. Elapse of %d seconds", (svdEndTime - svdStartTime) / 1000));

        final int numberOfEigensamples = determineNumberOfEigensamples(requestedNumberOfEigensamples, numberOfCountColumns, logNormalizedSVD, logger);
        logger.info(String.format("Including %d eigensamples in the reduced PoN", numberOfEigensamples));

        final double[] singularValues = logNormalizedSVD.getSingularValues();
        final RealMatrix reducedCounts = logNormalizedSVD.getU().getSubMatrix(0, logNormalizedCounts.getRowDimension()-1, 0, numberOfEigensamples - 1).copy();
        reducedCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) { return singularValues[column]*value; }
        });

        // The Pseudo inverse comes nearly for free from having run the SVD decomposition.
        final RealMatrix logNormalizedPseudoInverse = logNormalizedSVD.getPinv();

        logger.info("Calculating the reduced PoN inverse matrix...");
        final long riStartTime = System.currentTimeMillis();
        final RealMatrix reducedCountsPseudoInverse = SVDFactory.createSVD(reducedCounts, ctx).getPinv();
        final long riEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished calculating the reduced PoN inverse matrix. Elapse of %d seconds", (riEndTime - riStartTime) / 1000));
        return new ReductionResult(logNormalizedPseudoInverse, reducedCounts, reducedCountsPseudoInverse, logNormalizedSVD.getSingularValues());
    }

    /**
     * Determine the variance for each target in the PoN (panel targets).
     *
     * @return      array of doubles where each double corresponds to a target in the PoN (panel targets)
     */
    private static double[] calculateTargetVariances(final ReadCountCollection normalizedCounts,
                                                     final List<String> panelTargetNames,
                                                     final ReductionResult reduction,
                                                     final JavaSparkContext ctx) {
        Utils.nonNull(panelTargetNames);
        Utils.nonNull(normalizedCounts);
        Utils.nonNull(reduction);
        final PCATangentNormalizationResult allNormals =
                PCATangentNormalizationUtils.tangentNormalizeNormalsInPoN(normalizedCounts, panelTargetNames, reduction.getReducedCounts(), reduction.getReducedPseudoInverse(), ctx);
        final RealMatrix allSampleProjectedTargets = allNormals.getTangentNormalized().counts();

        return MatrixSummaryUtils.getRowVariances(allSampleProjectedTargets);
    }

    /**
     * Calculate the target factors (median coverage per target).
     *
     * @param readCounts    the input read counts
     * @return              never {@code null}, with as many elements as targets in {@code readCounts}
     */
    private static double[] calculateTargetFactors(final ReadCountCollection readCounts) {
        Utils.nonNull(readCounts, "ReadCountCollection cannot be null.");
        return MatrixSummaryUtils.getRowMedians(readCounts.counts());
    }

    /**
     * Calculates the absolute number of zeroes tolerated given the user requested percentage
     * and the number of counts involved.
     *
     * @param totalCounts   a non-negative number of counts to check
     * @param percentage    a value in [0, 100]
     * @return              a value between 0 and {@code totalCounts}
     */
    private static int calculateMaximumZerosCount(final int totalCounts, final double percentage) {
        ParamUtils.isPositiveOrZero(totalCounts, "Total counts must be non-negative.");
        ParamUtils.inRange(percentage, 0, 100, "Percentage must be in [0, 100].");
        return (int) Math.ceil(totalCounts * percentage / 100.0);
    }

    /**
     * Calculate the number of eigensamples given the user preferences and the
     * result of the log-normals SVD.
     *
     * @param requestedNumberOfEigensamples the user requested eigenvalues (empty means the user didn't specify any in particular).
     * @param numberOfCountColumns number of count columns in the original input.
     * @param logNormalizedSVD SVD results on the log-normalized counts.
     * @return always greater than 0.
     */
    @VisibleForTesting
    static int determineNumberOfEigensamples(final OptionalInt requestedNumberOfEigensamples, final int numberOfCountColumns, final SVD logNormalizedSVD, final Logger logger) {
        final int numberOfEigensamples;
        if (requestedNumberOfEigensamples.isPresent()) {
            if (requestedNumberOfEigensamples.getAsInt() > numberOfCountColumns) {
                logger.warn(String.format("The number of requested eigensamples (%d) is larger than the available number of read count columns after filtering (%d), thus we will have to use the latter.", requestedNumberOfEigensamples.getAsInt(), numberOfCountColumns));
            }
            numberOfEigensamples = Math.min(requestedNumberOfEigensamples.getAsInt(), numberOfCountColumns);
        } else {
            final double[] singularValues = logNormalizedSVD.getSingularValues();
            final double mean = MathUtils.mean(singularValues, 0, singularValues.length);
            final double eigenvalueCutoff = mean * JOLLIFES_RULE_MEAN_FACTOR; // Jollife's less strict version of Kaiser' Rule.
            numberOfEigensamples = (int) DoubleStream.of(singularValues).filter(v -> v > eigenvalueCutoff).count();
            logger.info(String.format("Jollife's rule produced %d eigensamples out of %d possibles for the reduced PoN", numberOfEigensamples, numberOfCountColumns));
        }
        return numberOfEigensamples;
    }
}
