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
     * {@link HDF5PCACoveragePoNCreationUtils#removeColumnsWithTooManyZeros},
     * {@link HDF5PCACoveragePoNCreationUtils#removeTargetsWithTooManyZeros},
     * {@link HDF5PCACoveragePoNCreationUtils#removeColumnsWithExtremeMedianCounts},
     * {@link HDF5PCACoveragePoNCreationUtils#imputeZeroCountsAsTargetMedians}, and
     * {@link HDF5PCACoveragePoNCreationUtils#truncateExtremeCounts}.
     */
    private static ReadCountCollection cleanNormalizedCounts(final ReadCountCollection readCounts, final Logger logger,
                                                             double maximumPercentageZeroColumns,
                                                             double maximumPercentageZeroTargets,
                                                             double extremeColumnMedianCountPercentileThreshold,
                                                             double countTruncatePercentile) {
        // Remove column and targets with too many zeros and targets with extreme median coverage:
        final int maximumColumnZerosCount = calculateMaximumZerosCount(readCounts.targets().size(), maximumPercentageZeroColumns);
        final int maximumTargetZerosCount = calculateMaximumZerosCount(readCounts.columnNames().size(), maximumPercentageZeroTargets);
        ReadCountCollection cleanedCounts = removeColumnsWithTooManyZeros(readCounts, maximumColumnZerosCount, logger);
        cleanedCounts = removeTargetsWithTooManyZeros(cleanedCounts, maximumTargetZerosCount, logger);
        cleanedCounts = removeColumnsWithExtremeMedianCounts(cleanedCounts, extremeColumnMedianCountPercentileThreshold, logger);

        // Impute zero counts to be same as the median for the same target.
        // This happens in-place.
        imputeZeroCountsAsTargetMedians(cleanedCounts, logger);

        // Truncate extreme count values.
        // This happens in-place.
        truncateExtremeCounts(cleanedCounts, countTruncatePercentile, logger);
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
        final RealMatrix counts = readCounts.counts();

        final Set<String> columnsToKeep = IntStream.range(0, counts.getColumnDimension())
                .filter(i -> countZeroes(counts.getColumn(i)) <= maximumColumnZeros)
                .mapToObj(i -> readCounts.columnNames().get(i)).collect(Collectors.toSet());

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
        final RealMatrix counts = readCounts.counts();

        final Set<Target> targetsToKeep = IntStream.range(0, counts.getRowDimension()).boxed()
                .filter(i -> countZeroes(counts.getRow(i)) <= maximumTargetZeros)
                .map(i -> readCounts.targets().get(i)).collect(Collectors.toSet());

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

    private static long countZeroes(final double[] data) {
        return DoubleStream.of(data).filter(d -> d == 0.0).count();
    }

    /**
     * Creates a new read-count collection that is a copy of the input but dropping columns with extreme median coverage.
     *
     * @param readCounts the input read-counts.
     * @param extremeColumnMedianCountPercentileThreshold bottom percentile to use as an exclusion threshold.
     * @return never {@code null}. It might be a reference to the input read-counts if
     * there are no columns to be dropped.
     */
    @VisibleForTesting
    static ReadCountCollection removeColumnsWithExtremeMedianCounts(final ReadCountCollection readCounts, final double extremeColumnMedianCountPercentileThreshold, final Logger logger) {
        final List<String> columnNames = readCounts.columnNames();
        final RealMatrix counts = readCounts.counts();
        final double[] columnMedians = MatrixSummaryUtils.getColumnMedians(counts);

        // Calculate thresholds:
        final double bottomExtremeThreshold = new Percentile(extremeColumnMedianCountPercentileThreshold).evaluate(columnMedians);
        final double topExtremeThreshold = new Percentile(100 - extremeColumnMedianCountPercentileThreshold).evaluate(columnMedians);

        // Determine kept and dropped column sets.
        final Set<String> columnsToKeep = new HashSet<>(readCounts.columnNames().size());
        final int initialMapSize = ((int) (2. * extremeColumnMedianCountPercentileThreshold / 100.)) * readCounts.columnNames().size();
        final Map<String, Double> columnsToDrop = new LinkedHashMap<>(initialMapSize);
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
     * Impute zero counts to the median of non-zero values in the enclosing target row.
     *
     * <p>The imputation is done in-place, thus the input matrix is well be modified as a result of this call.</p>
     *
     * @param readCounts the input and output read-count matrix.
     */
    @VisibleForTesting
    static void imputeZeroCountsAsTargetMedians(final ReadCountCollection readCounts, final Logger logger) {

        final RealMatrix counts = readCounts.counts();
        final int targetCount = counts.getRowDimension();

        final Median medianCalculator = new Median();
        int totalCounts = counts.getColumnDimension() * counts.getRowDimension();

        // Get the number of zeroes contained in the counts.
        final long totalZeroCounts = IntStream.range(0, targetCount)
                .mapToLong(t -> DoubleStream.of(counts.getRow(t))
                        .filter(c -> c == 0.0).count()).sum();

        // Get the median of each row, not including any entries that are zero.
        final double[] medians = IntStream.range(0, targetCount)
                .mapToDouble(t -> medianCalculator.evaluate(
                        DoubleStream.of(counts.getRow(t))
                                .filter(c -> c != 0.0)
                                .toArray()
                )).toArray();

        // Change any zeros in the counts to the median for the row.
        counts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) {
                return value != 0 ? value : medians[row];
            }
        });

        if (totalZeroCounts > 0) {
            logger.info(String.format("Some 0.0 counts, %d of %d (%.2f%%), were imputed to their enclosing target non-zero median", totalZeroCounts, totalZeroCounts, 100.0 * (totalZeroCounts / totalCounts)));
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
     * <p>The imputation is done in-place, thus the input matrix is modified as a result of this call.</p>
     *
     * @param readCounts the input and output read-count matrix.
     */
    @VisibleForTesting
    static void truncateExtremeCounts(final ReadCountCollection readCounts, final double percentile, final Logger logger) {

        final RealMatrix counts = readCounts.counts();
        final int targetCount = counts.getRowDimension();
        final int columnCount = counts.getColumnDimension();

        // Create a row major array of the counts.
        final double[] values = Doubles.concat(counts.getData());

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
