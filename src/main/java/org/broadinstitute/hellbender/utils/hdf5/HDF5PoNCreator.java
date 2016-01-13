package org.broadinstitute.hellbender.utils.hdf5;

import com.google.api.client.repackaged.com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.google.common.primitives.Doubles;
import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
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
 * This class writes HDF5 PoN files.
 */
final public class HDF5PoNCreator {

    public static final double CURRENT_PON_VERSION = 5.0;

    public static final double JOLLIFES_RULE_MEAN_FACTOR = 0.7;

    private static final double INV_LN_2 = 1.0 / Math.log(2);

    public static final double EPSILON = 10e-10;

    private static final Logger logger = LogManager.getLogger(HDF5PoNCreator.class);

    private HDF5PoNCreator() {}

    /**
     * Create an HDF5 PoN file with the output from {@link CombineReadCounts}.
     *
     * When using this method in testing, see {@link CreatePanelOfNormals} for good default values.
     *
     * @param ctx If no spark context is available, specify {@code null}
     * @param inputPCovFile  pcov file from {@link CombineReadCounts} with samples as columns.  Never {@code null}
     * @param numberOfEigenSamples  desired number of eigen samples in the final PoN.  If missing, will apply a rule to determine.
     * @param sampleNameBlacklist  sample names in {@code inputPCovFile} that should be ignored.  Names must match exactly.  Never {@code null}.  Entires in this parameter that are not in {@code inputPCovFile} are ignored.
     * @param outputHDF5Filename  final filename for the output PoN HDF5 file
     * @param isDryRun  whether this is a dry run.  If you wish to do diagnostics and skip the actual writing of a PoN file, set this to true.
     * @param initialTargets  the set of targets that was used to generate the {@code inputPCovFile}.
     * @param targetFactorPercentileThreshold  the percentile of extreme target factor values to cut.
     * @param extremeColumnMedianCountPercentileThreshold Percentile to cut columns after they are ranked by median value
     * @param countTruncatePercentile percentile (on either end) to truncate extreme values (i.e. extreme high or low values will be set to the high or low count truncation value, which is based on this percentlie)
     * @param maximumPercentageZeroTargets  the maximum percentage of zero values in a target (across samples) before the target is filtered.
     * @param maximumPercentageZeroColumns  the maximum percentage of zero values in a sample (across targets) before the sample is filtered.
     */
    public static void createPoN(final JavaSparkContext ctx, final File inputPCovFile, final OptionalInt numberOfEigenSamples,
                                 final List<String> sampleNameBlacklist, final File outputHDF5Filename,
                                 final boolean isDryRun, final TargetCollection<Target> initialTargets,
                                 final double targetFactorPercentileThreshold,
                                 final double extremeColumnMedianCountPercentileThreshold,
                                 final double countTruncatePercentile, final double maximumPercentageZeroTargets,
                                 final double maximumPercentageZeroColumns) {

        Utils.nonNull(inputPCovFile);
        ParamUtils.inRange(targetFactorPercentileThreshold, 0, 100, "Target factor percentile threshold must be in range [0, 100]");
        ParamUtils.inRange(maximumPercentageZeroTargets, 0, 100, "Maximum percentage of zero-targets must be in range [0, 100]");
        ParamUtils.inRange(maximumPercentageZeroColumns, 0, 100, "Maximum percentage of zero-columns must be in range [0, 100]");
        ParamUtils.inRange(extremeColumnMedianCountPercentileThreshold, 0, 50, "Extreme column median percentile threshold must be in range [0, 50]");
        ParamUtils.inRange(countTruncatePercentile, 0, 50, "Count truncation threshold percentile threshold must be in range [0, 50]");
        Utils.nonNull(sampleNameBlacklist, "Blacklist sample list is null.  Please use empty list, instead.");
        Utils.nonNull(initialTargets, "Target collection is null, which is invalid.");

        // Remove low coverage targets:
        final Pair<ReadCountCollection, double[]> inputSubsetByUsableTargets = subsetReadCountsToUsableTargets(readReadCountsFromFile(inputPCovFile, initialTargets), targetFactorPercentileThreshold, logger);

        // Remove samples listed on the blacklist
        ReadCountCollection readCounts = inputSubsetByUsableTargets.getLeft();
        readCounts = readCounts.subsetColumns(Sets.difference(new HashSet<>(readCounts.columnNames()), new HashSet<>(sampleNameBlacklist)) );
        final double[] targetFactors = inputSubsetByUsableTargets.getRight();

        // Normalize read-counts by removing the target factor component:
        normalizeReadCountsByTargetFactors(readCounts, targetFactors);

        // Create and write the first part of the PoN HDF5 file.
        if (!isDryRun) {
            writeTargetFactorNormalizeReadCountsAndTargetFactors(outputHDF5Filename, readCounts, targetFactors, initialTargets.targets());
        }

        // Remove column and targets with too many zeros and targets with extreme median coverage:
        final int maximumColumnZerosCount = calculateMaximumZerosCount(readCounts.targets().size(), maximumPercentageZeroColumns);
        final int maximumTargetZerosCount = calculateMaximumZerosCount(readCounts.columnNames().size(), maximumPercentageZeroTargets);
        readCounts = removeColumnsWithTooManyZeros(readCounts, maximumColumnZerosCount, logger);
        readCounts = removeTargetsWithTooManyZeros(readCounts, maximumTargetZerosCount, logger);
        readCounts = removeColumnsWithExtremeMedianCounts(readCounts, extremeColumnMedianCountPercentileThreshold, logger);

        // Impute zero counts to be same as the median for the same target.
        // This happens in-place.
        imputeZeroCountsAsTargetMedians(readCounts, logger);

        // Truncate extreme count values.
        // This happens in-place.
        truncateExtremeCounts(readCounts, countTruncatePercentile, logger);

        // Normalize by the median and log scale the read counts.
        normalizeAndLogReadCounts(readCounts, logger);
        subtractBGSCenter(readCounts, logger);

        final ReductionResult reduction = calculateReducedPanelAndPInverses(readCounts, numberOfEigenSamples, logger, ctx);

        // Write the log-normals and reduced panel matrices.
        if (!isDryRun) {
            writeLogNormalsReducedPanel(outputHDF5Filename, readCounts, reduction, ctx);
        }
    }

    // Writes normalized and reduced panel matrices to the output HDF5 file.  Also, writes the target variance...
    private static void writeLogNormalsReducedPanel(final File outFile, final ReadCountCollection readCounts, final ReductionResult reduction, final JavaSparkContext ctx) {
        try (final HDF5File file = new HDF5File(outFile, HDF5File.OpenMode.READ_WRITE)) {
            final HDF5PoN pon = new HDF5PoN(file);
            pon.setPanelSampleNames(readCounts.columnNames());
            pon.setPanelTargetNames(readCounts.targets().stream()
                    .map(Target::getName)
                    .collect(Collectors.toList()));
            pon.setLogNormalCounts(readCounts.counts());
            pon.setLogNormalPInverseCounts(reduction.getPseudoInverse());
            pon.setReducedPanelCounts(reduction.getReducedCounts());
            pon.setReducedPanelPInverseCounts(reduction.getReducedInverse());
            pon.setPanelTargets(readCounts.targets());
            pon.setVersion(CURRENT_PON_VERSION);
            logger.info("Calculating and writing the target variances...");
            final double[] targetVariances = createTargetVariance(pon, ctx);
            pon.setTargetVariances(targetVariances);
        }
    }

    private static void writeTargetFactorNormalizeReadCountsAndTargetFactors(final File outFile, final ReadCountCollection readCounts, final double[] targetFactors, final List<Target> rawTargets) {
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
            pon.setTargets(readCounts.targets());
            pon.setRawTargets(rawTargets);
            pon.setRawTargetNames(rawTargets.stream().map(Target::getName).collect(Collectors.toList()));
        }
    }

    /**
     *  Creates a new PoN file using values from a given PoN file, except for reduction.  Does the reduction step a
     *   second time and populates a new PoN file.
     *
     * Please note that this code will redo the reduction even if the number of eigensamples is the same in both the input and desired output.
     *
     * @param ctx  {@code null} is okay if not using Spark
     * @param newNumberOfEigenSamples  desired number of eigensamples in the new PoN reduction.
     * @param inputHDF5Filename  input PoN file
     * @param outputHDF5Filename  output PoN file, which will be the same as the input, except for fields regarding the PoN reduction
     */
    public static void redoReduction(final JavaSparkContext ctx, final OptionalInt newNumberOfEigenSamples, final File inputHDF5Filename, final File outputHDF5Filename) {
        if (inputHDF5Filename.getAbsolutePath().equals(outputHDF5Filename.getAbsolutePath())) {
            throw new UserException.CouldNotCreateOutputFile(outputHDF5Filename, "Cannot create a new PoN overwriting an old one.");
        }

        Utils.regularReadableUserFile(inputHDF5Filename);

        try (final HDF5File ponReader = new HDF5File(inputHDF5Filename, HDF5File.OpenMode.READ_ONLY)) {
            final PoN inputPoN = new HDF5PoN(ponReader);
            final ReadCountCollection logNormalizedCounts = new ReadCountCollection(SetUniqueList.setUniqueList(inputPoN.getPanelTargets()), SetUniqueList.setUniqueList(new ArrayList<>(inputPoN.getPanelSampleNames())), inputPoN.getLogNormalizedCounts());
            final ReadCountCollection coverageProfile = new ReadCountCollection(SetUniqueList.setUniqueList(inputPoN.getTargets()), SetUniqueList.setUniqueList(new ArrayList<>(inputPoN.getSampleNames())), inputPoN.getNormalizedCounts());

            final ReductionResult newReduction = calculateReducedPanelAndPInverses(logNormalizedCounts, newNumberOfEigenSamples, logger, ctx);

            final PoN newPoN = new RamPoN(inputPoN, newReduction);
            writeTargetFactorNormalizeReadCountsAndTargetFactors(outputHDF5Filename, coverageProfile, newPoN.getTargetFactors().getColumn(0), inputPoN.getRawTargets());
            writeLogNormalsReducedPanel(outputHDF5Filename, logNormalizedCounts, newReduction, ctx);
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

        final double[] singularValues = logNormalsSVD.getSingularValues();
        final RealMatrix reducedCounts = logNormalsSVD.getU().getSubMatrix(0, logNormalCounts.getRowDimension()-1, 0, numberOfEigenSamples - 1).copy();
        reducedCounts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) { return singularValues[column]*value; }
        });

        // The Pseudo inverse comes nearly for free from having run the SVD decomposition.
        final RealMatrix logNormalsPseudoInverse = logNormalsSVD.getPinv();

        logger.info("Calculating the reduced PoN inverse matrix...");
        final long riStartTime = System.currentTimeMillis();
        final RealMatrix reducedCountsPseudoInverse = SVDFactory.createSVD(reducedCounts, ctx).getPinv();
        final long riEndTime = System.currentTimeMillis();
        logger.info(String.format("Finished calculating the reduced PoN inverse matrix. Elapse of %d seconds", (riEndTime - riStartTime) / 1000));
        return new ReductionResult(logNormalsPseudoInverse, reducedCounts, reducedCountsPseudoInverse, logNormalsSVD.getSingularValues());
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
        final Median medianCalculator = new Median();
        final double[] columnMedians = MatrixSummaryUtils.getColumnMedians(counts);

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
        final Median medianCalculator = new Median();
        final double[] columnMedians = IntStream.range(0, counts.getColumnDimension())
                .mapToDouble(col -> medianCalculator.evaluate(counts.getColumn(col))).toArray();

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

    private static long countZeroes(final double[] data) {
       return DoubleStream.of(data).filter(d -> d == 0.0).count();
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

        final Set<String> columnsToKeep = IntStream.range(0, counts.getColumnDimension()).boxed()
                .filter(i -> countZeroes(counts.getColumn(i)) <= maximumColumnZeros)
                .map(i -> readCounts.columnNames().get(i)).collect(Collectors.toSet());

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
     * Normalizes read-counts values by the target factors in-place.
     * <p>
     *    This basically consists in dividing each count by the corresponding target factor
     *    (which in turn is typically the median across columns although this is not assured by this method.).
     * </p>
     * @param readCounts read-counts to normalize.
     * @param targetFactors the target factors
     */
    @VisibleForTesting
    static void normalizeReadCountsByTargetFactors(final ReadCountCollection readCounts, final double[] targetFactors) {
        final RealMatrix counts = readCounts.counts();
        if (targetFactors.length != counts.getRowDimension()) {
            throw new GATKException("Number of target factors does not correspond to the number of rows.");
        }
        // Divide all counts by the target factor for the row.
        counts.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(final int row, final int column, final double value) {
                return value / targetFactors[row];
            }
        });
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

    /**
     * Calculate the target factors:
     *
     * <p>These are the median coverage per target</p>.
     * @param readCounts the input read counts.
     * @return never {@code null}, with as many elements as targets in {@code readCounts}.
     */
    private static double[] calculateTargetFactors(final ReadCountCollection readCounts) {
        return MatrixSummaryUtils.getRowMedians(readCounts.counts());
    }

    /**
     * Calculates the absolute number of zeroes tolerated given the user requested percentage
     * and the number of counts involved.
     * @param totalCounts number of counts to check.
     * @param percentage a value from 0 to 100, otherwise a exception will be thrown.
     * @return a value between 0 and {@code totalCounts}.
     */
    private static int calculateMaximumZerosCount(final int totalCounts, final double percentage) {
        return (int) Math.ceil(totalCounts * percentage / 100.0);
    }

    /**
     * Determine the variance for each target in the PoN (panel targets).
     *
     * @param pon This pon must have the normalizedCounts populated.  Never {@code null}
     * @return array of doubles where each double corresponds to a target in the PoN (panel targets)
     */
    private static double[] createTargetVariance(final PoN pon, final JavaSparkContext ctx) {
        final TangentNormalizationResult allNormals = TangentNormalizer.tangentNormalizeNormalsInPoN(pon, ctx);
        final RealMatrix allSampleProjectedTargets = allNormals.getTangentNormalized().counts();

        return calculateRowVariances(allSampleProjectedTargets);
    }

    /**
     * Calculate the variance for each row in a matrix.
     *
     * @param m Not {@code null}
     * @return Not {@code null}.  double array with length {@code m.getRowDimension()}
     */
    static double[] calculateRowVariances(final RealMatrix m) {
        Utils.nonNull(m, "Cannot calculate variances for a null matrix.");
        return MatrixSummaryUtils.getRowVariances(m);
    }

    /**
     * Create an array of target weights that has the same length as the PoN panel Targets.  These will correspond to
     * the Targets returned by {@link PoN :: getPanelTargets}
     *
     * @param pon Never {@code null}
     * @return Never {@code null}
     */
    public static double[] calculateTargetWeights(final PoN pon) {
        Utils.nonNull(pon, "PoN cannot be null");
        final double[] targetVariances = pon.getTargetVariances();
        final double[] targetWeights = new double[targetVariances.length];
        for (int i = 0; i < targetVariances.length; i++) {
            targetWeights[i] = 1/targetVariances[i];
        }
        return targetWeights;
    }
}
