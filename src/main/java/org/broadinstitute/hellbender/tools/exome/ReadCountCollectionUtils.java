package org.broadinstitute.hellbender.tools.exome;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.primitives.Doubles;
import htsjdk.samtools.util.Locatable;
import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.samplenamefinder.SampleNameFinder;
import org.broadinstitute.hellbender.utils.MatrixSummaryUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

/**
 * Reads {@link ReadCountCollection} instances from a tab-separated text file and implements a number of read count filters.
 * <p>
 * The tab separated file consist of a header and body with the data.
 * </p>
 * <p>
 * The header consist of at least a line with column names optionally preceded by comment lines (starting with '#').
 * A part from target coordinates and name columns there should be at least on actual count column (sample, read-group or cohort).
 * but there could be more than one.
 * </p>
 * <p>
 * The body are the coordinates and counts for each target.
 * </p>
 * <p>
 * Example:
 * </p>
 * <pre>
 *     ##comment-line1  (optional)
 *     ##comment-line2  (optional)
 *     CONTIG   START   END NAME    SAMPLE1 SAMPLE2 SAMPLE3
 *     1    1000    1100    tgt_0   5   2   10
 *     1    2000    2200    tgt_1   1   2   2
 *     ...
 *     X    21300   21400   tgt_2311    10  3   7
 * </pre>
 * <p>
 * You may omit either the target name column (NAME) or some of the genomic interval columns (CONTIG, START and END)
 * but not both at the same time.
 * </p>
 * <p>
 * If the source omits the target name, a exonCollection should be provided in order to resolve the name given its coordinates
 * using {@link #parse(File, TargetCollection, boolean)}.
 * </p>
 * <p>
 * This class will check whether the content of the input file is well formatted and consistent
 * (e.g. counts are double values, each row have the same number of values, on for each column in the header,
 * and so forth).
 * </p>
 * <p>
 * If there is any formatting problems the appropriate exception will be thrown
 * as described in {@link #parse}.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class ReadCountCollectionUtils {

    // Prevents instantiation of the class.
    private ReadCountCollectionUtils() {}

    /**
     * Writes the content of a collection into a writer.
     *
     * @param writer         the output writer.
     * @param collection     the output collection.
     * @param headerComments header comments.
     * @throws IllegalArgumentException if any of the input parameters is {@code null}
     *                                  or {@code collection} has a mixture of targets with and without intervals
     *                                  defined.
     * @throws IOException              if there is some IO issue when writing into the output file.
     */
    public static void write(final Writer writer, final ReadCountCollection collection, final String... headerComments) throws IOException {
        Utils.nonNull(collection, "input collection cannot be null");
        Utils.nonNull(headerComments, "header comments cannot be null");
        final boolean withIntervals = collection.targets().stream().anyMatch(t -> t.getInterval() != null);
        final TableWriter<ReadCountRecord> tableWriter = withIntervals
                ? writerWithIntervals(writer, collection.columnNames()) : writerWithoutIntervals(writer, collection);
        performWriting(collection, tableWriter, headerComments);
    }

    /**
     * Writes the content of a collection into a file.
     *
     * @param file           the output file.
     * @param collection     the output collection.
     * @param headerComments header comments.
     * @throws IllegalArgumentException if any of the input parameters is {@code null}
     *                                  or {@code collection} has a mixture of targets with and without intervals
     *                                  defined.
     * @throws IOException              if there is some IO issue when writing into the output file.
     */
    public static void write(final File file, final ReadCountCollection collection, final String... headerComments) throws IOException {
        Utils.nonNull(file, "output file cannot be null");
        try (final Writer writer = new FileWriter(file)) {
            write(writer, collection, headerComments);
        }
    }

    private static void performWriting(ReadCountCollection collection, TableWriter<ReadCountRecord> tableWriter, String[] headerComments) throws IOException {
        // print the header comments
        for (final String comment : headerComments) {
            tableWriter.writeComment(comment);
        }
        final List<Target> targets = collection.targets();
        final RealMatrix counts = collection.counts();
        for (int i = 0; i < targets.size(); i++) {
            tableWriter.writeRecord(new ReadCountRecord(targets.get(i), counts.getRow(i)));
        }
    }

    /**
     * Creates a new table writer that will output the target intervals.
     * @param writer where to output the table formatted content.
     * @param countColumnNames list of count column names.
     * @return never {@code null}.
     * @throws IOException if there is some low level IO problem creating the writer.
     * @throws IllegalArgumentException if {@code countColumnNames} is {@code null}, contains
     *  {@code null} or a non valid count column name (e.g. a reserved word).
     */
    public static TableWriter<ReadCountRecord> writerWithIntervals(final Writer writer, final List<String> countColumnNames) throws IOException {

        final List<String> columnNames = new ArrayList<>();

        columnNames.add(TargetTableColumn.CONTIG.toString());
        columnNames.add(TargetTableColumn.START.toString());
        columnNames.add(TargetTableColumn.END.toString());
        columnNames.add(TargetTableColumn.NAME.toString());
        columnNames.addAll(Utils.nonNull(countColumnNames));
        final TableColumnCollection columns = new TableColumnCollection(columnNames);

        return new TableWriter<ReadCountRecord>(writer, columns) {
            @Override
            protected void composeLine(final ReadCountRecord record, final DataLine dataLine) {
                final SimpleInterval interval = record.getTarget().getInterval();
                if (interval == null) {
                    throw new IllegalStateException("invalid combination of targets with and without intervals defined");
                }
                dataLine.append(interval.getContig())
                        .append(interval.getStart())
                        .append(interval.getEnd())
                        .append(record.getTarget().getName());
                record.appendCountsTo(dataLine);
            }
        };
    }

    private static TableWriter<ReadCountRecord> writerWithoutIntervals(final Writer writer,
                                                                       final ReadCountCollection collection) throws IOException {

        final List<String> columnNames = new ArrayList<>();
        columnNames.add(TargetTableColumn.NAME.toString());
        columnNames.addAll(collection.columnNames());
        return createReadCountRecordTableWriterWithoutIntervals(writer, columnNames);
    }

    private static TableWriter<ReadCountRecord> createReadCountRecordTableWriterWithoutIntervals(final Writer writer, List<String> columnNames) throws IOException {
        final TableColumnCollection columns = new TableColumnCollection(columnNames);
        return new TableWriter<ReadCountRecord>(writer, columns) {

            @Override
            protected void composeLine(final ReadCountRecord record, final DataLine dataLine) {
                dataLine.append(record.getTarget().getName());
                record.appendCountsTo(dataLine);
            }
        };
    }

    /**
     * Reads the content from a reader into a {@link ReadCountCollection}.
     *
     * @param sourceReader the source reader
     * @param sourceName source name.
     * @return never {@code null}.
     * @throws IOException            if there was some problem reading the file contents.
     * @throws UserException.BadInput if there is some formatting issue in the source file contents. This includes
     *                                lack of target names in the source file.
     */
    public static ReadCountCollection parse(final Reader sourceReader, final String sourceName) throws IOException {
        return parse(sourceReader, sourceName, null, false);
    }

    /**
     * Reads the content of a file into a {@link ReadCountCollection}.
     *
     * @param file the source file.
     * @return never {@code null}.
     * @throws IOException            if there was some problem reading the file contents.
     * @throws UserException.BadInput if there is some formatting issue win the source file contents. This includes
     *                                lack of target names in the source file.
     */
    public static ReadCountCollection parse(final File file) throws IOException {
        return parse(file, null, false);
    }

    /**
     * Reads the content of a file into a {@link ReadCountCollection}.
     * <p>
     * If no target name is included in the input but intervals are present, the {@code exons} collection provided
     * will be utilized to resolve those names.
     * </p>
     *
     * @param file  the source file.
     * @param targets collection of exons (targets). This parameter can be {@code null}, to indicate that no exon
     *              collection is to be considered.
     * @param ignoreMissingTargets whether we ignore read counts that make reference to targets that are not present in
     *                             the input target collection {@code targets}.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     * @throws IOException              if there was any problem reading the content of {@code file}.
     * @throws UserException.BadInput   if there is some formatting issue with the file. This includes inability to
     *                                  resolve a target name based on the input file content and the target collection
     *                                  provided as long as {@code ignoreMissingTargets} is {@code false}.
     */
    public static ReadCountCollection parse(final File file, final TargetCollection<Target> targets,
                                                final boolean ignoreMissingTargets) throws IOException {
        Utils.nonNull(file, "the input file cannot be null");
        final ReadCountsReader reader = new ReadCountsReader(file, targets, ignoreMissingTargets);
        return readCounts(file.getPath(), reader, reader.getCountColumnNames());
    }

    /**
     * Reads the content of a source reader into a {@link ReadCountCollection}.
     * <p>
     * If no target name is included in the input but intervals are present, the {@code exons} collection provided
     * will be utilized to resolve those names.
     * </p>
     *
     * @param sourceReader  the source reader.
     * @param sourceName  source name.
     * @param targets collection of exons (targets). This parameter can be {@code null}, to indicate that no exon
     *              collection is to be considered.
     * @param ignoreMissingTargets whether we ignore read counts that make reference to targets that are not present in
     *                             the input target collection {@code targets}.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     * @throws IOException              if there was any problem reading the content of {@code file}.
     * @throws UserException.BadInput   if there is some formatting issue with the file. This includes inability to
     *                                  resolve a target name based on the input file content and the target collection
     *                                  provided as long as {@code ignoreMissingTargets} is {@code false}.
     */
    public static ReadCountCollection parse(final Reader sourceReader, final String sourceName,
                                                final TargetCollection<Target> targets, final boolean ignoreMissingTargets) throws IOException {
        Utils.nonNull(sourceReader, "the input source reader cannot be null");
        Utils.nonNull(sourceName, "the input source name be null");

        final ReadCountsReader readCountsReader = new ReadCountsReader(sourceReader, targets, ignoreMissingTargets);
        return readCounts(sourceName, readCountsReader, readCountsReader.getCountColumnNames());
    }

    /**
     * Reads the counts section of the file and create the resulting collection.
     *
     * @param sourceName  the source name (used in error messages).
     * @param tableReader the source table-reader.
     * @param columnNames the name of the columns.
     * @return never {@code null}.
     * @throws IOException if there is a low level IO error.
     */
    private static ReadCountCollection readCounts(final String sourceName,
                                                  final TableReader<ReadCountRecord> tableReader,
                                                  final List<String> columnNames) throws IOException {
        final Buffer buffer = new Buffer();

        ReadCountRecord record;
        while ((record = tableReader.readRecord()) != null) {
            final Target target = record.getTarget();
            final double[] lineCounts = record.getDoubleCounts();
            if (!buffer.add(target, lineCounts)) {
                throw new UserException.BadInput(String.format("duplicated target with name %s in %s", target.getName(), sourceName));
            }
        }
        if (buffer.getTargets().isEmpty()) {
            throw new UserException.BadInput("there is no counts (zero targets) in the input source " + sourceName);
        }
        return new ReadCountCollection(buffer.getTargets(), columnNames, new Array2DRowRealMatrix(buffer.getCounts(), false));
    }

    /**
     * Retrieve the sample names from a ReadCountCollection file (e.g. TangentNormalization file)
     * @param readCountsFile targets and coverage
     * @return list of strings that contain the sample names.  I.e. returns list of all columns that do not describe the
     * targets themselves.
     */
    public static List<String> retrieveSampleNamesFromReadCountsFile(final File readCountsFile) {
        try  {
            return new ReadCountsReader(readCountsFile).getCountColumnNames();
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(readCountsFile, e);
        }
    }

    /**
     * Retrieve the sample names from a ReadCountCollection file (e.g. TangentNormalization file)
     * @param readCountsReader targets and coverage reader
     * @param readCountsName name of the read count collection
     * @return list of strings that contain the sample names.  I.e. returns list of all columns that do not describe the
     * targets themselves.
     */
    public static List<String> retrieveSampleNamesFromReadCountsReader(final Reader readCountsReader,
                                                                       final String readCountsName) {
        try  {
            return new ReadCountsReader(readCountsReader).getCountColumnNames();
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(readCountsName, e);
        }
    }

    /** Extract the sample name from a read counts file.
     *
     * Throws an exception if there is not exactly one sample in the file.
     *
     * @param readCountsFile Never {@code null}
     * @return a single sample name.
     */
    public static String getSampleNameForCLIsFromReadCountsFile(final File readCountsFile) {
        String sampleName;
        final List<String> sampleNames = SampleNameFinder.determineSampleNamesFromReadCountsFile(readCountsFile);
        if (sampleNames.size() == 1) {
            sampleName = sampleNames.get(0);
        } else if (sampleNames.size() > 1) {
            throw new UserException.BadInput("Input file must contain data for only one sample.  Found samples: " +
                    StringUtils.join(sampleNames, ", "));
        } else {
            throw new UserException.BadInput("Input file must contain data for only one sample.  Could not find any sample information.");
        }
        return sampleName;
    }

    /** Extract the sample name from a {@link ReadCountCollection}.
     *
     * Throws an exception if there is not exactly one sample in the file.
     *
     * @param counts Never {@code null}
     * @return a single sample name.
     */
    public static String getSampleNameFromReadCounts(final ReadCountCollection counts) {
        Utils.nonNull(counts);
        final List<String> sampleNames = counts.columnNames();
        Utils.validateArg(sampleNames.size() == 1, "read counts must have exactly one sample.");
        return sampleNames.get(0);
    }

    /**
     * Write a read counts file of targets with coverage to file with dummy names
     * @param outFile File to write targets with coverage. Never {@code null}
     * @param sampleName Name of sample being written. Never {@code null}
     * @param byKeySorted Map of simple-intervals to their copy-ratio. Never {@code null}
     * @param comments Comments to add to header of coverage file.
     */
    public static <N extends Number> void writeReadCountsFromSimpleInterval(final File outFile, final String sampleName,
                                                                            final SortedMap<SimpleInterval, N> byKeySorted,
                                                                            final String[] comments) {
        try (final Writer outWriter = new FileWriter(outFile)) {
            writeReadCountsFromSimpleInterval(outWriter, outFile.getPath(), sampleName, byKeySorted, comments);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outFile, e);
        }
    }

    /**
     * Write a read counts file of targets with coverage to file with dummy names
     * @param outWriter Writer to write targets with coverage. Never {@code null}
     * @param outName Name of the output writer. Never {@code null}
     * @param sampleName Name of sample being written. Never {@code null}
     * @param byKeySorted Map of simple-intervals to their copy-ratio. Never {@code null}
     * @param comments Comments to add to header of coverage file.
     */
    public static <N extends Number> void writeReadCountsFromSimpleInterval(final Writer outWriter, final String outName,
                                                                            final String sampleName, final SortedMap<SimpleInterval, N> byKeySorted,
                                                                            final String[] comments) {

        Utils.nonNull(outWriter, "Output writer cannot be null.");
        Utils.nonNull(sampleName, "Sample name cannot be null.");
        Utils.nonNull(byKeySorted, "Targets cannot be null.");
        Utils.nonNull(comments, "Comments cannot be null.");

        final boolean areTargetIntervalsAllPopulated = byKeySorted.keySet().stream().allMatch(t -> t != null);
        if (!areTargetIntervalsAllPopulated) {
            throw new UserException("Cannot write target coverage file with any null intervals.");
        }

        try (final TableWriter<ReadCountRecord> writer =
                     writerWithIntervals(outWriter, Collections.singletonList(sampleName))) {
            for (String comment : comments) {
                writer.writeComment(comment);
            }

            for (final Locatable locatable : byKeySorted.keySet()) {
                final SimpleInterval interval = new SimpleInterval(locatable);
                final double coverage = byKeySorted.get(locatable).doubleValue();
                writer.writeRecord(new ReadCountRecord.SingleSampleRecord(new Target(interval), coverage));
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outName, e);
        }

    }

    /**
     * Impute zero counts to the median of non-zero values in the enclosing target row.
     *
     * <p>The imputation is done in-place, thus the input matrix is well be modified as a result of this call.</p>
     *
     * @param readCounts the input and output read-count matrix.
     */
    public static void imputeZeroCountsAsTargetMedians(final ReadCountCollection readCounts, final Logger logger) {

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
            final double totalZeroCountsPercentage = 100.0 * (totalZeroCounts / totalCounts);
            logger.info(String.format("Some 0.0 counts (%d out of %d, %.2f%%) were imputed to their enclosing target " +
                    "non-zero median", totalZeroCounts, totalZeroCounts, totalZeroCountsPercentage));
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
    public static void truncateExtremeCounts(final ReadCountCollection readCounts, final double percentile, final Logger logger) {

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
            logger.info(String.format("None of the %d counts were truncated as they all fall in the non-extreme range " +
                    "[%.2f, %.2f]", totalCounts, bottomPercentileThreshold, topPercentileThreshold));
        } else {
            final double truncatedPercentage = ((double)(topTruncatedCounts + bottomTruncatedCounts) / totalCounts) * 100;
            logger.info(String.format("Some counts (%d out of %d, %.2f%%) were truncated as they fall out of the " +
                    "non-extreme range [%.2f, %.2f]", topTruncatedCounts + bottomTruncatedCounts, totalCounts,
                    truncatedPercentage, bottomPercentileThreshold, topPercentileThreshold));
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
    public static ReadCountCollection removeColumnsWithExtremeMedianCounts(final ReadCountCollection readCounts,
                                                                           final double extremeColumnMedianCountPercentileThreshold,
                                                                           final Logger logger) {
        final List<String> columnNames = readCounts.columnNames();
        final RealMatrix counts = readCounts.counts();
        final double[] columnMedians = MatrixSummaryUtils.getColumnMedians(counts);

        // Calculate thresholds:
        final double bottomExtremeThreshold = new Percentile(extremeColumnMedianCountPercentileThreshold).evaluate(columnMedians);
        final double topExtremeThreshold = new Percentile(100 - extremeColumnMedianCountPercentileThreshold).evaluate(columnMedians);

        // Determine kept and dropped column sets.
        final Set<String> columnsToKeep = new LinkedHashSet<>(readCounts.columnNames().size());
        final int initialMapSize = ((int) (2.0 * extremeColumnMedianCountPercentileThreshold / 100.0)) * readCounts.columnNames().size();
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
            logger.info(String.format("No column dropped due to extreme counts outside [%.10f, %.10f]",
                    bottomExtremeThreshold, topExtremeThreshold));
            return readCounts;
        } else {
            final double droppedPercentage = ((double)(columnsToDrop.size()) / columnNames.size()) * 100;
            logger.info(String.format("Some columns dropped (%d out of %d, %.2f%%) as they are classified as having extreme " +
                    "median counts across targets outside [%.10f, %.10f]: e.g. %s",
                    columnsToDrop.size(), columnNames.size(), droppedPercentage, bottomExtremeThreshold, topExtremeThreshold,
                    columnsToDrop.entrySet().stream().limit(10).map(kv -> kv.getKey() + " (" + kv.getValue() + ")")
                            .collect(Collectors.joining(", "))));
            return readCounts.subsetColumns(columnsToKeep);
        }
    }

    /**
     * Remove columns with NaNs and infinities
     *
     * @param readCounts
     * @param logger
     * @return
     */
    public static ReadCountCollection removeColumnsWithBadValues(final ReadCountCollection readCounts,
                                                                 final Logger logger) {
        final List<String> columnNames = readCounts.columnNames();

        // Determine kept and dropped column sets.
        final Set<String> columnsToKeep = new LinkedHashSet<>(readCounts.columnNames().size());
        final Set<String> columnsToDrop = new HashSet<>();
        for (int i = 0; i < columnNames.size(); i++) {
            if (Arrays.stream(readCounts.getColumn(i))
                    .filter(d -> Double.isNaN(d) || Double.isInfinite(d))
                    .count() == 0) {
                columnsToKeep.add(columnNames.get(i));
            } else {
                columnsToDrop.add(columnNames.get(i));
            }
        }

        // log and drop columns if it applies
        if (columnsToKeep.isEmpty()) {
            throw new UserException.BadInput("No column count left after removing those with NaN and infinities");
        } else if (columnsToKeep.size() == columnNames.size()) {
            logger.info("No column dropped due to bad values");
            return readCounts;
        } else {
            final double droppedPercentage = ((double)(columnsToDrop.size()) / columnNames.size()) * 100;
            logger.info(String.format("Some columns dropped (%d out of %d, %.2f%%) as they contained bad values: %s",
                    columnsToDrop.size(), columnNames.size(), droppedPercentage,
                    columnsToDrop.stream().collect(Collectors.joining(", ", "[", "]"))));
            return readCounts.subsetColumns(columnsToKeep);
        }
    }

    private static long countZeroes(final double[] data, final boolean roundToInteger) {
        if (roundToInteger) {
            return DoubleStream.of(data).filter(d -> (int) FastMath.round(d) == 0).count();
        } else {
            return DoubleStream.of(data).filter(d -> d == 0.0).count();
        }
    }

    /**
     * Remove targets that have too many counts equal to 0.
     * <p>
     *     It will return a copy of the input read-count collection with such targets dropped.
     * </p>
     *
     * @param readCounts the input read counts.
     * @param maximumTargetZeros maximum number of counts equal to 0 per target tolerated.
     * @return never {@code null}. It might be a reference to the input read-counts if there is
     *   is no target to be dropped.
     */
    public static ReadCountCollection removeTargetsWithTooManyZeros(final ReadCountCollection readCounts,
                                                                    final int maximumTargetZeros, final boolean roundToInteger,
                                                                    final Logger logger) {
        final RealMatrix counts = readCounts.counts();

        final Set<Target> targetsToKeep = IntStream.range(0, counts.getRowDimension()).boxed()
                .filter(i -> countZeroes(counts.getRow(i), roundToInteger) <= maximumTargetZeros)
                .map(i -> readCounts.targets().get(i))
                .collect(Collectors.toCollection(LinkedHashSet::new));

        final int targetsToDropCount = readCounts.targets().size() - targetsToKeep.size();
        if (targetsToDropCount == 0) {
            logger.info(
                    String.format("There are no targets with large number of columns with zero counts (<= %d of %d) to drop",
                            maximumTargetZeros, readCounts.columnNames().size()));
            return readCounts;
        } else if (targetsToDropCount == readCounts.targets().size()) {
            throw new UserException.BadInput("the number of zeros per target in the input is too large resulting " +
                    "in all targets being dropped");
        } else {
            final double droppedPercentage = ((double)(targetsToDropCount) / readCounts.targets().size()) * 100;
            logger.info(String.format("Some targets dropped (%d out of %d, %.2f%%) as they had too many zeros (> %d of %d).",
                    targetsToDropCount, readCounts.targets().size(), droppedPercentage, maximumTargetZeros,
                    readCounts.columnNames().size()));
            return readCounts.subsetTargets(targetsToKeep);
        }
    }

    /**
     * Remove targets that are fully uncovered in the collection
     * <p>
     *     It will return a copy of the input read-count collection with such targets dropped.
     * </p>

     * @param readCounts the input read counts.
     * @param logger instance of logger.
     * @return never {@code null}. It might be a reference to the input read-counts if there is
     *   is no target to be dropped.
     */
    public static ReadCountCollection removeTotallyUncoveredTargets(final ReadCountCollection readCounts, final Logger logger) {
        return removeTargetsWithTooManyZeros(readCounts, readCounts.columnNames().size() - 1, true, logger);
    }

    /**
     * Remove targets that have read counts strictly lower than {@code minReadCount} for all
     * samples.
     * <p>
     *     It will return a copy of the input read-count collection with such targets dropped.
     * </p>
     *
     * TODO github/gatk-protected issue #843 -- write unit test
     *
     * @param readCounts the input read counts.
     * @param minCoverage minimum read count
     * @return never {@code null}. It might be a reference to the input read-counts if there is
     *   is no target to be dropped.
     */
    public static ReadCountCollection removeTargetsWithUniformlyLowCoverage(final ReadCountCollection readCounts,
                                                                            final int minCoverage,
                                                                            final Logger logger) {
        /* filter targets to keep, put them in an ordered set */
        final RealMatrix counts = readCounts.counts();
        final List<Target> originalTargets = readCounts.targets();
        final Set<Target> targetsToKeep = IntStream.range(0, counts.getRowDimension())
                .filter(ti -> Arrays.stream(counts.getRow(ti))
                        .anyMatch(coverage -> coverage >= minCoverage))
                .mapToObj(originalTargets::get)
                .collect(Collectors.toCollection(LinkedHashSet::new));

        final int targetsToDropCount = originalTargets.size() - targetsToKeep.size();
        if (targetsToDropCount == 0) {
            logger.info(String.format("There are no targets with coverage uniformly lower than %d. Keeping all" +
                    " %d targets.", minCoverage, originalTargets.size()));
            return readCounts;
        } else if (targetsToDropCount == readCounts.targets().size()) {
            throw new UserException.BadInput(String.format("The coverage on all targets and samples lower than %d," +
                    " resulting in all targets being dropped", minCoverage));
        } else {
            final double droppedPercentage = ((double)(targetsToDropCount) / originalTargets.size()) * 100;
            logger.info(String.format("Some targets dropped (%d out of %d, %.2f%%) as their coverage were uniformly" +
                    " lower than %d.", targetsToDropCount, originalTargets.size(), droppedPercentage, minCoverage));
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
     * @param maximumColumnZeros maximum number of counts equal to 0 per column tolerated.
     * @return never {@code null}. It might be a reference to the input read-counts if there is
     *   is no column to be dropped.
     */
    @VisibleForTesting
    public static ReadCountCollection removeColumnsWithTooManyZeros(final ReadCountCollection readCounts,
                                                                    final int maximumColumnZeros,
                                                                    final boolean roundToInteger,
                                                                    final Logger logger) {

        final RealMatrix counts = readCounts.counts();

        final Set<String> columnsToKeep = IntStream.range(0, counts.getColumnDimension()).boxed()
                .filter(i -> countZeroes(counts.getColumn(i), roundToInteger) <= maximumColumnZeros)
                .map(i -> readCounts.columnNames().get(i))
                .collect(Collectors.toCollection(LinkedHashSet::new));

        final int columnsToDropCount = readCounts.columnNames().size() - columnsToKeep.size();
        if (columnsToDropCount == 0) {
            logger.info(
                    String.format("There were no columns with a large number of targets with zero counts " +
                            "(<= %d of %d) to drop", maximumColumnZeros, readCounts.targets().size()));
            return readCounts;
        } else if (columnsToDropCount == readCounts.columnNames().size()) {
            throw new UserException.BadInput("The number of zeros per count column is too large resulting in all count " +
                    "columns to be dropped");
        } else {
            final double droppedPercentage = ((double)(columnsToDropCount) / readCounts.columnNames().size()) * 100;
            logger.info(String.format("Some counts columns dropped (%d out of %d, %.2f%%) as they had too many targets with zeros (> %d of %d)",
                    columnsToDropCount, readCounts.columnNames().size(), droppedPercentage, maximumColumnZeros,
                    readCounts.targets().size()));
            return readCounts.subsetColumns(columnsToKeep);
        }
    }

    /**
     * Helper class used to accumulate read counts, target names and intervals as they are read
     * from the source file.
     * <p>
     * Its capacity auto-extends as more targets are added to it.
     * </p>
     */
    private static class Buffer {

        /**
         * Set of targets indexed by their names.
         * <p>
         * The correspondence between targets and columns in {@code counts} is through this map iteration order.
         * </p>
         */
        private SetUniqueList<Target> targets;

        /**
         * Contains the counts so far.
         */
        private List<double[]> counts;

        /**
         * Creates a new buffer
         */
        private Buffer() {
            targets = SetUniqueList.setUniqueList(new ArrayList<>());
            counts = new ArrayList<>();
        }

        /**
         * Adds a new target and counts to the buffer.
         * <p>This call will do anything if the target is already in the buffer returning {@code false}.</p>
         *
         * @param target the target.
         * @param values the counts for that target.
         * @return true iff {@code target} was new to the buffer, {@code false} otherwise.
         */
        private boolean add(final Target target, final double[] values) {
            if (targets.add(target)) {
                counts.add(values);
                return true;
            } else {
                return false;
            }
        }

        /**
         * Returns a live modifiable unique list to the targets already in the buffer.
         *
         * @return never {@code null}.
         */
        private SetUniqueList<Target> getTargets() {
            return targets;
        }

        /**
         * Returns an array representation of the counts in the buffer.
         * <p>Each element of array corresponds the the ith target in the buffer.</p>
         * <p>The result array can be modified at will without altering the buffer, yet the element sub-arrays are
         * live objects and modifications will change the counts in the buffer</p>
         */
        private double[][] getCounts() {
            return counts.toArray(new double[counts.size()][]);
        }
    }
}
