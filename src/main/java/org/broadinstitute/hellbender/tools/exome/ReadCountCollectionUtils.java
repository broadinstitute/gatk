package org.broadinstitute.hellbender.tools.exome;

import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Reads {@link ReadCountCollection} instances from a tab-separated text file.
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
 */
public final class ReadCountCollectionUtils {

    // Prevents instantiation of the class.
    private ReadCountCollectionUtils() {}

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
        Utils.nonNull(collection, "input collection cannot be null");
        Utils.nonNull(headerComments, "header comments cannot be null");
        try (final Writer writer = new FileWriter(file)) {
            final boolean withIntervals = collection.targets().stream().anyMatch(t -> t.getInterval() != null);
            final TableWriter<ReadCountRecord> tableWriter = withIntervals
                    ? writerWithIntervals(writer, collection.columnNames()) : writerWithoutIntervals(writer, collection);
            performWriting(collection, tableWriter, headerComments);
        }
    }

    /**
     * Writes the content of a collection into a file as if it was target coverage.  Intervals will always be included,
     *  even if non-existent
     *
     * @param outFile           the output file.
     * @param collection     the output collection.  MUST only contain one sample.
     * @param headerComments header comments.
     * @throws IllegalArgumentException if any of the input parameters is {@code null}
     *                                  or {@code collection} has no intervals
     *                                  defined or {@code collection} has more than one sample.
     * @throws IOException              if there is some IO issue when writing into the output file.
     *
     */
    public static void writeAsTargetCoverage(final File outFile, final ReadCountCollection collection, final String... headerComments) throws IOException {
        Utils.nonNull(outFile, "output file cannot be null");
        Utils.nonNull(collection, "input collection cannot be null");
        Utils.nonNull(headerComments, "header comments cannot be null.  Use 'new String [0]', if you would like no comments.");
        if (collection.columnNames().size() != 1) {
            throw new UserException.BadInput("Attempting to write a read count collection with more than one sample as a target coverage.  This is currently not supported.  Please use an input with one sample.");
        }

        final List<TargetCoverage> targetCollection = convertToTargetCoverageList(collection);
        TargetCoverageUtils.writeTargetsWithCoverage(outFile, collection.columnNames().get(0), targetCollection, headerComments);
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

        columnNames.add(TargetTableColumns.CONTIG.toString());
        columnNames.add(TargetTableColumns.START.toString());
        columnNames.add(TargetTableColumns.END.toString());
        columnNames.add(TargetTableColumns.NAME.toString());
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
        columnNames.add(TargetTableColumns.NAME.toString());
        columnNames.addAll(collection.columnNames());
        return createReadCountRecordTableWriterWithoutIntervals(writer, columnNames);
    }

    /**
     *  Convert a ReadCountCollection into a List of TargetCoverage.
     *
     *  <p>Caveats, please read: </p>
     *  <ul>
     *      <li>A list of TargetCoverage implies that only one sample is supported.  An exception will be thrown if this method
     *      receives more than one sample in the ReadCountCollection</li>
     *      <li>TargetCoverage does not currently support sample name.  Therefore, this information will be lost upon conversion.</li>
     *  </ul>
     *
     * @param collection -- ReadCountCollection that contains one sample only
     * @return never {@code null}.
     */
    private static List<TargetCoverage> convertToTargetCoverageList(final ReadCountCollection collection) {

        Utils.nonNull(collection, "Cannot convert a null ReadCountCollection to TargetCoverage.");

        if (collection.counts().getColumnDimension() != 1) {
            throw new UserException("Cannot convert ReadCountCollection with multiple samples to TargetCoverageCollection.  Please use input containing one sample.");
        }

        // Please note that the "0" in getEntry(i,0) is the assumption of one sample only
        return IntStream.range(0, collection.targets().size())
                .mapToObj(i -> new TargetCoverage(collection.targets().get(i).getName(), collection.targets().get(i).getInterval(), collection.counts().getEntry(i,0)))
                .collect(Collectors.toList());
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
     * If no target name is include in the input but intervals are present, the {@code exons} collection provided
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
     *                                  resolve a target name based on the input file content and the target collection provided as long as
     *                                  {@code ignoreMissingTargets} is {@code false}.
     */
    public static <E> ReadCountCollection parse(final File file, final TargetCollection<E> targets, final boolean ignoreMissingTargets) throws IOException {
        Utils.nonNull(file, "the input file cannot be null");

        final List<String> countColumnNames = new ArrayList<>();

        final TableReader<ReadCountRecord> tableReader = new TableReader<ReadCountRecord>(file) {

            private Function<DataLine, ReadCountRecord> recordExtractor;

            @Override
            protected void processColumns(final TableColumnCollection columns) {
                countColumnNames.clear();
                countColumnNames.addAll(columns.names().stream()
                        .filter(name -> !TargetTableColumns.isStandardTargetColumnName(name))
                        .collect(Collectors.toList()));

                @SuppressWarnings("all")
                final Function<DataLine, SimpleInterval> intervalExtractor = intervalExtractor(columns, (message) -> formatException(message));
                final Function<DataLine, String> targetNameExtractor = targetNameExtractor(columns);
                final Function<DataLine, double[]> countExtractor = countExtractor(columns);
                recordExtractor = composeRecordExtractor(intervalExtractor, targetNameExtractor, countExtractor, targets);
            }

            private Function<DataLine, ReadCountRecord> composeRecordExtractor(final Function<DataLine, SimpleInterval> intervalExtractor,
                        final Function<DataLine, String> targetNameExtractor, final Function<DataLine, double[]> countExtractor,
                        final TargetCollection<E> targets) {
                if (targetNameExtractor == null && (targets == null || intervalExtractor == null)) {
                    throw formatException("the input files does not contain a target name column and no target file was provided");
                } else if (targetNameExtractor != null && intervalExtractor != null) {
                    if (targets == null) {
                        return (v) -> new ReadCountRecord(new Target(targetNameExtractor.apply(v), intervalExtractor.apply(v)), countExtractor.apply(v));
                    } else {
                        return (v) -> extractRecordWithTargetNameAndIntervalExtractorsAndTargetCollection(v, intervalExtractor, targetNameExtractor, countExtractor);
                    }
                } else if (targetNameExtractor != null) {
                    if (targets == null) {
                        return (v) -> new ReadCountRecord(new Target(targetNameExtractor.apply(v), null), countExtractor.apply(v));
                    } else {
                        return (v) -> extractRecordWithTargetNameExtractorAndTargetCollection(v, targetNameExtractor, countExtractor);
                    }
                } else { // at this point targets != null && intervalExtractor != null.
                    return (v) -> extractRecordWithIntervalExtractorAndTargetCollection(v, intervalExtractor, countExtractor);
                }
            }

            private ReadCountRecord extractRecordWithTargetNameExtractorAndTargetCollection(final DataLine v,
                    final Function<DataLine, String> targetNameExtractor, final Function<DataLine, double[]> countExtractor) {
                final String name = targetNameExtractor.apply(v);
                final E nameTarget = targets.target(name);
                if (nameTarget == null) {
                    if (ignoreMissingTargets) {
                        return null;
                    } else {
                        throw formatException(String.format("unknown target '%s' not present in the target collection", name));
                    }
                } else {
                    return new ReadCountRecord(new Target(name, targets.location(nameTarget)), countExtractor.apply(v));
                }
            }

            private ReadCountRecord extractRecordWithIntervalExtractorAndTargetCollection(final DataLine v,
                    final Function<DataLine, SimpleInterval> intervalExtractor, final Function<DataLine, double[]> countExtractor) {
                final SimpleInterval interval = intervalExtractor.apply(v);
                final E intervalTarget = targets.target(interval);
                if (intervalTarget == null) {
                    if (ignoreMissingTargets) {
                        return null;
                    } else {
                        throw formatException(String.format("unknown target with interval %s", interval));
                    }
                } else if (!targets.location(intervalTarget).equals(interval)) {
                    throw formatException(String.format("mismatching yet overlapping intervals in the input (%s) and the target collection (%s)", interval, targets.location(intervalTarget)));
                } else {
                    return new ReadCountRecord(new Target(targets.name(intervalTarget), interval), countExtractor.apply(v));
                }
            }

            private ReadCountRecord extractRecordWithTargetNameAndIntervalExtractorsAndTargetCollection(final DataLine data,
                        final Function<DataLine, SimpleInterval> intervalExtractor, final Function<DataLine, String> targetNameExtractor,
                        final Function<DataLine, double[]> countExtractor) {
                    final String name = targetNameExtractor.apply(data);
                    final SimpleInterval interval = intervalExtractor.apply(data);
                    final E nameTarget = targets.target(name);
                    final E intervalTarget = targets.target(interval);
                    if (!Objects.equals(nameTarget, intervalTarget)) {
                        throw formatException(String.format("conflicting target resolution from the name (%s) and interval (%s) provided", name, interval));
                    } else if (nameTarget == null) {
                        return ignoreMissingTargets ? null : new ReadCountRecord(new Target(name, interval), countExtractor.apply(data));
                    } else if (!targets.location(intervalTarget).equals(interval)) {
                        throw formatException(String.format("mismatch interval from input (%s) and target collection (%s) for target '%s'", interval, targets.location(intervalTarget), name));
                    } else {
                        return new ReadCountRecord(new Target(name, interval), countExtractor.apply(data));
                    }
            }

            @Override
            protected ReadCountRecord createRecord(final DataLine dataLine) {
                return recordExtractor.apply(dataLine);
            }

            private Function<DataLine, double[]> countExtractor(final TableColumnCollection columns) {
                final int[] countColumnIndexes = IntStream.range(0, columns.columnCount())
                        .filter(i -> !TargetTableColumns.isStandardTargetColumnName(columns.nameAt(i))).toArray();
                return (v) -> {
                    final double[] result = new double[countColumnIndexes.length];
                    for (int i = 0; i < countColumnIndexes.length; i++) {
                        result[i] = v.getDouble(countColumnIndexes[i]);
                    }
                    return result;
                };
            }
        };

        return readCounts(file, tableReader, countColumnNames);

    }

    /**
     * Reads the counts section of the file and create the resulting collection.
     *
     * @param file        the source file name (used in error messages).
     * @param tableReader the source table-reader.
     * @param columnNames the name of the columns.
     * @return never {@code null}.
     * @throws IOException if there is a low level IO error.
     */
    private static ReadCountCollection readCounts(final File file,
                                                  final TableReader<ReadCountRecord> tableReader,
                                                  final List<String> columnNames) throws IOException {
        final Buffer buffer = new Buffer();

        ReadCountRecord record;
        while ((record = tableReader.readRecord()) != null) {
            final Target target = record.getTarget();
            final double[] lineCounts = record.getDoubleCounts();
            if (!buffer.add(target, lineCounts)) {
                throw new UserException.BadInput(String.format("duplicated target with name %s in file %s", target.getName(), file));
            }
        }
        if (buffer.getTargets().size() == 0) {
            throw new UserException.BadInput("there is no counts (zero targets) in the input file " + file);
        }
        return new ReadCountCollection(buffer.getTargets(), SetUniqueList.setUniqueList(columnNames),
                new Array2DRowRealMatrix(buffer.getCounts(),false));
    }

    /**
     * Composes a lambda to extract the name of the target given a row of values from the input read-count file.
     * <p>
     * This method will return {@code null} if it is not possible to extract the target name from the input directly; for
     * example the input only contain the coordinates of the target and not the target name itself (
     * (i.e. the {@link TargetTableColumns#NAME NAME} column is missing).
     * </p>
     *
     * @param columns the column-name array for that file.
     * @return non-{@code null} iff is not possible to extract the target name from the input directly.
     */
    private static Function<DataLine, String> targetNameExtractor(final TableColumnCollection columns) {
        final int nameColumnIndex = columns.indexOf(TargetTableColumns.NAME.toString());
        return nameColumnIndex < 0 ? null : (v) -> v.get(nameColumnIndex);
    }

    /**
     * Constructs an per line interval extractor given the header column names.
     *
     * @param columns               the header column names.
     * @param errorExceptionFactory the error handler to be called when there is any problem resoling the interval.
     * @return never {@code null} if there is enough columns to extract the coordinate information, {@code null} otherwise.
     */
    private static Function<DataLine, SimpleInterval> intervalExtractor(final TableColumnCollection columns,
                                                                        final Function<String, RuntimeException> errorExceptionFactory) {

        final int contigColumnNumber = columns.indexOf(TargetTableColumns.CONTIG.toString());
        final int startColumnNumber = columns.indexOf(TargetTableColumns.START.toString());
        final int endColumnNumber = columns.indexOf(TargetTableColumns.END.toString());
        return composeIntervalBuilder(contigColumnNumber, startColumnNumber, endColumnNumber, errorExceptionFactory);
    }

    /**
     * Returns a function that translate an source line string value array into
     * into a interval.
     *
     * @param contigColumnNumber    the number of the input column that contains the
     *                              contig name. {@code -1} if missing.
     * @param startColumnNumber     the number of the input column that contains the
     *                              start position. {@code -1} if missing.
     * @param endColumnNumber       the number of the input column that contains the
     *                              end position. {@code -1} if missing.
     * @param errorExceptionFactory instantiates the exception to thrown in case
     *                              of a formatting error. cannot be {@code null}.
     * @return not a {@code null} if there is enough information to find out the
     * sample intervals, {@code null} if it is insufficient.
     */
    private static Function<DataLine, SimpleInterval> composeIntervalBuilder(final int contigColumnNumber,
                                                                             final int startColumnNumber,
                                                                             final int endColumnNumber,
                                                                             final Function<String, RuntimeException> errorExceptionFactory) {
        if (contigColumnNumber == -1 || startColumnNumber == -1 || endColumnNumber == -1) {
            return null;
        }

        return (v) -> {
            final String contig = v.get(contigColumnNumber);
            final int start = v.getInt(startColumnNumber);
            final int end = v.getInt(endColumnNumber);
            if (start <= 0) {
                throw errorExceptionFactory.apply(String.format("start position must be greater than 0: %d", start));
            } else if (start > end) {
                throw errorExceptionFactory.apply(String.format("end position '%d' must equal or greater than the start position '%d'", end, start));
            } else {
                return new SimpleInterval(contig, start, end);
            }
        };

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
            this.targets = SetUniqueList.setUniqueList(new ArrayList<>());
            this.counts = new ArrayList<>();
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
