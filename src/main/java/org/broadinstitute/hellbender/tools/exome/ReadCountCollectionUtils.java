package org.broadinstitute.hellbender.tools.exome;

import htsjdk.samtools.util.Locatable;
import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.exome.samplenamefinder.SampleNameFinder;
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
import java.util.Arrays;
import java.util.List;
import java.util.SortedMap;

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
    public static ReadCountCollection parse(final File file, final TargetCollection<Target> targets, final boolean ignoreMissingTargets) throws IOException {
        Utils.nonNull(file, "the input file cannot be null");
        final ReadCountsReader reader = new ReadCountsReader(file, targets, ignoreMissingTargets);
        return readCounts(file, reader, reader.getCountColumnNames());
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
        return new ReadCountCollection(buffer.getTargets(), columnNames, new Array2DRowRealMatrix(buffer.getCounts(),false));
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
            throw new UserException.BadInput("Input file must contain data for only one sample.  Found samples: " + StringUtils.join(sampleNames, ", "));
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

        Utils.nonNull(outFile, "Output file cannot be null.");
        Utils.nonNull(sampleName, "Sample name cannot be null.");
        Utils.nonNull(byKeySorted, "Targets cannot be null.");
        Utils.nonNull(comments, "Comments cannot be null.");

        final boolean areTargetIntervalsAllPopulated = byKeySorted.keySet().stream().allMatch(t -> t != null);
        if (!areTargetIntervalsAllPopulated) {
            throw new UserException("Cannot write target coverage file with any null intervals.");
        }

        try (final TableWriter<ReadCountRecord> writer =
                     writerWithIntervals(new FileWriter(outFile), Arrays.asList(sampleName))) {
            for (String comment : comments) {
                writer.writeComment(comment);
            }

            for (final Locatable locatable : byKeySorted.keySet()) {
                final SimpleInterval interval = new SimpleInterval(locatable);
                final double coverage = byKeySorted.get(locatable).doubleValue();
                writer.writeRecord(new ReadCountRecord.SingleSampleRecord(new Target(interval), coverage));
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outFile, e);
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
