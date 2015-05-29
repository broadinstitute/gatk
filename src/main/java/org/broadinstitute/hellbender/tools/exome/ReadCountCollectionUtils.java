package org.broadinstitute.hellbender.tools.exome;

import htsjdk.tribble.bed.BEDFeature;
import org.apache.commons.collections4.list.SetUniqueList;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.function.Function;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Reads {@link ReadCountCollection} instances from a tab-separated text file.
 *
 * <p>
 *     The tab separated file consist of a header and body with the data.
 * </p>
 *
 * <p>
 *     The header consist of at least a line with column names optionally preceded by comment lines (starting with '#').
 *     A part from target coordinates and name columns there should be at least on actual count column (sample, read-group or cohort).
 *     but there could be more than one.
 * </p>
 *
 * <p>
 *     The body are the coordinates and counts for each target.
 * </p>
 * <p>
 *     Example:
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
 *
 * <p>
 *     You may omit either the target name column (NAME) or some of the genomic interval columns (CONTIG, START and END)
 *     but not both at the same time.
 * </p>
 *
 * <p>
 *     If the source omits the target name, a exonCollection should be provided in order to resovle the name given its coordinates
 *     using {@link #parse(File, ExonCollection)}.
 * </p>
 *
 * <p>
 *     This class will check whether the content of the input file is well formatted and consistent
 *     (e.g. counts are double values, each row have the same number of values, on for each column in the header,
 *      and so forth).
 * </p>
 *
 * <p>
 *     If there is any formatting problems the appropriate exception will be thrown
 *     as described in {@link #parse}.
 * </p>
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class ReadCountCollectionUtils {

    public static final String COMMENT_PREFIX = "#";
    public static final String COLUMN_SEPARATOR = "\t";

    // Prevents instantiation of the class.
    private ReadCountCollectionUtils() {}

    /**
     * Writes the content of a collection into a file.
     * @param file the output file.
     * @param collection the output collection.
     * @param headerComments header comments.
     *
     * @throws IllegalArgumentException if any of the input parameters is {@code null}
     * or {@code collection} has a mixture of targets with and without intervals
     * defined.
     *
     * @throws IOException if there is some IO issue when writing into the output file.
     */
    public static void write(final File file, final ReadCountCollection collection, final String ... headerComments) throws IOException {
        Utils.nonNull(file,"output file cannot be null");
        Utils.nonNull(collection,"input collection cannot be null");
        Utils.nonNull(headerComments, "header comments cannot be null");
        try (final PrintWriter writer = new PrintWriter(new FileWriter(file))) {
            // print the header comments
            for (final String comment : headerComments) {
                writer.println(COMMENT_PREFIX + Utils.nonNull(comment, "header comments cannot contain nulls"));
            }
            if (collection.targets().stream().anyMatch(t -> t.getInterval() != null)) {
                writeContentWithIntervals(writer, collection);
            } else {
                writeContentWithoutIntervals(writer, collection);
            }
        }
    }

    /**
     * Prints column header and count content when without target intervals.
     * @param writer the output writer.
     * @param collection the input collection.
     * @throws IOException if there is some IO issue when writing into the output file.
     */
    private static void writeContentWithoutIntervals(final PrintWriter writer, final ReadCountCollection collection)
            throws IOException {
        final String countColumnString = String.join(COLUMN_SEPARATOR,collection.columnNames());
        writer.println(String.join(COLUMN_SEPARATOR, ReadCountsSpecialColumns.NAME.name(), countColumnString));
        final List<Target> targets = collection.targets();
        final RealMatrix counts = collection.counts();
        for (int i = 0; i < targets.size(); i++) {
            final Target target = targets.get(i);
            final double[] values = counts.getRow(i);
            writer.println(String.join(COLUMN_SEPARATOR,target.getName(),Utils.join(COLUMN_SEPARATOR,values)));
        }
    }

    /**
     * Prints column header and count content when with target intervals.
     * @param writer the output writer.
     * @param collection the input collection.
     * @throws IllegalStateException if any target in the collection does not have an interval defined.
     * @throws IOException if there is some IO issue when writing into the output file.
     */
    private static void writeContentWithIntervals(final PrintWriter writer, final ReadCountCollection collection)
            throws IOException {
        final String countColumnString = String.join(COLUMN_SEPARATOR,collection.columnNames());
        writer.println(String.join(COLUMN_SEPARATOR, ReadCountsSpecialColumns.CONTIG.name(),
                ReadCountsSpecialColumns.START.name(),
                ReadCountsSpecialColumns.END.name(),
                ReadCountsSpecialColumns.NAME.name(),
                countColumnString));
        final List<Target> targets = collection.targets();
        final RealMatrix counts = collection.counts();
        for (int i = 0; i < targets.size(); i++) {
            final Target target = targets.get(i);
            final SimpleInterval interval = target.getInterval();
            if (interval == null) {
                throw new IllegalStateException("invalid combination of targets with and without intervals defined");
            }
            final double[] values = counts.getRow(i);
            writer.println(String.join(COLUMN_SEPARATOR,interval.getContig(),
                    Integer.toString(interval.getStart()),
                    Integer.toString(interval.getEnd()),
                    target.getName(),
                    Utils.join(COLUMN_SEPARATOR,values)));
        }
    }

    /**
     * Reads the content of a file into a {@link ReadCountCollection}.
     *
     * @param file the source file.
     * @return never {@code null}.
     * @throws IOException if there was some problem reading the file contents.
     * @throws UserException.BadInput if there is some formatting issue win the source file contents. This includes
     *    lack of target names in the source file.
     */
    public static ReadCountCollection parse(final File file) throws IOException {
        return parse(file,null);
    }

    /**
     * Reads the content of a file into a {@link ReadCountCollection}.
     *
     * <p>
     *     If no target name is include in the input but intervals are present, the {@code exons} collection provided
     *     will be utilized to resolve those names.
     * </p>
     *
     * @param file the source file.
     * @param exons collection of exons (targets). This parameter can be {@code null}, to indicate that no exon
     *              collection is to be considered.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code file} is {@code null}.
     * @throws IOException if there was any problem reading the content of {@code file}.
     * @throws UserException.BadInput if there is some formatting issue with the file. This includes inability to
     *   resolve a target name based on the input file content and the exon collection provided.
     */
    public static <E extends BEDFeature> ReadCountCollection parse(final File file, final ExonCollection<E> exons) throws IOException {
        Utils.nonNull(file, "the input file cannot be null");
        try (final LineNumberReader reader =
                     new LineNumberReader(new FileReader(file))) {
            final Function<String,RuntimeException> errorExceptionFactory = (message) -> new UserException.BadInput(locateErrorMessage(file,reader,message));
            final String headerColumnNamesString = readHeader(file, reader, errorExceptionFactory);
            final String[] headerColumnNames = headerColumnNamesString.split(COLUMN_SEPARATOR);
            if (Stream.of(headerColumnNames).distinct().count() != headerColumnNames.length) {
                throw errorExceptionFactory.apply("there are repeated count column names");
            }
            final SetUniqueList<String> countColumnNames = extractCountColumnNames(headerColumnNames);
            if (countColumnNames.size() == 0) {
                throw errorExceptionFactory.apply("there is no count columns after removing all special column names: " + ReadCountsSpecialColumns.allSpecialColumnNameString);
            }

            final Function<String[], SimpleInterval> intervalExtractor = intervalExtractor(headerColumnNames,
                    (message) -> new UserException.BadInput(locateErrorMessage(file, reader, message)));
            final Function<String[], String> targetNameExtractor = targetNameExtractor(headerColumnNames);
            if (targetNameExtractor == null && exons == null) {
                throw new UserException.BadInput(locateErrorMessage(file, reader, "the input files does not contain a target name column and no target file was provided"));
            }

            final Function<String[], Target> targetExtractor = targetExtractor(exons, targetNameExtractor, intervalExtractor, errorExceptionFactory);
            if (intervalExtractor == null && targetNameExtractor == null) {
                throw errorExceptionFactory.apply("missing special header column names: " + ReadCountsSpecialColumns.allSpecialColumnNameString);
            }
            final Function<String[], double[]> countExtractor = countValuesExtractor(headerColumnNames, errorExceptionFactory);
            return readCounts(file, reader, headerColumnNames.length, countColumnNames, targetExtractor, countExtractor);
        }
    }

    /**
     * Composes a lambda to extract the target information from the input row.
     *
     * @param exons the exon collection to default to if there is no target names but coordinates in the input.
     * @param targetNameExtractor the target name extractor.
     * @param intervalExtractor the genomic interval extractor.
     * @param exceptionFactory factory to use in order to create the exception to thrown in case of an format error.
     * @return never {@code null}.
     */
    private static <E extends BEDFeature> Function<String[],Target> targetExtractor(final ExonCollection<E> exons,
                                                             final Function<String[], String> targetNameExtractor,
                                                             final Function<String[], SimpleInterval> intervalExtractor,
                                                             final Function<String,RuntimeException> exceptionFactory) {
        return (values) -> {
                final SimpleInterval interval = intervalExtractor == null ? null : intervalExtractor.apply(values);
                final String name = targetNameExtractor == null ? exons.name(exons.exon(interval)) : targetNameExtractor.apply(values);
                if (name == null) {
                    throw exceptionFactory.apply("cannot resolve the target name for interval " + interval);
                }
                return new Target(name,interval);
        };
    }

    /**
     * Reads the counts section of the file and create the resulting collection.
     *
     * @param file the source file name (used in error messages).
     * @param reader the source reader.
     * @param totalColumnCount total number of columns.
     * @param countColumnNames name of count containing columns.
     * @param targetExtractor lambda to extract interval information from the parsed and split line.
     * @param countExtractor lambda to extract counts from the parsed and split line.
     *
     * @return never {@code null}.
     * @throws IOException if there is a low level IO error.
     */
    private static ReadCountCollection readCounts(final File file, final LineNumberReader reader,
                                           final int totalColumnCount,
                                           final SetUniqueList<String> countColumnNames,
                                           final Function<String[], Target> targetExtractor,
                                           final Function<String[], double[]> countExtractor) throws IOException {
        final Buffer buffer = new Buffer();

        String line;
        while ((line = reader.readLine()) != null) {
            String[] values = line.split(COLUMN_SEPARATOR);
            if (values.length != totalColumnCount) {
                throw new UserException.BadInput(locateErrorMessage(file, reader,
                        String.format("number of elements '%d' does not match the number of header columns '%d'",
                                values.length, totalColumnCount)));
            }
            final Target target = targetExtractor.apply(values);
            final double[] lineCounts = countExtractor.apply(values);
            if (!buffer.add(target,lineCounts)) {
                throw new UserException.BadInput(locateErrorMessage(file, reader,
                        String.format("duplicated target with name %s",target.getName())));
            }
        }
        if (buffer.getTargets().size() == 0) {
            throw new UserException.BadInput(locateErrorMessage(file, reader, "there is no counts (zero targets) in the input file"));
        }
        return new ReadCountCollection(buffer.getTargets(), countColumnNames, new Array2DRowRealMatrix(buffer.getCounts()));
    }

    /**
     * Reads out the commented out header and returns the line that contains the column names.
     *
     * @param file the source file.
     * @param reader the reader to use in order to read lines from the source file.
     * @param errorExceptionFactory to use to instantiate the exception to throw in case there is a formatting error.
     * @return never {@code null}.
     * @throws IOException if there is a IO error accessing the contents for the source {@code file} throw the input {@code reader}.
     * @throws RuntimeException if there is a formatting error, the exception type is determine by the factory lambda {@code errorExceptionFactory}.
     */
    private static String readHeader(final File file, final LineNumberReader reader, final Function<String,RuntimeException> errorExceptionFactory) throws IOException {
        String line;
        while ((line = reader.readLine()) != null) {
            if (!line.startsWith("#")) {
                break;
            }
        }
        if (line == null) {
            throw errorExceptionFactory.apply(String.format("file '%s': premature end of file column names",file));
        }
        return line;
    }

    /**
     * Prefixes location information to an error message.
     *
     * <p>
     *     It adds the file name and the line number.
     * </p>
     *
     * @param file the source file.
     * @param reader the reader used to read the source file.
     * @param message the message to extend.
     * @return never {@code null}.
     */
    private static String locateErrorMessage(final File file, final LineNumberReader reader, final String message) {
        return String.format("file '%s' line '%d': %s",file,reader.getLineNumber(),message);
    }

    /**
     * Returns the count column names amongst an input array with all column names in the header.
     *
     * <p>
     *     Basically any column name that is not an special one in this enum is considered a count
     *     column.
     * </p>
     *
     * <p>
     *     The result contains the count column names in the same order as they appear in the
     *     input array.
     * </p>
     * @param headerColumnNames all column names.
     * @return never {@code null}
     */
    private static SetUniqueList<String> extractCountColumnNames(final String[] headerColumnNames) {
        final SetUniqueList<String> result = SetUniqueList.setUniqueList(new ArrayList<>());
        for (final String columnName : headerColumnNames) {
            if (!ReadCountsSpecialColumns.isSpecialColumnName(columnName)) {
                result.add(columnName);
            }
        }
        return result;
    }

    /**
     * Creates a lambda function that extracts the counts value given header column names.
     * <p>
     *     It assumes that any non-special column (as returned by {@link #extractCountColumnNames(String[])}),
     *     is a count column. Thus the resulting extractor
     *     will parse the values at those columns into a double array.
     * </p>
     *
     * @param headerColumnNames  header column names.
     * @param errorExceptionFactory called when there seem to be an error in the input data; an human friendly explanation
     *                     of the issue is provided as an argument.
     * @throws RuntimeException if any is thrown by the {@code errorHandler} provided.
     * @return never {@code null}.
     */
    private static Function<String[], double[]> countValuesExtractor(final String[] headerColumnNames,
                                                                                           final Function<String,RuntimeException> errorExceptionFactory)  {
        final int[] countColumnIndexes = IntStream.range(0, headerColumnNames.length)
                .filter(i -> !ReadCountsSpecialColumns.isSpecialColumnName(headerColumnNames[i])).toArray();
        return (v) -> {
            final double[] result = new double[countColumnIndexes.length];
            for (int i = 0; i < countColumnIndexes.length; i++) {
                try {
                    result[i] = Double.parseDouble(v[countColumnIndexes[i]]);
                } catch (final NumberFormatException ex) {
                    throw errorExceptionFactory.apply(String.format("element %d '%s' cannot be parsed into a double value", countColumnIndexes[i], v[countColumnIndexes[i]]));
                }
            }
            return result;
        };
    }

    /**
     * Composes a lambda to extract the name of the target given a row of values from the input read-count file.
     *
     * <p>
     *     This method will return {@code null} if it is not possible to extract the target name from the input directly; for
     *     example the input only contain the coordinates of the target and not the target name itself (
     *     (i.e. the {@link ReadCountsSpecialColumns#NAME NAME} column is missing).
     * </p>
     *
     * @param headerNames the column-name array for that file.
     * @return non-{@code null} iff is not possible to extract the target name from the input directly.
     */
    private static Function<String[], String> targetNameExtractor(final String[] headerNames) {
        for (int i = 0; i < headerNames.length; i++) {
            final String columnName = headerNames[i];
            if (ReadCountsSpecialColumns.NAME.name().equals(columnName)) {
                final int nameIndex = i;
                return (v) -> v[nameIndex];
            }
        }
        return null;
    }

    /**
     * Constructs an per line interval extractor given the header column names.
     *
     * @param headerColumnNames the header column names.
     * @param errorExceptionFactory the error handler to be called when there is any problem resoling the interval.
     * @return never {@code null} if there is enough columns to extract the coordinate information, {@code null} otherwise.
     */
    private static Function<String[], SimpleInterval> intervalExtractor(final String[] headerColumnNames,
                                                                        final Function<String,RuntimeException> errorExceptionFactory) {

        int contigColumnNumber = -1;
        int startColumnNumber = -1;
        int endColumnNumber = -1;
        for (int i = 0; i < headerColumnNames.length; i++) {
            final String columnName = headerColumnNames[i];
            if (ReadCountsSpecialColumns.isSpecialColumnName(columnName)) {
                switch (ReadCountsSpecialColumns.valueOf(columnName)) {
                    case CONTIG:
                        contigColumnNumber = i;
                        break;
                    case START:
                        startColumnNumber = i;
                        break;
                    case END:
                        endColumnNumber = i;
                        break;
                    default:
                }
            }
        }
        return composeIntervalBuilder(contigColumnNumber, startColumnNumber, endColumnNumber, errorExceptionFactory);
    }

    /**
     * Returns a function that translate an source line string value array into
     * into a interval.
     *
     * @param contigColumnNumber the number of the input column that contains the
     *                           contig name. {@code -1} if missing.
     * @param startColumnNumber  the number of the input column that contains the
     *                           start position. {@code -1} if missing.
     * @param endColumnNumber    the number of the input column that contains the
     *                           end position. {@code -1} if missing.
     * @param errorExceptionFactory instantiates the exception to thrown in case
     *                              of a formatting error. cannot be {@code null}.
     * @return not a {@code null} if there is enough information to find out the
     *         sample intervals, {@code null} if it is insufficient.
     */
    private static Function<String[], SimpleInterval> composeIntervalBuilder(final int contigColumnNumber,
                                                                             final int startColumnNumber,
                                                                             final int endColumnNumber,
                                                                             final Function<String,RuntimeException> errorExceptionFactory) {
        if (contigColumnNumber == -1 || startColumnNumber == -1 || endColumnNumber == -1) {
            return null;
        }

        return (v) -> {
            final String contig = Objects.requireNonNull(v[contigColumnNumber]);
            final String startString = Objects.requireNonNull(v[startColumnNumber]);
            final String endString = Objects.requireNonNull(v[endColumnNumber]);
            final int start, end;
            try {
                start = Integer.parseInt(startString);
            } catch (final NumberFormatException ex) {
                throw errorExceptionFactory.apply(String.format("start position '%s' cannot be parsed into a integer", startString));
            }
            try {
                end = Integer.parseInt(endString);
            } catch (final NumberFormatException ex) {
                throw errorExceptionFactory.apply(String.format("end position '%s' cannot be parsed into an integer", endString));
            }
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
     *     Its capacity auto-extends as more targets are added to it.
     * </p>
     */
    private static class Buffer {

        /**
         * Set of targets indexed by their names.
         * <p>
         *     The correspondence between targets and columns in {@code counts} is through this map iteration order.
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
         *
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
         * @return never {@code null}.
         */
        private SetUniqueList<Target> getTargets() {
            return targets;
        }

        /**
         *  Returns an array representation of the counts in the buffer.
         *
         *  <p>Each element of array corresponds the the ith target in the buffer.</p>
         *  <p>The result array can be modified at will without altering the buffer, yet the element sub-arrays are
         *  live objects and modifications will change the counts in the buffer</p>
         */
        private double[][] getCounts() {
            return counts.toArray(new double[counts.size()][]);
        }
    }
}
