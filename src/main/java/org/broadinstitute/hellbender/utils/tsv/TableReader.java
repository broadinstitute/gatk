package org.broadinstitute.hellbender.utils.tsv;

import com.opencsv.CSVReader;
import java.nio.file.Path;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Reads the contents of a tab separated value formatted text input into
 * records of an arbitrary type {@link R}.
 * <h3>Format description</h3>
 * <p>
 * Tab separated values may contain any number of <i>comment lines</i> (started with {@link TableUtils#COMMENT_PREFIX}),
 * a column name containing line (aka. the <i>header line</i>) and any number of <i>data lines</i> one per record.
 * </p>
 * <p>While comment lines can contain any sequence of characters, the header and data lines are divided in
 * columns using exactly one {@link TableUtils#COLUMN_SEPARATOR} character.</p>
 * <p>Blank lines are treated is having a single column with the empty string as the only value (or column name)</p>
 * <p>
 * The header line is the first non-comment line, whereas any other non-comment line after that is
 * considered a data line. Comment lines can appear anywhere in the file and their
 * present is ignored by the reader.
 * </p>
 * <p>
 * The header line values, the column names, must all be different (otherwise a formatting exception will be thrown), and
 * all data lines have to have as many values as the header line.
 * </p>
 * <p>Values can be quoted using {@link TableUtils#QUOTE_CHARACTER} becoming necessary when the value contain
 * any special formatting characters like a new-line, the quote character itself, the column separator character or
 * the escape character {@link TableUtils#ESCAPE_CHARACTER}.</p>
 * <p>Within quotes, especial characters must be escaped using the {@link TableUtils#ESCAPE_CHARACTER}</p>
 * <h3>Implementing your own reader</h3>
 * <p>
 * Implementations control how instances of {@link R} are instantiated by extending
 * {@link #createRecord(DataLine) createRecord}. This method is passed an non-null nor null containing
 * {@link DataLine} with exactly as many elements are columns where passed to the
 * {@link #processColumns(TableColumnCollection)} earlier on in the execution.
 * </p>
 * <p>
 * The i-th element in the input array represent the value for the i-th column for that record.
 * </p>
 * <p>
 * The exact list (array) of column names are always accessible through {@link #columns}.
 * </p>
 * <p>
 * Implementations can also override {@link #processColumns} (that by default does nothing) in order to
 * get prepared to received data lines following the that format or simply to verify that
 * that the sequence of column names is expected, throwing an exception if not.
 * </p>
 * <p>
 * For the later, extending classes must use {@link #formatException(String)} to create that exception including an explanation
 * message describing the format violation.
 * </p>
 * <p>
 * Example:
 * <pre>
 *         public class Person {
 *             public final String name;
 *             public final int age;
 *             public final double netWorth;
 *         }
 *
 *         public class PeopleTableReader extends TableReader&lt;Person&gt; {
 *             // ...
 *
 *             // If you don't trust the columns that you are given,
 *             // you can check them here (see also {@link TableUtils#checkMandatoryColumns}):
 *             &#64;Override
 *             public void processColumns(final TableColumns columns) {
 *                 if (!columns.containsExactly("name","age","net.worth"))
 *                     throw formatException("invalid column names")
 *             }
 *
 *             &#64;Override
 *             protected Person createRecord(final DataLine dataLine) {
 *                  return new Person(
 *                      dataLine.get("name"),
 *                      dataLine.getInt("age"),
 *                      dataLine.getDouble("net.worth")
 *                  );
 *             }
 *         }
 *     </pre>
 * </p>
 *
 * @param <R> the record type for the reader.
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public abstract class TableReader<R> implements Closeable, Iterable<R> {

    /**
     * Name of the input source.
     * <p>It can be {@code null} indicating that no name was provided at construction</p>.
     */
    private final String source;

    /**
     * Input text reader.
     * <p>
     * Keeps track of the last line number read for error reporting purposes.
     * </p>
     */
    private final LineNumberReader reader;

    /**
     * Holds a reference to the column names.
     */
    private TableColumnCollection columns;

    private Map<String, String> metadata = new HashMap<>();

    /**
     * Holds a reference to the csvReader object use to read and parse the input
     * into {@link String} arrays.
     */
    private CSVReader csvReader;

    /**
     * Indicates whether the reader has tried to fetch the next record.
     * <p>If {@code true} the content of {@link #nextRecord} represent the next record to be returned
     * by {@link #readRecord} ({@code null} if we reached the end of the table), otherwise {@link #nextRecord} reference
     * is invalid and the next record must be fetched using {@link #fetchNextRecord()}.</p>
     */
    private boolean nextRecordFetched = false;

    /**
     * Holds a reference to the next record.
     * <p>
     * This is {@code null} when we reached the end of the source.
     * </p>
     */
    private R nextRecord;

    /**
     * Creates a new table reader given the input file path.
     * <p>
     * This operation will read the first lines of the input file until the
     * column name header line is found.
     * </p>
     * <p>
     * The source's name used in error reporting is the file's path as returned by
     * {@link File#getPath}.
     * </p>
     *
     * @param path the input file path.
     * @throws IllegalArgumentException if {@code name} is {@code null}.
     * @throws IOException              if any is raised when accessing the file.
     */
    public TableReader(final Path path) throws IOException {
        this(
            Utils.nonNull(path, "the input file cannot be null").toString(),
            IOUtils.makeReaderMaybeGzipped(path),
            new TableReaderOptions()
        );
    }

    public TableReader(final Path path, final TableReaderOptions tableReaderOptions) throws IOException {
        this(
            Utils.nonNull(path, "the input file cannot be null").toString(),
            IOUtils.makeReaderMaybeGzipped(path),
            tableReaderOptions
        );
    }

    /**
     * Creates a new table reader given an input {@link Reader}.
     * <p>
     * This operation will read the first lines of the input file until the
     * column name header line is found.
     * </p>
     *
     * @param sourceReader the source text reader.
     * @throws IOException if any is raised when reading from {@code sourceReader}.
     */
    public TableReader(final Reader sourceReader) throws IOException {
        this(null, sourceReader, new TableReaderOptions());
    }

    public TableReader(final Reader sourceReader, final TableReaderOptions tableReaderOptions) throws IOException {
        this(null, sourceReader, tableReaderOptions);
    }

    protected TableReader(final String sourceName, final Reader sourceReader) throws IOException {
        this(sourceName, sourceReader, new TableReaderOptions());
    }

    /**
     * Creates a new table reader given an input {@link Reader}.
     * <p>
     * It assigns an arbitrary
     * </p>
     *
     * @param sourceName   name of the source to use in error messages. It can be {@code null}, indicating that is anonymous.
     * @param sourceReader reader to the text to process.
     * @throws IllegalArgumentException if {@code sourceReader} is {@code null}.
     * @throws IOException              if is raised when reading from the source.
     */
    protected TableReader(final String sourceName, final Reader sourceReader,
                          final TableReaderOptions tableReaderOptions) throws IOException {
        Utils.nonNull(sourceReader, "the reader cannot be null");

        this.source = sourceName;
        this.reader = sourceReader instanceof LineNumberReader ? (LineNumberReader) sourceReader : new LineNumberReader(sourceReader);
        this.csvReader = new CSVReader(this.reader, TableUtils.COLUMN_SEPARATOR, TableUtils.QUOTE_CHARACTER, TableUtils.ESCAPE_CHARACTER);
        findAndProcessHeaderLine(tableReaderOptions);
        this.nextRecordFetched = false;
    }

    /**
     * Process the first lines of the input source until the header line.
     *
     * @throws IOException            if an {@link IOException} occurred when reading from the source.
     * @throws UserException.BadInput if there is formatting error in the input.
     */
    protected void findAndProcessHeaderLine(final TableReaderOptions tableReaderOptions) throws IOException {
        // get header line, and remap any fields that have been requested
        final String[] line = Arrays.stream(skipCommentLines(tableReaderOptions.headerIsLastComment))
            .map(columnName -> tableReaderOptions.columnRenamer.getOrDefault(columnName, columnName))
            .toArray(String[]::new);
        if (line == null) {
            throw formatException("premature end of table: header line not found");
        } else {
            TableColumnCollection.checkNames(line, UserException.BadInput::new);
            columns = new TableColumnCollection(line);
            processColumns(columns);
        }
    }

    /**
     * Checks whether a line is a comment line or not.
     *
     * @param line input line already split into line-values.
     * @return {@code true} if {@code line} seems to be a comment line.
     */
    protected boolean isCommentLine(final String[] line) {
        return line.length > 0 && line[0].startsWith(TableUtils.COMMENT_PREFIX);
    }

    /**
     * Composes the exception to be thrown due to a formatting error.
     * <p>
     * The input {@code message} can be omitted by providing a {@code null} value.
     * </p>
     *
     * @param message custom error message.
     * @return never {@code null}.
     */
    protected final UserException.BadInput formatException(final String message) {
        return new UserException.BadInput(formatExceptionMessageWithLocationInfo(message));
    }

    /**
     * Composes the exception to be thrown due to a formatting error.
     * <p>
     * The input {@code message} can be omitted by providing a {@code null} value.
     * </p>
     *
     * @param message custom error message.
     * @return never {@code null}.
     */
    protected final UserException.BadInput formatExceptionWithoutLocation(final String message) {
        return new UserException.BadInput(formatExceptionMessageWithoutLocationInfo(message));
    }

    /**
     * Composes the error exception message string.
     * <p>
     * The input {@code message} can be omitted by providing a {@code null} value.
     * </p>
     *
     * @param message custom error message.
     * @return never {@code null}.
     */
    private String formatExceptionMessageWithLocationInfo(final String message) {
        final String explanation = message == null ? "" : ": " + message;
        if (source == null) {
            return String.format("format error at line %d" + explanation, reader.getLineNumber());
        } else {
            return String.format("format error in '%s' at line %d" + explanation, source, reader.getLineNumber());
        }
    }

    /**
     * Composes the error exception message string.
     * <p>
     * The input {@code message} can be omitted by providing a {@code null} value.
     * </p>
     *
     * @param message custom error message.
     * @return never {@code null}.
     */
    private String formatExceptionMessageWithoutLocationInfo(final String message) {
        return "format error: " + message;
    }

    /**
     * Process the header line's column names.
     * <p>
     * Implementations must use {@link #formatException(String)} to create the exception to throw in case
     * there is any formatting issue.
     * </p>
     *
     * @param tableColumns columns found in the input. It is guarantee not to be
     *                     a {@code null} and not to contain any {@code null} values.
     * @throws UserException.BadInput if there is a formatting issue.
     */
    protected void processColumns(@SuppressWarnings("unused") final TableColumnCollection tableColumns) {
        // nothing by default.
    }

    /**
     * Returns the column collection for this reader.
     * @throws IllegalStateException if this methods is invoked before the table column can be determined
     *  (e.g. when processing a comment line before the header extending {@link #processCommentLine(String, long)}).
     * @return never {@code null}.
     */
    public TableColumnCollection columns() {
        Utils.validate(columns != null, "columns are null");
        return columns;
    }

    /**
     * Returns the next record form the source.
     *
     * @return {@code null} if there is no more record in the input.
     * @throws IOException if a {@link IOException} was thrown when reading from the input.
     */
    public final R readRecord() throws IOException {
        if (nextRecordFetched == false) {
            nextRecord = fetchNextRecord();
        }
        if (nextRecord != null) {
            nextRecordFetched = false;
            return nextRecord;
        } else {
            return null;
        }
    }

    /**
     * Reads the record from a string rather than from the input reader.
     *
     * @return {@code null} for comment or header lines, a non-null record otherwise.
     */
    public final R readRecord(final String line) {
        try {
            final String[] fields = csvReader.getParser().parseLine(line);
            if (isCommentLine(fields) || isHeaderLine(fields)) {
                return null;
            } else if (fields.length != columns.columnCount()) {
                throw formatExceptionWithoutLocation("invalid number of columns");
            } else {
                return createRecord(new DataLine(fields, columns, this::formatExceptionWithoutLocation));
            }
        } catch (final IOException ex) {
            throw new GATKException("the single line input is in fact a multi-line entry");
        }
    }

    /**
     * Fetch the next record from the source.
     *
     * @return {@code null} if there is no more record in the input.
     * @throws IOException if a {@link IOException} was thrown when reading from the input.
     */
    private R fetchNextRecord() throws IOException {
        nextRecordFetched = true;
        String[] line;
        while ((line = csvReader.readNext()) != null) {
            if (isCommentLine(line)) {
                processCommentLine(line, reader.getLineNumber());
            } else if (!isHeaderLine(line)) {
                if (line.length != columns.columnCount()) {
                    throw formatException(String.format("mismatch between number of values in line (%d) and number of columns (%d)", line.length, columns.columnCount()));
                } else {
                    final R result = createRecord(new DataLine(reader.getLineNumber(), line, columns, this::formatException));
                    if (result != null) {
                        return result;
                    }
                }
            }
        }
        return null;
    }

    private void processCommentLine(final String[] line, final long lineNumber) {
        final StringBuilder builder = new StringBuilder();
        builder.append(line[0].substring(TableUtils.COMMENT_PREFIX.length()));
        for (int i = 1; i < line.length; i++)
            builder.append(TableUtils.COLUMN_SEPARATOR_STRING).append(line[i]);
        processCommentLine(builder.toString(), lineNumber);
    }

    /**
     * Called with the content of the comment line every time one is found in the input.
     * <p>
     *     The comment prefix string ({@link TableUtils#COMMENT_PREFIX}) is not included in the
     *     input string.
     * </p>
     * <p>
     *     Notice that since comments might be present before the header line, so the result of invoking {@link #columns()} is undefined.
     * </p>
     * @param commentText the input comment string.
     * @param lineNumber the source line number that contains the comment line.
     */
    protected void processCommentLine(final String commentText, final long lineNumber) {
        if (commentText.startsWith(TableWriter.METADATA_TAG)) {
            final String[] keyAndValue = commentText.substring(TableWriter.METADATA_TAG.length()).split("=");
            metadata.put(keyAndValue[0], keyAndValue[1]);
        }
        // do nothing with non-metadata lines by default
    }

    /**
     * Determines whether a line is a repetition of the header.
     *
     * <p>
     *     By default, input lines that match the header exactly are ignored. By overriding this method
     *     extending classes may change what is interpretated as a repetition of the header (e.g. just treat such
     *     lines as regular data line)
     * </p>
     * @param line the input line.
     * @return {@code true} if the input line is a header line and it should be ignored.
     */
    protected boolean isHeaderLine(final String[] line) {
        return columns.matchesExactly(line);
    }

    /**
     * Skip comment lines from the output.
     * <p>
     * It returns the contents of the first non comment line found.  As a side effect, it builds a metadata map
     * of key-value pairs from comment lines with the <METADATA> tag from TableWriter
     * </p>
     *
     * @return {@code null} if we reached the end of the source, the next non-comment line content otherwise.
     * @throws IOException if it was raised when reading for the source.
     */
    private String[] skipCommentLines(final boolean headerIsLastComment) throws IOException {
        String[] line;
        if(headerIsLastComment) {
            String[] maybeHeader = null;
            while((line = csvReader.readNext()) != null) {
                if(maybeHeader != null) {
                    // the previous comment was really a comment, process it
                    processCommentLine(maybeHeader, reader.getLineNumber() - 1);
                }
                if(isCommentLine(line)) {
                    maybeHeader = line;
                } else {
                    break;
                }
            }
            // remove the first comment and return header
            if(maybeHeader == null || maybeHeader.length == 0) {
                throw formatException("premature end of table: header line empty or not found");
            }
            maybeHeader[0] = maybeHeader[0].substring(1 + maybeHeader[0].indexOf(TableUtils.COMMENT_PREFIX)).trim();
            if(maybeHeader[0].length() == 0) {
                // remove empty first "word" which was only comment
                maybeHeader = Arrays.stream(maybeHeader).skip(1).toArray(String[]::new);
                if(maybeHeader.length == 0) {
                    throw formatException("header line was empty");
                }
            }
            return maybeHeader;
        } else {
            while ((line = csvReader.readNext()) != null) {
                if (isCommentLine(line)) {
                    processCommentLine(line, reader.getLineNumber());
                } else {
                    return line;
                }
            }
            throw formatException("premature end of table: header line not found");
        }
    }

    /**
     * Transforms a data-line column values into a record.
     * <p>
     * Implementation should use {@link #formatException(String)} to indicate a formatting error
     * that makes impossible to create a instance of {@link R} given the input line values {@code dataLine}.
     * </p>
     *
     * @param dataLine values corresponding to the column names that was passed earlier to {@link #processColumns(TableColumnCollection)}}.
     *                 it is guaranteed to not be {@code null}, contain no {@code null} value and have the same columns
     *                 as this readers' (accessible through {@link #columns}).
     * @return never {@code null}.
     */
    protected abstract R createRecord(final DataLine dataLine);

    @Override
    public void close() throws IOException {
        csvReader.close();
    }

    /**
     * Returns an iterator on the remaining records in
     * the input.
     *
     * @return never {@code null}.
     */
    @Override
    public Iterator<R> iterator() {
        return new Iterator<R>() {

            @Override
            public boolean hasNext() {
                if (!nextRecordFetched) {
                    try {
                        nextRecord = fetchNextRecord();
                    } catch (final IOException ex) {
                        throw new UncheckedIOException(ex);
                    }
                }
                return nextRecord != null;
            }

            @Override
            public R next() {
                if (!nextRecordFetched) {
                    try {
                        nextRecord = fetchNextRecord();
                    } catch (final IOException ex) {
                        throw new UncheckedIOException(ex);
                    }
                }
                if (nextRecord == null) {
                    throw new NoSuchElementException("there is no more record in the input");
                } else {
                    nextRecordFetched = false;
                    return nextRecord;
                }
            }
        };
    }

    @Override
    public Spliterator<R> spliterator() {
        return Spliterators.spliteratorUnknownSize(iterator(), Spliterator.ORDERED | Spliterator.NONNULL);
    }

    /**
     * Returns an stream on the reaming records in the source.
     * <p>
     * Notice that the returned stream will consume records as if you were calling {@link #readRecord} directly.
     * </p>
     * <p>
     * Any {@link IOException} raised when using the stream will be propagated up wrapped in a {@link UncheckedIOException}.
     * </p>
     * <p>
     * Any format exception will still be indicated with a {@link UserException.BadInput}.
     * </p>
     *
     * @return never {@code null}.
     */
    public Stream<R> stream() {
        return Utils.stream(this);
    }


    /**
     * Read the remaining records into a list.
     * <p>
     *     Notice that this operation does not close the reader.
     * </p>
     *
     * @return never {@code null}, but potentially empty.
     */
    public List<R> toList() {
        return stream().collect(Collectors.toList());
    }

    public Map<String, String> getMetadata() {
        return metadata;
    }

    /**
     * Returns the reader source name.
     *
     * @return null if the source name cannot be determined.
     */
    public String getSource() {
        return source;
    }
}
