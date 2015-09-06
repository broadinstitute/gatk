package org.broadinstitute.hellbender.utils.tsv;

import com.opencsv.CSVReader;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;
import java.util.*;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

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
 *             // you can check them here:
 *             &#64Override
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

    /**
     * Holds a reference to the csvReader object use to read and parse the input
     * into {@link String} arrays.
     */
    private CSVReader csvReader;

    /**
     * Holds a reference to the next record.
     * <p>
     * This is {@code null} when we reached the end of the source.
     * </p>
     */
    private R nextRecord;

    /**
     * Creates a new table reader given the input file name.
     * <p>
     * This operation will read the first lines of the input file until the
     * column name header line is found.
     * </p>
     * <p>
     * The source's name used in error reporting is the file's path as returned by
     * {@link File#getPath}.
     * </p>
     *
     * @param file the input file name.
     * @throws IllegalArgumentException if {@code name} is {@code null}.
     * @throws IOException              if any is raised when accessing the file.
     */
    public TableReader(final File file) throws IOException {
        this(Utils.nonNull(file, "the input file cannot be null").getPath(), new FileReader(file));
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
        this(null, sourceReader);
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
    protected TableReader(final String sourceName, final Reader sourceReader) throws IOException {
        Utils.nonNull(sourceReader, "the reader cannot be null");
        this.source = sourceName;
        this.reader = sourceReader instanceof LineNumberReader ? (LineNumberReader) sourceReader : new LineNumberReader(sourceReader);
        this.csvReader = new CSVReader(this.reader, TableUtils.COLUMN_SEPARATOR, TableUtils.QUOTE_CHARACTER, TableUtils.ESCAPE_CHARACTER);
        findAndProcessHeaderLine();
        this.nextRecord = fetchNextRecord();
    }

    /**
     * Process the first lines of the input source until the header line.
     *
     * @throws IOException            if an {@link IOException} occurred when reading from the source.
     * @throws UserException.BadInput if there is formatting error in the input.
     */
    protected void findAndProcessHeaderLine() throws IOException {
        final String[] line = skipCommentLines();
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
    private boolean isCommentLine(final String[] line) {
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
        return new UserException.BadInput(formatExceptionMessage(message));
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
    private String formatExceptionMessage(final String message) {
        final String explanation = message == null ? "" : ": " + message;
        if (source == null) {
            return String.format("format error at line %d" + explanation, reader.getLineNumber());
        } else {
            return String.format("format error in '%s' at line %d" + explanation, source, reader.getLineNumber());
        }
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
     * @return never {@code null}.
     */
    public TableColumnCollection columns() {
        return columns;
    }

    /**
     * Returns the next record form the source.
     *
     * @return {@code null} if there is no more record in the input.
     * @throws IOException if a {@link IOException} was thrown when reading from the input.
     */
    public final R readRecord() throws IOException {
        final R result = nextRecord;
        if (result != null) {
            nextRecord = fetchNextRecord();
        }
        return result;
    }

    /**
     * Fetch the next record from the source.
     *
     * @return {@code null} if there is no more record in the input.
     * @throws IOException if a {@link IOException} was thrown when reading from the input.
     */
    private R fetchNextRecord() throws IOException {

        String[] line;
        while ((line = csvReader.readNext()) != null) {
            if (isCommentLine(line)) {
                continue;
            } else if (line.length != columns.columnCount()) {
                throw formatException(String.format("mismatch between number of values in line (%d) and number of columns (%d)", line.length, columns.columnCount()));
            } else {
                final R result = createRecord(new DataLine(line, columns, this::formatException));
                if (result != null) {
                    return result;
                }
            }
        }
        return null;
    }

    /**
     * Skip comment lines from the output.
     * <p>
     * It returns the contents of the first non comment line found.
     * </p>
     *
     * @return {@code null} if we reached the end of the source, the next non-comment line content otherwise.
     * @throws IOException if it was raised when reading for the source.
     */
    private String[] skipCommentLines() throws IOException {
        String[] line;
        while ((line = csvReader.readNext()) != null) {
            if (!isCommentLine(line)) {
                break;
            }
        }
        return line;
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
                return nextRecord != null;
            }

            @Override
            public R next() {
                if (nextRecord == null) {
                    throw new NoSuchElementException("there is no more record in the input");
                } else {
                    final R result = nextRecord;
                    try {
                        nextRecord = fetchNextRecord();
                    } catch (final IOException ex) {
                        throw new UncheckedIOException(ex);
                    }
                    return result;
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
        return StreamSupport.stream(spliterator(), false);
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
