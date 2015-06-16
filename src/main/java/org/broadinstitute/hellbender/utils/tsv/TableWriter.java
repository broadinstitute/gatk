package org.broadinstitute.hellbender.utils.tsv;

import com.opencsv.CSVWriter;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.*;

/**
 * Class to write tab separated value files.
 * <p>
 * The column (and they names) are passed in the constructor parameter along the output {@link File file}
 * or {@link Writer writer}.
 * </p>
 * <p>
 * Extending classes must indicate how we can transcribe row record or type {@link R} to the corresponding
 * record data-line in the output by overriding {@link #dataLine(R)}.
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
 *         public class PeopleTableWriter extends TableWriter&lt;Person&gt; {
 *
 *             public MyRecordWriter() {
 *                 super(new TableColumns("name","age","net.worth"));
 *             }
 *
 *             &#64;Override
 *             protected DataLine dataLine(final Person person) {
 *                  return dataLine(person.name, "" + person.age, "" + person.netWorth);
 *             }
 *         }
 *     </pre>
 * </p>
 * <p>
 * You must instantiate the result {@link DataLine} using the {@link #dataLine()} method.
 * </p>
 * <p>
 * Instead of passing all the values as converted string in column order you may opt to use {@link DataLine#set}
 * method family to set values one by one using the column index or column name like so:
 * </p>
 * <p>
 * Example (using the column index):
 * <pre>
 *          &#64;Override
 *          protected DataLine dataLine(final Person person) {
 *              return dataLine()
 *                  .set(0,person.name)
 *                  .set(1,person.age)
 *                  .set(2,person.netWorth);
 *          }
 *      </pre>
 * </p>
 * <p>
 * Example (using column names):
 * <pre>
 *          &#64;Override
 *          protected DataLine dataLine(final Person person) {
 *              return dataLine()
 *                  .set("name",person.name)
 *                  .set("age",person.age)
 *                  .set("net.worth",person.netWorth);
 *          }
 * </pre>
 * Notice that you don't need to explicitly convert neither the age nor the net-worth into a
 * string thanks to {@link DataLine#set set} overloads.
 * </p>
 * <p>
 * Alternatively, if you know the column order, that should quite often the case, you can avoid indexing all together
 * using {@link DataLine#append append} operations instead:
 * <pre>
 *         &#64;Override
 *          protected DataLine dataLine(final Person person) {
 *              return dataLine()
 *                  .append(person.name)
 *                  .append(person.age)
 *                  .append(person.netWorth);
 *          }
 *     </pre>
 * </p>
 * <p>
 * At any time the implementation can query the correspondence between column names and position within the data-line
 * by querying the {@link TableColumns} object directly referenced by the {@link #columns} field.
 * </p>
 *
 * @param <R> the row record type.
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public abstract class TableWriter<R> implements Closeable {

    /**
     * Csv writer use to do the actual writing.
     */
    private final CSVWriter writer;

    /**
     * The table column names.
     */
    private final TableColumns columns;

    /**
     * Whether the header column name line has been written or not.
     */
    private boolean headerWritten = false;

    /**
     * Creates a new table writer given the file and column names.
     *
     * @param file         the destination file.
     * @param tableColumns the table column names.
     * @throws IllegalArgumentException if either {@code file} or {@code tableColumns} are {@code null}.
     * @throws IOException              if one was raised when opening the the destination file for writing.
     */
    public TableWriter(final File file, final TableColumns tableColumns) throws IOException {
        this(new FileWriter(Utils.nonNull(file, "the file cannot be null")), tableColumns);
    }

    /**
     * Creates a new table writer given the destination writer and column names.
     *
     * @param writer  the destination writer.
     * @param columns the table column names.
     * @throws IllegalArgumentException if either {@code writer} or {@code columns} are {@code null}.
     * @throws IOException              if one was raised when opening the the destination file for writing.
     */
    public TableWriter(final Writer writer, final TableColumns columns) throws IOException {

        this.columns = Utils.nonNull(columns, "the columns cannot be null");
        this.writer = new CSVWriter(Utils.nonNull(writer, "the input writer cannot be null"),
                TableConstants.COLUMN_SEPARATOR, TableConstants.QUOTE_CHARACTER, TableConstants.ESCAPE_CHARACTER);
    }

    /**
     * Writes a comment into the output.
     * <p>
     * This can be invoked at any time; comment lines can be present anywhere in the file.
     * </p>
     * <p>
     * Comments written before any record, will be output
     * </p>
     *
     * @param comment the comment to write out.
     * @throws IllegalArgumentException if {@code comment} is {@code null}.
     * @throws IOException              if any was raised by this operation.
     */
    public final void writeComment(final String comment) throws IOException {
        Utils.nonNull(comment, "the comment cannot be null");
        writer.writeNext(new String[]{TableConstants.COMMENT_PREFIX + comment}, false);
    }

    /**
     * Writes a new record.
     *
     * @param record the record to write.
     * @throws IOException              if it was raised when writing the record.
     * @throws ClassCastException       if {@code record} is of the correct type
     *                                  for this writer.
     * @throws IllegalArgumentException if {@code record} is {@code null} or it is not a valid record
     *                                  as per the implementation of this writer (see {@link #dataLine}).
     */
    public final void writeRecord(final R record) throws IOException {
        Utils.nonNull(record, "the record cannot be null");
        writeHeaderIfApplies();
        writer.writeNext(dataLine(record).unpack(), false);
    }

    /**
     * Write all the records in a {@link Iterable}.
     * <p>
     * Records are written in the order they appear in the input {@link Iterable}.
     * </p>
     *
     * @param records to write.
     * @throws IOException              if any raised when writing any of the records.
     * @throws ClassCastException       if {@code record} is of the correct type
     *                                  for this writer.
     * @throws IllegalArgumentException if {@code records} is {@code null} or it contains
     *                                  some values that would cause such an exception when {@link #writeRecord} is call on
     *                                  that value. Previous record in the iterable would have been already written by then.
     */
    public final void writeAllRecords(final Iterable<R> records) throws IOException {
        Utils.nonNull(records, "the record iterable cannot be null");
        for (final R record : records) {
            writeRecord(record);
        }
    }

    @Override
    public final void close() throws IOException {
        writeHeaderIfApplies();
        writer.close();
    }

    /**
     * Writes the header if it has not been written already.
     * <p>
     * The header is written automatically before the first record is written or when the writer is closed
     * and no record was written.
     * </p>
     * <p>
     * Comments written using {@link #writeComment} before any record will precede the header
     * unless you invoke your method first.
     * </p>
     * <p>
     * Once the header line has been written, invoking this method does not have any effect.
     * </p>
     *
     * @throws IOException if any raised when writing into the destination writer.
     */
    public void writeHeaderIfApplies() throws IOException {
        if (!headerWritten) {
            writer.writeNext(columns.names().toArray(new String[columns.columnCount()]), false);
        }
        headerWritten = true;
    }

    /**
     * Returns the array of values that represent a record in the output text format.
     * <p>
     * This method must not return {@code null} {@link DataLine data-line}, and all its
     * values must have been set.
     * </p>
     * <p>
     * Also the first element cannot contain the {@link TableConstants#COMMENT_PREFIX comment prefix}.
     * If that is a genuine valid value for the first column you shall consider to re-order the columns or
     * change the encoding of the first column to avoid this issue.
     * </p>
     *
     * @return never {@code null}.
     * @throws ClassCastException       if {@code record} is of the correct type
     *                                  for this writer.
     * @throws IllegalArgumentException if there is some conversion issue that does
     *                                  not allow the current write to generate a valid string array to encode the record.
     */
    protected abstract DataLine dataLine(final R record);

    /**
     * Instantiates a record {@link DataLine} appropriate for this writer column set.
     *
     * @return never {@code null}.
     */
    protected DataLine dataLine() {
        return new DataLine(columns, IllegalArgumentException::new);
    }
}
