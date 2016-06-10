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
 * record data-line in the output by overriding {@link #composeLine(R,DataLine)}.
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
 *             public MyRecordWriter(final File file) {
 *                 super(file, new TableColumnCollection("name","age","net.worth"));
 *             }
 *
 *             &#64;Override
 *             protected void dataLine(final Person person, final DataLine dataLine) {
 *                  dataLine.setAll(person.name, "" + person.age, "" + person.netWorth);
 *             }
 *         }
 *     </pre>
 * </p>
 * <p>
 * You must use the {@link DataLine} instance passed and no other.
 * </p>
 * <p>
 * Instead of passing all the values as converted string in column order you may opt to use {@link DataLine#set}
 * method family to set values one by one using the column index or column name like so:
 * </p>
 * <p>
 * Example (using the column index):
 * <pre>
 *          &#64;Override
 *          protected void composeLine(final Person person, final DataLine dataLine) {
 *              dataLine
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
 *          protected void composeLine(final Person person, final DataLine dataLine) {
 *              dataLine
 *                  .set("name",person.name)
 *                  .set("age",person.age)
 *                  .set("net.worth",person.netWorth);
 *          }
 * </pre>
 * Notice that you don't need to explicitly convert neither the age nor the net-worth into a
 * string thanks to {@link DataLine#set set} various overloads.
 * </p>
 * <p>
 * Alternatively, if you know the column order, that should quite often the case, you can avoid
 * indexing all together using {@link DataLine#append append} operations instead:
 * <pre>
 *         &#64;Override
 *          protected void composeLine(final Person person, final DataLine dataLine) {
 *              dataLine
 *                  .append(person.name)
 *                  .append(person.age)
 *                  .append(person.netWorth);
 *          }
 *     </pre>
 * </p>
 * <p>
 * At any time the implementation can query the correspondence between column names and position within the data-line
 * by querying the {@link TableColumnCollection} object directly that can be obtained from the dataLine's {@link #columns} field.
 * </p>
 * <p>
 * Example (using column names):
 * <pre>
 *          &#64;Override
 *          protected void composeLine(final Person person, final DataLine dataLine) {
 *              dataLine
 *                .set("name",person.name)
 *                .set("age",person.age);
 *
 *              if (dataLine.columns().contains("net.worth"))
 *                dataLine.set("net.worth",person.netWorth);
 *          }
 * </pre>
 * </p>
 *
 * @param <R> the row record type.
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public abstract class TableWriter<R> implements Closeable {

    private long lineNumber;

    /**
     * Csv writer use to do the actual writing.
     */
    private final CSVWriter writer;

    /**
     * The table column names.
     */
    private final TableColumnCollection columns;

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
    public TableWriter(final File file, final TableColumnCollection tableColumns) throws IOException {
        this(new FileWriter(Utils.nonNull(file, "The file cannot be null.")), tableColumns);
    }

    /**
     * Creates a new table writer given the destination writer and column names.
     *
     * @param writer  the destination writer.
     * @param columns the table column names.
     * @throws IllegalArgumentException if either {@code writer} or {@code columns} are {@code null}.
     * @throws IOException              if one was raised when opening the the destination file for writing.
     */
    public TableWriter(final Writer writer, final TableColumnCollection columns) throws IOException {

        this.columns = Utils.nonNull(columns, "The columns cannot be null.");
        this.writer = new CSVWriter(Utils.nonNull(writer, "the input writer cannot be null"),
                TableUtils.COLUMN_SEPARATOR, TableUtils.QUOTE_CHARACTER, TableUtils.ESCAPE_CHARACTER);
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
        Utils.nonNull(comment, "The comment cannot be null.");
        writer.writeNext(new String[]{TableUtils.COMMENT_PREFIX + comment}, false);
        lineNumber++;
    }

    /**
     * Writes a new record.
     *
     * @param record the record to write.
     * @throws IOException              if it was raised when writing the record.
     * @throws ClassCastException       if {@code record} is of the correct type
     *                                  for this writer.
     * @throws IllegalArgumentException if {@code record} is {@code null} or it is not a valid record
     *                                  as per the implementation of this writer (see {@link #composeLine}).
     */
    public void writeRecord(final R record) throws IOException {
        Utils.nonNull(record, "The record cannot be null.");
        writeHeaderIfApplies();
        final DataLine dataLine = new DataLine(lineNumber + 1, columns,IllegalArgumentException::new);
        composeLine(record,dataLine);
        writer.writeNext(dataLine.unpack(), false);
        lineNumber++;
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
        Utils.nonNull(records, "The record iterable cannot be null.");
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
            lineNumber++;
        }
        headerWritten = true;
    }

    /**
     * Composes the data-line to write into the output to represent a given record
     * <p>
     * Also the first element cannot contain the {@link TableUtils#COMMENT_PREFIX comment prefix}.
     * If that is a genuine valid value for the first column you shall consider to re-order the columns or
     * change the encoding of the first column to avoid this issue.
     * </p>
     * <p>
     * Both inputs, {@code record} and {@code dataLine} are guaranteed not to be {@code null}s.
     * </p>
     *
     * @param record the record to write into the data-line.
     * @param dataLine the destination data-line object.
     * @throws ClassCastException       if {@code record} is of the correct type
     *                                  for this writer.
     * @throws IllegalArgumentException if there is some conversion issue that does
     *                                  not allow the current write to generate a valid string array to encode the record.
     */
    protected abstract void composeLine(final R record, final DataLine dataLine);
}
