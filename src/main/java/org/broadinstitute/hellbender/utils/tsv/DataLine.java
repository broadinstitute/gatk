package org.broadinstitute.hellbender.utils.tsv;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.function.Function;

/**
 * Table data-line string array wrapper.
 * <p>
 * This wrapper includes convenience methods to {@link #get get} and {@link #set} values on <i>data-line</i> string array from within
 * record type conversion methods:
 * {@link TableReader#record(DataLine) TableReader.record} and {@link TableWriter#dataLine TableWriter.dataLine}.
 * </p>
 * <p>
 * Apart for {@link #set set} operations, it includes convenience methods to set column values in order
 * of appearance when this is known using {@link #append append}.
 * </p>
 * <p>
 * There are various method overload for type conversion to and from primitive type {@code int} and {@code double}.
 * </p>
 */
public class DataLine {

    /**
     * Holds the values for the data line in construction.
     */
    private final String[] values;

    /**
     * Next appending index used by {@link #append append} methods.
     */
    private int nextIndex = 0;

    /**
     * Reference to the enclosing table's columns.
     */
    private final TableColumns columns;

    /**
     * Reference to the format error exception factory.
     */
    private Function<String, RuntimeException> formatErrorFactory;

    /**
     * Creates a new data-line instance.
     * <p>
     * The value array passed is not copied and will be used directly to store the data-line values.
     * </p>
     * <p>
     * Therefore later modification of its content after creating this {@link #DataLine} may result in
     * breaking the consistency of this instance.
     * </p>
     *
     * @param values             the value array.
     * @param columns            the columns of the table that will enclose this data-line instance.
     * @param formatErrorFactory to be used when there is a column formatting error based on the requested data-type.
     * @throws IllegalArgumentException if {@code columns} or {@code formatErrorFactory} are {@code null}.
     */
    DataLine(final String[] values, final TableColumns columns, final Function<String, RuntimeException> formatErrorFactory) {
        this.values = Utils.nonNull(values, "the value array cannot be null");
        this.columns = Utils.nonNull(columns, "the columns cannot be null");
        this.formatErrorFactory = Utils.nonNull(formatErrorFactory, "the format error factory cannot be null");
        if (values.length != columns.columnCount()) {
            throw new IllegalArgumentException("mismatching value length and column count");
        }
    }

    /**
     * Creates a new data-line instance.
     *
     * @param columns            the columns of the table that will enclose this data-line instance.
     * @param formatErrorFactory to be used when there is a column formatting error based on the requested data-type.
     * @throws IllegalArgumentException if {@code columns} or {@code formatErrorFactory} are {@code null}.
     */
    DataLine(final TableColumns columns, final Function<String, RuntimeException> formatErrorFactory) {
        this(new String[Utils.nonNull(columns, "the columns cannot be null").columnCount()], columns, formatErrorFactory);
    }

    /**
     * Returns a reference to the data-line values after making sure that they are all defined.
     * <p>
     * The value returned is a direct reference to the array and any modification of its contents
     * may result in invalidating the state of this {@link DataLine} object.
     * </p>
     *
     * @return never {@code null} and with no {@code null} elements.
     */
    String[] unpack() {
        for (int i = 0; i < values.length; i++) {
            if (values[i] == null) {
                throw new IllegalStateException(String.format("some data line value remains undefined: e.g. column '%s' index %d", columns.nameAt(i), i));
            }
        }
        if (values[0].startsWith(TableConstants.COMMENT_PREFIX)) {
            throw new IllegalStateException(String.format("the first column value cannot start with the comment prefix string %s; consider to change the record encoding or re-order columns", TableConstants.COMMENT_PREFIX));
        }
        return values;
    }

    /**
     * Sets the string value in the data-line that correspond to a column by its name.
     *
     * @param name  the column name.
     * @param value the new value.
     * @return reference to this builder.
     * @throws IllegalArgumentException if {@code name} is {@code null} or it does not match an actual column name.
     */
    public DataLine set(final String name, final String value) {
        values[columnIndex(name)] = value;
        return this;
    }

    /**
     * Sets the int value in the data-line that correspond to a column by its name.
     *
     * @param name  the column name.
     * @param value the new value.
     * @return reference to this builder.
     * @throws IllegalArgumentException if {@code name} is {@code null} or it does not match an actual column name.
     */
    public DataLine set(final String name, final int value) {
        return set(name, "" + value);
    }

    /**
     * Sets the double value in the data-line that correspond to a column by its name.
     *
     * @param name  the column name.
     * @param value the new value.
     * @return reference to this builder.
     * @throws IllegalArgumentException if {@code name} is {@code null} or it does not match an actual column name.
     */
    public DataLine set(final String name, final double value) {
        return set(name, "" + value);
    }

    /**
     * Sets the string value of a column given its index.
     *
     * @param index the target column index.
     * @param value the new value for that column.
     * @return reference to this builder.
     * @throws IllegalArgumentException if {@code index} is not a valid column index.
     */
    public DataLine set(final int index, final String value) {
        Utils.validIndex(index, values.length);
        values[index] = value;
        return this;
    }

    /**
     * Sets the int value of a column given its index.
     *
     * @param index the target column index.
     * @param value the new value for that column.
     * @return reference to this builder.
     * @throws IllegalArgumentException if {@code index} is not a valid column index.
     */
    public DataLine set(final int index, final int value) {
        return set(index, Integer.toString(value));
    }

    /**
     * Sets the int value of a column given its index.
     *
     * @param index the target column index.
     * @param value the new value for that column.
     * @return reference to this builder.
     * @throws IllegalArgumentException if {@code index} is not a valid column index.
     */
    public DataLine set(final int index, final double value) {
        return set(index, Double.toString(value));
    }

    /**
     * Returns the string value in a column by its index.
     *
     * @param index target column index.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code index} is not a valid column index.
     * @throws IllegalStateException    if the value for that column is undefined ({@code null}).
     */
    public String get(final int index) {
        Utils.validIndex(index, values.length);
        if (values[index] == null) {
            throw new IllegalStateException("requested column value at " + index + " has not been initialized yet");
        }
        return values[index];
    }

    /**
     * Returns the int value in a column by its index.
     *
     * @param index the target column index.
     * @return any int value.
     * @throws IllegalArgumentException if {@code index} is not valid.
     * @throws IllegalStateException    if {@code index} has not been initialized and contains a {@code null}.
     * @throws RuntimeException         if the value at that target column cannot be transform into an integer.
     *                                  The exact class of the exception will depend on the exception factory provided when creating this
     *                                  {@link #DataLine}.
     */
    public int getInt(final int index) {
        try {
            return Integer.parseInt(get(index));
        } catch (final NumberFormatException ex) {
            throw formatErrorFactory.apply(String.format("expected int value for column %s but found %s", columns.nameAt(index), get(index)));
        }
    }

    /**
     * Returns the double value in a column by its index.
     *
     * @param index the target column index.
     * @return any double value.
     * @throws IllegalArgumentException if {@code index} is not valid.
     * @throws IllegalStateException    if {@code index} has not been initialized and contains a {@code null}.
     * @throws RuntimeException         if the value at that target column cannot be transform into a double.
     *                                  The exact class of the exception will depend on the exception factory provided when creating this
     *                                  {@link #DataLine}.
     */
    public double getDouble(final int index) {
        try {
            return Double.parseDouble(get(index));
        } catch (final NumberFormatException ex) {
            throw formatErrorFactory.apply(String.format("expected int value for column %s but found %s", columns.nameAt(index), get(index)));
        }
    }

    /**
     * Returns the string value in a column by its index.
     *
     * @param columnName the target column name.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code columnName} is {@code null} or an unknown column name.
     * @throws IllegalStateException    if that column values is undefined ({@code null}).
     */
    public String get(final String columnName) {
        final int index = columnIndex(columnName);
        if (values[index] == null) {
            throw new IllegalStateException(String.format("the value for column '%s' is undefined", columnName));
        } else {
            return values[index];
        }
    }

    /**
     * Returns the index of a column by its name or fails if invalid or unknown.
     *
     * @param columnName the column name.
     * @return a valid index from 0 to <code>{@link #columns}.size() - 1</code>
     * @throws IllegalArgumentException if {@code columnName} is {@code null} or it not a known column name.
     */
    private int columnIndex(final String columnName) {
        final int index = columns.indexOf(columnName);
        if (index < 0) {
            throw new IllegalArgumentException("there is no such a column: " + columnName);
        }
        return index;
    }

    /**
     * Returns the int value in a column by its index.
     *
     * @param columnName the target column name.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code columnName} is {@code null} or an unknown column name.
     * @throws IllegalStateException    if that column values is undefined ({@code null}).
     * @throws RuntimeException         if the value at that target column cannot be transform into an integer.
     *                                  The exact class of the exception will depend on the exception factory provided when creating this
     *                                  {@link #DataLine}.
     */
    public int getInt(final String columnName) {
        return getInt(columnIndex(columnName));
    }

    /**
     * Returns the double value in a column by its index.
     *
     * @param columnName the target column name.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code columnName} is {@code null} or an unknown column name.
     * @throws IllegalStateException    if that column values is undefined ({@code null}).
     * @throws RuntimeException         if the value at that target column cannot be transform into a double.
     *                                  The exact class of the exception will depend on the exception factory provided when creating this
     *                                  {@link #DataLine}.
     */
    public double getDouble(final String columnName) {
        return getDouble(columnIndex(columnName));
    }

    /**
     * Sets the next string value in the data-line that correspond to a column.
     * <p>
     * The next column index advances so that the following append operations will change the value of
     * the following columns and so forth.
     * </p>
     *
     * @param value the new value.
     * @return reference to this builder.
     * @throws IllegalStateException if the next column to set is beyond the last column.
     */
    public DataLine append(final String value) {
        if (nextIndex == values.length) {
            throw new IllegalStateException("gone beyond of the end of the data-line");
        }
        values[nextIndex++] = value;
        return this;
    }

    /**
     * Sets the next int value in the data-line that correspond to a column.
     * <p>
     * The next column index advances so that the following {@link #append append} will change the value of
     * the following column and so forth.
     * </p>
     *
     * @param value the new value.
     * @return reference to this builder.
     * @throws IllegalStateException if the next column to set is beyond the last column.
     */
    public DataLine append(final int value) {
        return append(Integer.toString(value));
    }

    /**
     * Sets the next double value in the data-line that correspond to a column.
     * <p>
     * The next column index advances so that the following {@link #append append} will change the value of
     * the following column and so forth.
     * </p>
     *
     * @param value the new value.
     * @return reference to this builder.
     * @throws IllegalStateException if the next column to set is beyond the last column.
     */
    public DataLine append(final double value) {
        return append(Double.toString(value));
    }

    /**
     * Sets the next int values in the data-line that correspond to next few columns.
     * <p>
     * The next column index advances so that the following {@link #append append} will change the value of
     * the following column and so forth.
     * </p>
     *
     * @param values the new values.
     * @return reference to this builder.
     * @throws IllegalStateException if this operation goes beyond the last column index.
     */
    public DataLine append(final int... values) {
        for (final int i : Utils.nonNull(values, "the values cannot be null")) {
            append(i);
        }
        return this;
    }

    /**
     * Sets the next string values in the data-line that correspond to next few columns.
     * <p>
     * The next column index advances so that the following {@link #append append} will change the value of
     * the following column and so forth.
     * </p>
     *
     * @param values the new values.
     * @return reference to this builder.
     * @throws IllegalStateException if this operation goes beyond the last column index.
     */
    public DataLine append(final String... values) {
        for (final String v : Utils.nonNull(values, "the values cannot be null")) {
            append(v);
        }
        return this;
    }

    /**
     * Sets the next double values in the data-line that correspond to next few columns.
     * <p>
     * The next column index advances so that the following {@link #append append} will change the value of
     * the following column and so forth.
     * </p>
     * @param values the new values.
     * @return reference to this builder.
     * @throws IllegalStateException if this operation goes beyond the last column index.
     */
    public DataLine append(final double... values) {
        for (final double d : Utils.nonNull(values, "the values cannot be null")) {
            append(d);
        }
        return this;
    }

    /**
     * Changes the index of the next value to set using {@link #append append} operations.
     *
     * @param index the new index.
     * @return a reference to this builder.
     * @throws IllegalArgumentException if {@code index} is not a valid column index.
     */
    public DataLine seek(final int index) {
        // +1 as is valid to seek to the position after the last value.
        nextIndex = Utils.validIndex(index, values.length + 1);
        return this;
    }

    /**
     * Changes the index of the next value to set using {@link #append append} operations.
     *
     * @param columnName the name of the column to seek to.
     * @return a reference to this builder.
     * @throws IllegalArgumentException if {@code columnName} is {@code null} or an unknown column name.
     */
    public DataLine seek(final String columnName) {
        nextIndex = columnIndex(columnName);
        return this;
    }

    /**
     * Returns the current values in a string array.
     * <p>
     * The returned array is a copy that can be modified without changing the state of the data line.
     * </p>
     * @return never {@code null}, but it can contain {@code null}s.
     */
    public String[] toArray() {
        return values.clone();
    }
}
