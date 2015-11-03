package org.broadinstitute.hellbender.utils.tsv;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.function.Function;

/**
 * Table data-line string array wrapper.
 * <p>
 * This wrapper includes convenience methods to {@link #get get} and {@link #set} values on <i>data-line</i> string array from within
 * record type conversion methods:
 * {@link TableReader#createRecord(DataLine) TableReader.record} and {@link TableWriter#composeLine TableWriter.composeLine}.
 * </p>
 * <p>
 * Apart for {@link #set set} operations, it includes convenience methods to set column values in order
 * of appearance when this is known using {@link #append append}.
 * </p>
 * <p>
 * There are various method overloads for type conversion to and from primitive type {@code int} and {@code double}.
 * </p>
 * <p>
 * You can use {@link #columns()} to obtain the corresponding {@link TableColumnCollection} and query the presence of
 * and the index of columns.
 * </p>
 */
public final class DataLine {

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
    private final TableColumnCollection columns;

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
     * Therefore later modification of its content after creating this data-line may result in
     * breaking the consistency of this instance.
     * </p>
     *
     * @param values             the value array.
     * @param columns            the columns of the table that will enclose this data-line instance.
     * @param formatErrorFactory to be used when there is a column formatting error based on the requested data-type.
     * @throws IllegalArgumentException if {@code columns} or {@code formatErrorFactory} are {@code null}.
     */
    DataLine(final String[] values, final TableColumnCollection columns, final Function<String, RuntimeException> formatErrorFactory) {
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
    public DataLine(final TableColumnCollection columns, final Function<String, RuntimeException> formatErrorFactory) {
        this(new String[Utils.nonNull(columns, "the columns cannot be null").columnCount()], columns, formatErrorFactory);
    }

    /**
     * Returns the column collection for this data-line instance.
     *
     * @return never {@code null}.
     */
    public TableColumnCollection columns() {
        return columns;
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
        return values;
    }

    /**
     * Sets the string value in the data-line that correspond to a column by its name.
     *
     * @param name  the column name.
     * @param value the new value.
     * @return reference to this data-line.
     * @throws IllegalArgumentException if {@code name} is {@code null} or it does not match an actual column name.
     */
    public DataLine set(final String name, final String value) {
        return set(columnIndex(name), value);
    }

    /**
     * Sets the boolean value in the data-line that correspond to a column by its name.
     *
     * @param name  the column name.
     * @param value the new value.
     * @return reference to this data-line.
     * @throws IllegalArgumentException if {@code name} is {@code null} or it does not match an actual column name.
     */
    public DataLine set(final String name, final boolean value) {
        return set(columnIndex(name), Boolean.toString(value));
    }

    /**
     * Sets the int value in the data-line that correspond to a column by its name.
     *
     * @param name  the column name.
     * @param value the new value.
     * @return reference to this data-line.
     * @throws IllegalArgumentException if {@code name} is {@code null} or it does not match an actual column name.
     */
    public DataLine set(final String name, final int value) {
        return set(name, Integer.toString(value));
    }

    /**
     * Sets the long value in the data-line that correspond to a column by its name.
     *
     * @param name  the column name.
     * @param value the new value.
     * @return reference to this data-line.
     * @throws IllegalArgumentException if {@code name} is {@code null} or it does not match an actual column name.
     */
    public DataLine set(final String name, final long value) {
        return set(name, Long.toString(value));
    }

    /**
     * Sets the double value in the data-line that correspond to a column by its name.
     *
     * @param name  the column name.
     * @param value the new value.
     * @return reference to this data-line.
     * @throws IllegalArgumentException if {@code name} is {@code null} or it does not match an actual column name.
     */
    public DataLine set(final String name, final double value) {
        return set(name, Double.toString(value));
    }

    /**
     * Sets the string value of a column given its index.
     *
     * @param index the target column index.
     * @param value the new value for that column.
     * @return reference to this data-line.
     * @throws IllegalArgumentException if {@code index} is not a valid column index.
     */
    public DataLine set(final int index, final String value) {
        Utils.validIndex(index, values.length);
        if (index == 0 && value != null) {
            if (value.startsWith(TableUtils.COMMENT_PREFIX)) {
                throw new IllegalArgumentException("the value of the first column cannot start with the comment prefix: " + TableUtils.COMMENT_PREFIX);
            }
        }
        values[index] = value;
        return this;
    }

    /**
     * Sets the int value of a column given its index.
     *
     * @param index the target column index.
     * @param value the new value for that column.
     * @return reference to this data-line.
     * @throws IllegalArgumentException if {@code index} is not a valid column index.
     */
    public DataLine set(final int index, final int value) {
        return set(index, Integer.toString(value));
    }

    /**
     * Sets the long value of a column given its index.
     *
     * @param index the target column index.
     * @param value the new value for that column.
     * @return reference to this data-line.
     * @throws IllegalArgumentException if {@code index} is not a valid column index.
     */
    public DataLine set(final int index, final long value) {
        return set(index, Long.toString(value));
    }

    /**
     * Sets the int value of a column given its index.
     *
     * @param index the target column index.
     * @param value the new value for that column.
     * @return reference to this data-line.
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
     * Returns the long value in a column by its index.
     *
     * @param index the target column index.
     * @return any long value.
     * @throws IllegalArgumentException if {@code index} is not valid.
     * @throws IllegalStateException    if {@code index} has not been initialized and contains a {@code null}.
     * @throws RuntimeException         if the value at that target column cannot be transform into a long.
     *                                  The exact class of the exception will depend on the exception factory provided when creating this
     *                                  {@link #DataLine}.
     */
    public long getLong(final int index) {
        try {
            return Long.parseLong(get(index));
        } catch (final NumberFormatException ex) {
            throw formatErrorFactory.apply(String.format("expected long value for column %s but found %s", columns.nameAt(index), get(index)));
        }
    }

    /**
     * Returns the boolean value in a column by its index.
     *
     * @param index the target column index.
     * @return any boolean value.
     * @throws IllegalArgumentException if {@code index} is not valid.
     * @throws IllegalStateException    if {@code index} has not been initialized and contains a {@code null}.
     * @throws RuntimeException         if the value at that target column cannot be transform into a boolean.
     *                                  The exact class of the exception will depend on the exception factory provided when creating this
     *                                  {@link #DataLine}.
     */
    public boolean getBoolean(final int index) {
        try {
            String b = get(index);
            final String FALSE = Boolean.toString(false);
            final String TRUE = Boolean.toString(true);

            if (!(b.equals(TRUE)) && !(b.equals(FALSE))) {
                throw formatErrorFactory.apply(String.format("Boolean value must be '%s' or '%s' (case sensitive) for column %s but found %s", TRUE, FALSE,
                        columns.nameAt(index), get(index)));
            }
            return Boolean.parseBoolean(get(index));
        } catch (final NumberFormatException ex) {
            throw formatErrorFactory.apply(String.format("expected boolean value for column %s but found %s", columns.nameAt(index), get(index)));
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
     * Returns the string value in a column by its index.  If column name does not exist, returns the default value.
     *
     * @param columnName the target column name.
     * @param defaultValue default value to use if columnName not found.
     * @return {@code null} iff {@code defaultValue == null} and there is not such a column with name {@code columnName}.
     */
    public String get(final String columnName, final String defaultValue) {
        final int index = columns.indexOf(columnName);
        if (index < 0) {
            return defaultValue;
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
     * @return any int value.
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
     * Returns the long value in a column by its index.
     *
     * @param columnName the target column name.
     * @return any long value.
     * @throws IllegalArgumentException if {@code columnName} is {@code null} or an unknown column name.
     * @throws IllegalStateException    if that column values is undefined ({@code null}).
     * @throws RuntimeException         if the value at that target column cannot be transform into a long.
     *                                  The exact class of the exception will depend on the exception factory provided when creating this
     *                                  {@link #DataLine}.
     */
    public long getLong(final String columnName) {
        return getLong(columnIndex(columnName));
    }

    /**
     * Returns the boolean value in a column by its index.
     *
     * @param columnName the target column name.
     * @return any boolean value.
     * @throws IllegalArgumentException if {@code columnName} is {@code null} or an unknown column name.
     * @throws IllegalStateException    if that column values is undefined ({@code null}).
     * @throws RuntimeException         if the value at that target column cannot be transform into a boolean.
     *                                  The exact class of the exception will depend on the exception factory provided when creating this
     *                                  {@link #DataLine}.
     */
    public boolean getBoolean(final String columnName) {
        return getBoolean(columnIndex(columnName));
    }

    /**
     * Returns the double value in a column by its index.
     *
     * @param columnName the target column name.
     * @return any double value.
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
     * Returns the string value of a column by its name.
     * <p>
     * The target column name is resolved as the string returned by {@link Object#toString toString} applied
     * to the input enum.
     * </p>
     *
     * @param column the enum value that provides the name of the column.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code column} is {@code null} or it does not point to a
     *         known column name.
     * @throws IllegalStateException    if that column values is undefined ({@code null}).
     * @throws RuntimeException         if the value at that target column cannot be transform into a double.
     *                                  The exact class of the exception will depend on the exception factory provided when creating this
     *                                  {@link #DataLine}.
     */
    public String get(final Enum<?> column) {
        return get(Utils.nonNull(column).toString());
    }

    /**
     * Returns the int value of a column by its name.
     * <p>
     * The target column name is resolved as the string returned by {@link Object#toString toString} applied
     * to the input enum.
     * </p>
     *
     * @param column the enum value that provides the name of the column.
     * @return any int value.
     * @throws IllegalArgumentException if {@code column} is {@code null} or it does not point to a
     *         known column name.
     * @throws IllegalStateException    if that column values is undefined ({@code null}).
     * @throws RuntimeException         if the value at that target column cannot be transform into a double.
     *                                  The exact class of the exception will depend on the exception factory provided when creating this
     *                                  {@link #DataLine}.
     */
    public int getInt(final Enum<?> column) {
        return getInt(Utils.nonNull(column).toString());
    }

    /**
     * Returns the double value of a column by its name.
     * <p>
     * The target column name is resolved as the string returned by {@link Object#toString toString} applied
     * to the input enum.
     * </p>
     *
     * @param column the enum value that provides the name of the column.
     * @return any double value.
     * @throws IllegalArgumentException if {@code column} is {@code null} or it does not point to a
     *         known column name.
     * @throws IllegalStateException    if that column values is undefined ({@code null}).
     * @throws RuntimeException         if the value at that target column cannot be transform into a double.
     *                                  The exact class of the exception will depend on the exception factory provided when creating this
     *                                  {@link #DataLine}.
     */
    public double getDouble(final Enum<?> column) {
        return getDouble(Utils.nonNull(column).toString());
    }

    /**
     * Sets the next string value in the data-line that correspond to a column.
     * <p>
     * The next column index advances so that the following append operations will change the value of
     * the following columns and so forth.
     * </p>
     *
     * @param value the new value.
     * @return reference to this data-line.
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
     * @return reference to this data-line.
     * @throws IllegalStateException if the next column to set is beyond the last column.
     */
    public DataLine append(final int value) {
        return append(Integer.toString(value));
    }

    /**
     * Sets the next long value in the data-line that correspond to a column.
     * <p>
     * The next column index advances so that the following {@link #append append} will change the value of
     * the following column and so forth.
     * </p>
     *
     * @param value the new value.
     * @return reference to this data-line.
     * @throws IllegalStateException if the next column to set is beyond the last column.
     */
    public DataLine append(final long value) {
        return append(Long.toString(value));
    }

    /**
     * Sets the next long values in the data-line that correspond to next few columns.
     * <p>
     * The next column index advances so that the following {@link #append append} will change the value of
     * the following column and so forth.
     * </p>
     *
     * @param values the new values.
     * @return reference to this data-line.
     * @throws IllegalStateException if this operation goes beyond the last column index.
     */
    public DataLine append(final long... values) {
        for (final long l : Utils.nonNull(values, "the values cannot be null")) {
            append(l);
        }
        return this;
    }

    /**
     * Sets the next double value in the data-line that correspond to a column.
     * <p>
     * The next column index advances so that the following {@link #append append} will change the value of
     * the following column and so forth.
     * </p>
     *
     * @param value the new value.
     * @return reference to this data-line.
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
     * @return reference to this data-line.
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
     * @return reference to this data-line.
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
     * @return reference to this data-line.
     * @throws IllegalStateException if this operation goes beyond the last column index.
     */
    public DataLine append(final double... values) {
        for (final double d : Utils.nonNull(values, "the values cannot be null")) {
            append(d);
        }
        return this;
    }

    /**
     * Sets all the data-line values at once.
     * @param values the new values.
     * @throws IllegalArgumentException if {@code values} is {@code null}, its length does not match
     *         this data-line column count, or the first element starts like a comment line.
     * @return a reference to this data-line.
     */
    public DataLine setAll(final String... values) {
        Utils.nonNull(values,"the input values cannot be null");
        if (values.length != this.values.length) {
            throw new IllegalArgumentException("the input value length must be equal to the total number of columns ");
        }
        // No index error here this.values.length is guaranteed to be greater than 0, as 0 columns is not allowed.
        if (values[0] != null && values[0].startsWith(TableUtils.COMMENT_PREFIX)) {
            throw new IllegalArgumentException("first column value cannot start as a comment: " + TableUtils.COMMENT_PREFIX);
        }
        System.arraycopy(values,0,this.values,0,values.length);
        return this;
    }

    /**
     * Changes the index of the next value to set using {@link #append append} operations.
     *
     * @param index the new index.
     * @return a reference to this data-line.
     * @throws IllegalArgumentException if {@code index} is greater than the number of columns
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
     * @return a reference to this data-line.
     * @throws IllegalArgumentException if {@code columnName} is {@code null} or an unknown column name.
     */
    public DataLine seek(final String columnName) {
        nextIndex = columnIndex(columnName);
        return this;
    }

    /**
     * Changes the index of the next value to set using {@link #append append} operations.
     * <p>
     * The input enum's string conversion using {@link Object#toString toString} determines the name
     * of the target column.
     * </p>
     * @param column the enum value that makes reference to the target column.
     * @throws IllegalArgumentException if {@code column} is {@code null} or does not make reference to
     *    a known column.
     *
     */
    public DataLine seek(final Enum<?> column) {
        nextIndex = columnIndex(Utils.nonNull(column).toString());
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
