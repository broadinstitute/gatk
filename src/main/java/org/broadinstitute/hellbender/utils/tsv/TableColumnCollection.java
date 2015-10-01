package org.broadinstitute.hellbender.utils.tsv;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Represents a list of table columns.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class TableColumnCollection {

    /**
     * List of column names sorted by their index.
     */
    private final List<String> names;

    /**
     * Map from column name to its index, the first column has index 0, the second 1 and so forth.
     */
    private final Map<String, Integer> indexByName;

    /**
     * Creates a new table-column names collection.
     * <p>
     * The new instance will have its own copy of the input name array,
     * thus is safe to modify such an array after calling this constructor;
     * it won't modify the state of this object.
     * </p>
     *
     * @param names the column names.
     * @throws IllegalArgumentException if {@code names} is not a valid column name iterable
     *                                  of names as determined by {@link #checkNames checkNames}.
     */
    public TableColumnCollection(final Iterable<String> names) {
        this(StreamSupport.stream(Utils.nonNull(names, "the names cannot be null").spliterator(), false).toArray(String[]::new));
    }

    /**
     * Creates a new table-column names collection.
     * <p>
     * The new instance will have its own copy of the input name array,
     * thus is safe to modify such an array after calling this constructor;
     * it won't modify the state of this object.
     * </p>
     *
     * @param names the column names.
     * @throws IllegalArgumentException if {@code names} is not a valid column name array
     *                                  of names as determined by {@link #checkNames checkNames}.
     */
    public TableColumnCollection(final String... names) {
        this.names = Collections.unmodifiableList(Arrays.asList(checkNames(names.clone(), IllegalArgumentException::new)));
        this.indexByName = IntStream.range(0, names.length).boxed()
                .collect(Collectors.toMap(this.names::get, Function.identity()));
    }

    /**
     * Creates a new table-column names collection.
     * <p>
     * The new instance will have its own copy of the input name array,
     * thus is safe to modify such an array after calling this constructor;
     * it won't modify the state of this object.
     * </p>
     *
     * @param names the column names. They will be transformed into strings using {@link Object#toString() toString}.
     * @throws IllegalArgumentException if {@code names} is not a valid column name array
     *                                  of names as determined by {@link #checkNames checkNames}.
     */
    public TableColumnCollection(final Object... names) {
        this.names = Collections.unmodifiableList(Arrays.asList(checkNames(names, IllegalArgumentException::new)));
        this.indexByName = IntStream.range(0, names.length).boxed()
                .collect(Collectors.toMap(this.names::get, Function.identity()));
    }

    /**
     * Creates a new table-column collection from a set of enum constants.
     *
     * <p>
     * The new instance will have one column for enum value in {@code enumClass}, in their
     * ordinal order.
     * </p>
     *
     * <p>
     * The {@link Object#toString() toString} transformation of each constant is used
     * as the column name. Notice that this might differ from the constant name if
     * this method is overloaded by the enum.
     * </p>
     * @param enumClass the input enum class.
     * @throws IllegalArgumentException if {@code enumClass} is {@code null} or is actually
     *  not an enum class or the corresponding names are illegal.
     */
    public TableColumnCollection(final Class<? extends Enum<?>> enumClass) {
        Utils.nonNull(enumClass);
        // Despite the generic annotation, due to erasure this might not be a enum class
        // in run-time.
        if (!enumClass.isEnum()) {
            throw new IllegalArgumentException("the input class must be an enum class");
        }
        names = Collections.unmodifiableList(
                Stream.of(enumClass.getEnumConstants())
                .map(Object::toString)
                .collect(Collectors.toList()));
        checkNames(names.toArray(new String[names.size()]), IllegalArgumentException::new);
        indexByName = IntStream.range(0, names.size()).boxed()
                .collect(Collectors.toMap(names::get, Function.identity()));
    }

    /**
     * Returns the column names ordered by column index.
     *
     * @return never {@code null}, a unmodifiable view to this collection column names.
     */
    public List<String> names() {
        return names;
    }

    /**
     * Returns the name of a column by its index.
     * <p>
     * Column indexes are 0-based.
     * </p>
     *
     * @param index the target column index.
     * @return never {@code null}.
     * @throws IllegalArgumentException if {@code index} is not a valid column index.
     */
    public String nameAt(final int index) {
        Utils.validIndex(index, names.size());
        return names.get(index);
    }

    /**
     * Returns the index of a column by its name.
     *
     * @param name the query column name.
     * @return {@code -1} if there is not such a column, 0 or greater otherwise.
     */
    public int indexOf(final String name) {
        Utils.nonNull(name, "the column name cannot be null");
        return indexByName.getOrDefault(name, -1);
    }

    /**
     * Check whether there is such a column by name.
     *
     * @param name the query name.
     * @return {@code true} iff there is a column with such a {@code name}.
     * @throws IllegalArgumentException if {@code name} is {@code null}.
     */
    public boolean contains(final String name) {
        return indexByName.containsKey(Utils.nonNull(name, "cannot be null"));
    }

    /**
     * Checks whether columns contain all the names in an array.
     *
     * @param names the names to test.
     * @return {@code true} iff all the names in {@code names} correspond to columns in this collection.
     * @throws IllegalArgumentException if {@code names} is {@code null} or contains any {@code null}.
     */
    public boolean containsAll(final String... names) {
        return Stream.of(Utils.nonNull(names, "names cannot be null")).allMatch(this::contains);
    }

    /**
     * Checks whether columns contain all the names in an array.
     *
     * @param names the names to test.
     * @return {@code true} iff all the names in {@code names} correspond to columns in this collection.
     * @throws IllegalArgumentException if {@code names} is {@code null} or contains any {@code null}.
     */
    public boolean containsAll(final Iterable<String> names) {
        for (final String name : Utils.nonNull(names, "names cannot be null")) {
            if (!contains(name)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Checks whether columns contain all the names in an array and no other.
     *
     * @param names the names to test.
     * @return {@code true} iff all the names in {@code names} correspond to columns in this collection, and there
     * is no other column name.
     * @throws IllegalArgumentException if {@code names} is {@code null} or contains {@code null}.
     */
    public boolean containsExactly(final String... names) {
        return containsAll(names) && columnCount() == names.length;
    }

    /**
     * Checks a column names matches given one.
     * <p>
     * A match against a column index beyond the last column returns {@code false}.
     * </p>
     *
     * @param index the column index to test.
     * @param name  the expected column name.
     * @return {@code true} iff the the column at {@code index}'s name is equal to the input {@code name}.
     * @throws IllegalArgumentException if {@code name} is {@code null} or {@code index} is negative.
     */
    public boolean matches(final int index, final String name) {
        if (index >= names.size()) {
            return false;
        } else {
            return names.get(Utils.validIndex(index, names.size())).equals(Utils.nonNull(name, "name cannot be null"));
        }
    }

    /**
     * Checks an array of column matches the column names starting from a given index.
     *
     * @param offset the column index to match to the first name in {@code names} test.
     * @param names  the expected column names.
     * @return {@code true} iff the the column names starting from {@code index} is equal to the input {@code names}.
     * @throws IllegalArgumentException if {@code names} is {@code null}, contains any {@code null}
     *                                  or {@code index} is not a valid index.
     */
    public boolean matchesAll(final int offset, final String... names) {
        Utils.validIndex(offset, this.names.size() + 1);
        final int toIndex = offset + names.length;
        if (toIndex > this.names.size()) {
            return false;
        }
        for (int i = 0; i < names.length; i++) {
            if (!this.names.get(offset + i).equals(names[i])) {
                return false;
            }
        }
        return true;
    }

    /**
     * Checks whether columns contain all the names in an array and no other and in the very same order.
     *
     * @param names the names to test.
     * @return {@code true} iff all the names in {@code names} correspond to columns in this collection in the same
     * order, and there is no other column.
     * @throws IllegalArgumentException if {@code names} is {@code null} or contains {@code null}.
     */
    public boolean matchesExactly(final String... names) {
        return matchesAll(0, names) && names.length == this.names.size();
    }

    /**
     * Returns the number of columns.
     *
     * @return never {@code column}
     */
    public int columnCount() {
        return names.size();
    }


    /**
     * Checks that a column name, as objects, array is valid.
     * <p>
     * The actual column names to test are obtained by mapping its object to its string representation using
     * {@link Object#toString() toString}.
     * </p>
     *
     * <p>
     * Assuming that it is null-value free, a column name array is invalid if:
     * <ul>
     * <li>has length 0,</li>
     * <li>the first name contains the comment prefix,</li>
     * <li>or contains repeats</li>
     * </ul>
     * </p>
     * <p>When the input array is invalid, an exception is thrown using the exception factory function provided.</p>
     * <p>The message passed to the factory explains why the column name array is invalid</p>.
     * <p>Notice that a {@code null} array or a {@code null} containing array is considered a illegal argument caused by
     * a bug and a {@code IllegalArgumentException} will be thrown instead.</p>
     *
     * @throws IllegalArgumentException if {@code columnNames} is {@code null}, or it contains any {@code null}, or {@code exceptionFactory} is {@code null} or ir returns a {@code null}
     *                                  when invoked.
     * @throws RuntimeException         if {@code columnNames} does not contain a valid list of column names. The exact type will depend on the
     *                                  input {@code exceptionFactory}.
     * @return never {@code null}, the same reference as the input column name array {@code columnNames}.
     */
    public static String[] checkNames(final Object[] columnNames,
                                      final Function<String, RuntimeException> exceptionFactory) {
        Utils.nonNull(columnNames, "column names cannot be null");
        final String[] stringNames = new String[columnNames.length];
        for (int i = 0; i < columnNames.length; i++) {
            stringNames[i] = Utils.nonNull(columnNames[i],"no column name can be null: e.g. " + i + " element").toString();
        }
        return checkNames(stringNames, exceptionFactory);
    }

    /**
     * Checks that a column name array is valid.
     *
     * <p>
     * Assuming that it is null-value free, a column name array is invalid if:
     * <ul>
     * <li>has length 0,</li>
     * <li>the first name contains the comment prefix,</li>
     * <li>or contains repeats</li>
     * </ul>
     * </p>
     * <p>When the input array is invalid, an exception is thrown using the exception factory function provided.</p>
     * <p>The message passed to the factory explains why the column name array is invalid</p>.
     * <p>Notice that a {@code null} array or a {@code null} containing array is considered a illegal argument caused by
     * a bug an a {@code IllegalArgumentException} will be thrown instead.</p>
     *
     * @throws IllegalArgumentException if {@code columnNames} is {@code null}, or it contains any {@code null}, or {@code exceptionFactory} is {@code null} or ir returns a {@code null}
     *                                  when invoked.
     * @throws RuntimeException         if {@code columnNames} does not contain a valid list of column names. The exact type will depend on the
     *                                  input {@code exceptionFactory}.
     * @return never {@code null}, the same reference as the input column name array {@code columnNames}.
     */
    public static String[] checkNames(final String[] columnNames,
                                       final Function<String, RuntimeException> exceptionFactory) {
        Utils.nonNull(columnNames, "column names cannot be null");
        Utils.nonNull(exceptionFactory, "exception factory cannot be null");

        if (columnNames.length == 0) {
            throw Utils.nonNull(exceptionFactory.apply("there must be at least one column"));
        }
        final Set<String> columnNameSet = new HashSet<>(columnNames.length);
        for (int i = 0; i < columnNames.length; i++) {
            final String columnName = Utils.nonNull(columnNames[i],"no column name can be null: e.g. " + i + " element");
            if (!columnNameSet.add(columnName)) {
                throw Utils.nonNull(exceptionFactory.apply("more than one column have the same name: " + columnNames[i]), "exception factory produces null exceptions");
            }
        }
        if (columnNames[0].startsWith(TableUtils.COMMENT_PREFIX)) {
            throw Utils.nonNull(exceptionFactory.apply("the first column name cannot start with the comment prefix"), "exception factory produces null exceptions");
        }
        return columnNames;
    }
}
