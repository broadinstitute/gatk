package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Enumeration of column names that have an special meaning in read-count files.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
enum TargetColumns {
    NAME("name"), CONTIG("contig"), START("start"), END("stop");

    private String columnName;  //store the column names

    TargetColumns(String columnName) {this.columnName = columnName; }

    @Override
    public String toString() {
        return columnName;
    }

    public static final String[] COLUMN_NAME_ARRAY =
            Stream.of(values()).map(TargetColumns::toString).toArray(String[]::new);

    /**
     * Set of all strings that correspond to a special column name.
     * <p>
     * Note: this is needed to avoid try-catch the exception thrown by {@link Enum#valueOf}
     * in the implementation of {@code #isTargetColumnName}.
     * </p>
     */
    private static final Set<String> COLUMN_NAME_SET =
            Collections.unmodifiableSet(Stream.of(COLUMN_NAME_ARRAY).collect(Collectors.toSet()));
            //new HashSet<>(Stream.of(values()).map(TargetColumns::name).collect(Collectors.toList()));

    /**
     * Checks whether a string is a special column's name.
     *
     * @param name never {@code null}.
     * @return {@code true} iff {@code name} maps to any of the special column.
     */
    public static boolean isTargetColumnName(final String name) {
        return COLUMN_NAME_SET.contains(Utils.nonNull(name, "name cannot be null"));
    }
}