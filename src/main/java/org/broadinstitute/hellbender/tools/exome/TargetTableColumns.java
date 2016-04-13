package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.Collections;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Enumeration of column names that have an special meaning in read-count files.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public enum TargetTableColumns {
    NAME("name", true),
    CONTIG("contig", true),
    START("start", true),
    END("stop", true),
    GC_CONTENT("GC_CONTENT", false),
    REPEAT_FRACTION("REPEAT_FRACTION", false),
    MEAN_COVERAGE("MEAN_COVERAGE", false),
    COVERAGE_VARIANCE("COVERAGE_VARIANCE", false),
    COVERAGE_INTERQUARTILE_RANGE("COVERAGE_INTERQUARTILE_RANGE", false);

    private final String columnName;  //store the column names

    private boolean mandatory;

    TargetTableColumns(final String columnName, final boolean mandatory) {
        this.columnName = Utils.nonNull(columnName);
        this.mandatory = mandatory;
    }

    @Override
    public String toString() {
        return columnName;
    }

    public static final String[] COLUMN_NAME_ARRAY =
            Stream.of(values()).map(TargetTableColumns::toString).toArray(String[]::new);

    /**
     * Set of all strings that correspond to a special column name.
     * <p>
     * Note: this is needed to avoid try-catch the exception thrown by {@link Enum#valueOf}
     * in the implementation of {@code #isStandardTargetColumnName}.
     * </p>
     */
    private static final Set<String> COLUMN_NAME_SET =
            Collections.unmodifiableSet(Stream.of(COLUMN_NAME_ARRAY).collect(Collectors.toSet()));

    /**
     * Set of all mandatory columns.
     */
    public static final Set<TargetTableColumns> MANDATORY_COLUMN_SET =
            Collections.unmodifiableSet(Stream.of(values())
                    .filter(c -> c.mandatory)
                    .collect(Collectors.toSet()));

    /**
     * Array with all the mandatory column names.
     */
    public static final String[] MANDATORY_COLUMN_NAME_ARRAY =
            MANDATORY_COLUMN_SET.stream().map(TargetTableColumns::toString).toArray(String[]::new);
    /**
     * Checks whether a string is a special column's name.
     *
     * @param name never {@code null}.
     * @return {@code true} iff {@code name} maps to any of the special column.
     */
    public static boolean isStandardTargetColumnName(final String name) {
        return COLUMN_NAME_SET.contains(Utils.nonNull(name, "name cannot be null"));
    }
}