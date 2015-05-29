package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Enumeration of column names that have an special meaning in read-count files.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
enum ReadCountsSpecialColumns {

    CONTIG, START, END, NAME;

    /**
     * Set of all strings that correspond to a special column name.
     * <p>
     * Note: this is needed to avoid try-catch the exception thrown by {@link Enum#valueOf}
     * in the implementation of {@code #isSpecialColumnName}.
     * </p>
     */
    private static final Set<String> allSpecialColumnNames = new HashSet<>(
            Stream.of(values()).map(ReadCountsSpecialColumns::name).collect(Collectors.toList()));


    /**
     * String containing all special column names separated by a comma and space (", ").
     */
    protected static final String allSpecialColumnNameString =
            Stream.of(values()).map(ReadCountsSpecialColumns::name).collect(Collectors.joining(", "));

    /**
     * Checks whether a string is a special column's name.
     *
     * @param name never {@code null}.
     * @return {@code true} iff {@code name} maps to any of the special column.
     */
    public static boolean isSpecialColumnName(final String name) {
        return allSpecialColumnNames.contains(Utils.nonNull(name, "name cannot be null"));
    }
}