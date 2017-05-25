package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.util.stream.Stream;

/**
 * Enumeration of column names that have an special meaning in read-count files.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public enum TargetTableColumn {
    NAME("name", true), /* it was a terrible idea to change this from "NAME" to "name"! */
    CONTIG("contig", true), /* it was a terrible idea to change this from "CONTIG" to "contig"! */
    START("start", true), /* it was a terrible idea to change this from "START" to "start"! */
    END("stop", true), /* it was a terrible idea to change this from "END" to "stop"! */
    GC_CONTENT("GC_CONTENT", false),
    BAIT_COUNT("BAIT_COUNT", false);

    private final String columnName;  //store the column names

    private final boolean mandatory;

    TargetTableColumn(final String columnName, final boolean mandatory) {
        this.columnName = Utils.nonNull(columnName);
        this.mandatory = mandatory;
    }

    @Override
    public String toString() {
        return columnName;
    }

    public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());

    public static final TableColumnCollection MANDATORY_COLUMNS =
            new TableColumnCollection(Stream.of(values()).filter(c -> c.mandatory).map(TargetTableColumn::toString).toArray(String[]::new));

    /**
     * Checks whether a string is a special column's name.
     *
     * @param name never {@code null}.
     * @return {@code true} iff {@code name} maps to any of the special column.
     */
    public static boolean isStandardTargetColumnName(final String name) {
        return COLUMNS.contains(Utils.nonNull(name, "Column name cannot be null."));
    }
}