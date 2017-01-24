package org.broadinstitute.hellbender.tools.copynumber.coverage.readcount;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

/**
 * Enumeration of required columns in read count files
 *
 * @author Andrey Smirnov &lt;asmirnov@broadinstitute.org&gt;
 */
public enum ReadCountTableColumn {
    CONTIG("CONTIG"),
    START("START"),
    END("STOP");

    private final String columnName;  //store the column names

    ReadCountTableColumn(final String columnName) {
        this.columnName = Utils.nonNull(columnName);
    }

    @Override
    public String toString() {
        return columnName;
    }

    public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());

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
