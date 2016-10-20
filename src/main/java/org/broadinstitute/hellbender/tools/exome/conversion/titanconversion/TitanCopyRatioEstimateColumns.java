package org.broadinstitute.hellbender.tools.exome.conversion.titanconversion;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Table columns of a TITAN allelic count tab-separated file.
 *
 */
public enum TitanCopyRatioEstimateColumns {
    // chr	start	end	log2_TNratio_corrected
    CONTIG("chr"),
    START("start"),
    END("end"),
    TN("log2_TNratio_corrected");

    private final String columnName;  //store the column names

    TitanCopyRatioEstimateColumns(final String columnName) {  this.columnName = columnName; }

    @Override
    public String toString() {
        return columnName;
    }

    public static final List<String> FULL_COLUMN_NAME_ARRAY = createUnmodifiableList(
            CONTIG, START, END, TN);

    private static List<String> createUnmodifiableList(final TitanCopyRatioEstimateColumns... columns) {
        return Collections.unmodifiableList(Arrays.asList(columns).stream()
                .map(TitanCopyRatioEstimateColumns::toString).collect(Collectors.toList()));
    }
}
