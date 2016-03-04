package org.broadinstitute.hellbender.tools.exome;

import java.util.stream.Stream;

/**
 * Created by davidben on 11/30/15.
 */
public enum AllelicCountTableColumns {
    CONTIG("CONTIG"), POSITION("POS"), REF_COUNT("REF_COUNT"), ALT_COUNT("ALT_COUNT");

    private String columnName;

    AllelicCountTableColumns(String columnName) {this.columnName = columnName; }

    @Override
    public String toString() {
        return columnName;
    }

    public static final String[] COLUMN_NAME_ARRAY =
            Stream.of(values()).map(AllelicCountTableColumns::toString).toArray(String[]::new);
}