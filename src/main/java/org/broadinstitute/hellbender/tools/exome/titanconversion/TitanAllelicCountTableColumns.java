package org.broadinstitute.hellbender.tools.exome.titanconversion;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Table columns of a TITAN allelic count tab-separated file.
 *
 */
public enum TitanAllelicCountTableColumns {
    CONTIG("Chr"),
    POSITION("Position"),
    REF_NUCLEOTIDE("Ref"),
    REF_COUNT("RefCount"),
    ALT_NUCLEOTIDE("Nref"),
    ALT_COUNT("NrefCount"),
    NORM_QUALITY("NormQuality");

    private final String columnName;  //store the column names

    TitanAllelicCountTableColumns(final String columnName) {  this.columnName = columnName; }

    @Override
    public String toString() {
        return columnName;
    }

    public static final List<String> FULL_COLUMN_NAME_ARRAY = createUnmodifiableList(
            CONTIG, POSITION, REF_NUCLEOTIDE, REF_COUNT, ALT_NUCLEOTIDE, ALT_COUNT, NORM_QUALITY);

    private static List<String> createUnmodifiableList(final TitanAllelicCountTableColumns... columns) {
        return Collections.unmodifiableList(Arrays.asList(columns).stream()
                .map(TitanAllelicCountTableColumns::toString).collect(Collectors.toList()));
    }
}
