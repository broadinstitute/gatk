package org.broadinstitute.hellbender.tools.exome.titanconversion;

import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.util.stream.Stream;

/**
 * Table columns of a TITAN allelic count tab-separated file.
 */
public enum TitanAllelicCountTableColumn {
    CONTIG("Chr"),
    POSITION("Position"),
    REF_NUCLEOTIDE("Ref"),
    REF_COUNT("RefCount"),
    ALT_NUCLEOTIDE("Nref"),
    ALT_COUNT("NrefCount"),
    NORM_QUALITY("NormQuality");

    private final String columnName;  //store the column names

    TitanAllelicCountTableColumn(final String columnName) {  this.columnName = columnName; }

    @Override
    public String toString() {
        return columnName;
    }

    public static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
}
