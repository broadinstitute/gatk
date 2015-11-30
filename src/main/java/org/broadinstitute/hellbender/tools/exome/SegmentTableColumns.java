package org.broadinstitute.hellbender.tools.exome;



import java.util.stream.Stream;



/**
 * Created by davidben on 11/30/15.
 */
enum SegmentTableColumns {
    SAMPLE("Sample"), CONTIG("Chromosome"), START("Start"), END("End"),
        NUM_PROBES("Num_Probes"), MEAN("Segment_Mean"), CALL("Segment_Call");

    private String columnName;  //store the column names

    SegmentTableColumns(String columnName) {this.columnName = columnName; }

    @Override
    public String toString() {
        return columnName;
    }

    public static final String[] COLUMN_NAME_ARRAY =
            Stream.of(values()).map(SegmentTableColumns::toString).toArray(String[]::new);
}