package org.broadinstitute.hellbender.tools.exome;


import java.util.stream.Stream;



/**
 * Created by davidben on 11/30/15.
 */
enum SegmentTableColumns {
    SAMPLE("Sample"), CONTIG("Chromosome"), START("Start"), END("End"),
    NUM_PROBES("Num_Probes"), MEAN("Segment_Mean"), CALL("Segment_Call"),
    NUM_TARGETS("Num_Targets"), NUM_SNPS("Num_SNPs"),
    SEGMENT_MEAN_POSTERIOR_MEAN("Segment_Mean_Post_Mean"),
    SEGMENT_MEAN_POSTERIOR_STANDARD_DEVIATION("Segment_Mean_Post_Std"),
    MINOR_ALLELE_FRACTION_POSTERIOR_MEAN("MAF_Post_Mean"),
    MINOR_ALLELE_FRACTION_POSTERIOR_STANDARD_DEVIATION("MAF_Post_Std");

    private final String columnName;  //store the column names

    SegmentTableColumns(final String columnName) {  this.columnName = columnName; }

    @Override
    public String toString() {
        return columnName;
    }

    private static final SegmentTableColumns[] INTERVAL_COLUMN_NAME_ENUM_ARRAY =
            new SegmentTableColumns[]{SAMPLE, CONTIG, START, END};

    private static final SegmentTableColumns[] MEAN_AND_CALL_COLUMN_NAME_ENUM_ARRAY =
            new SegmentTableColumns[]{SAMPLE, CONTIG, START, END, NUM_PROBES, MEAN, CALL};

    private static final SegmentTableColumns[] NO_CALL_COLUMN_NAME_ENUM_ARRAY =
            new SegmentTableColumns[]{SAMPLE, CONTIG, START, END, NUM_PROBES, MEAN};

    private static final SegmentTableColumns[] NUM_TARGETS_AND_SNPS_ENUM_ARRAY =
            new SegmentTableColumns[]{SAMPLE, CONTIG, START, END, NUM_TARGETS, NUM_SNPS};

    private static final SegmentTableColumns[] ACS_MODELED_SEGMENT_ENUM_ARRAY =
            new SegmentTableColumns[]{SAMPLE, CONTIG, START, END, NUM_TARGETS, NUM_SNPS,
                    SEGMENT_MEAN_POSTERIOR_MEAN, SEGMENT_MEAN_POSTERIOR_STANDARD_DEVIATION,
                    MINOR_ALLELE_FRACTION_POSTERIOR_MEAN, MINOR_ALLELE_FRACTION_POSTERIOR_STANDARD_DEVIATION};

    public static final String[] INTERVAL_COLUMN_NAME_ARRAY = toStringArray(INTERVAL_COLUMN_NAME_ENUM_ARRAY);

    public static final String[] MEAN_AND_CALL_COLUMN_NAME_ARRAY = toStringArray(MEAN_AND_CALL_COLUMN_NAME_ENUM_ARRAY);

    public static final String[] NO_CALL_COLUMN_NAME_ARRAY = toStringArray(NO_CALL_COLUMN_NAME_ENUM_ARRAY);

    public static final String[] NUM_TARGETS_AND_SNPS_COLUMN_NAME_ARRAY = toStringArray(NUM_TARGETS_AND_SNPS_ENUM_ARRAY);

    public static final String[] ACS_MODELED_SEGMENT_COLUMN_NAME_ARRAY = toStringArray(ACS_MODELED_SEGMENT_ENUM_ARRAY);

    private static String[] toStringArray(final SegmentTableColumns[] enumArray) {
        return Stream.of(enumArray).map(SegmentTableColumns::toString).toArray(String[]::new);
    }
}