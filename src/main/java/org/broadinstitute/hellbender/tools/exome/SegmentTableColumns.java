package org.broadinstitute.hellbender.tools.exome;


import java.util.stream.Stream;



/**
 * Created by davidben on 11/30/15.
 */
enum SegmentTableColumns {
    SAMPLE("Sample"), CONTIG("Chromosome"), START("Start"), END("End"),
    NUM_PROBES("Num_Probes"), MEAN("Segment_Mean"), SD("Segment_Std"), CALL("Segment_Call"),
    NUM_TARGETS("Num_Targets"), NUM_SNPS("Num_SNPs"),
    SEGMENT_MEAN_POSTERIOR_MODE("Segment_Mean_Post_Mode"),
    SEGMENT_MEAN_POSTERIOR_LOWER("Segment_Mean_Post_Lo"),
    SEGMENT_MEAN_POSTERIOR_UPPER("Segment_Mean_Post_Hi"),
    SEGMENT_MEAN_POSTERIOR_10("Segment_Mean_Post_10"),
    SEGMENT_MEAN_POSTERIOR_20("Segment_Mean_Post_20"),
    SEGMENT_MEAN_POSTERIOR_30("Segment_Mean_Post_30"),
    SEGMENT_MEAN_POSTERIOR_40("Segment_Mean_Post_40"),
    SEGMENT_MEAN_POSTERIOR_50("Segment_Mean_Post_50"),
    SEGMENT_MEAN_POSTERIOR_60("Segment_Mean_Post_60"),
    SEGMENT_MEAN_POSTERIOR_70("Segment_Mean_Post_70"),
    SEGMENT_MEAN_POSTERIOR_80("Segment_Mean_Post_80"),
    SEGMENT_MEAN_POSTERIOR_90("Segment_Mean_Post_90"),
    MINOR_ALLELE_FRACTION_POSTERIOR_MODE("MAF_Post_Mode"),
    MINOR_ALLELE_FRACTION_POSTERIOR_LOWER("MAF_Post_Lo"),
    MINOR_ALLELE_FRACTION_POSTERIOR_UPPER("MAF_Post_Hi"),
    MINOR_ALLELE_FRACTION_POSTERIOR_10("MAF_Post_10"),
    MINOR_ALLELE_FRACTION_POSTERIOR_20("MAF_Post_20"),
    MINOR_ALLELE_FRACTION_POSTERIOR_30("MAF_Post_30"),
    MINOR_ALLELE_FRACTION_POSTERIOR_40("MAF_Post_40"),
    MINOR_ALLELE_FRACTION_POSTERIOR_50("MAF_Post_50"),
    MINOR_ALLELE_FRACTION_POSTERIOR_60("MAF_Post_60"),
    MINOR_ALLELE_FRACTION_POSTERIOR_70("MAF_Post_70"),
    MINOR_ALLELE_FRACTION_POSTERIOR_80("MAF_Post_80"),
    MINOR_ALLELE_FRACTION_POSTERIOR_90("MAF_Post_90"),

    // Germline Segments quality scores:

    /**
     * The Phred posterior probability of an event with the same sign as the call that covers all the targets in the
     * segment without interruptions.
     */
    EXACT_QUALITY("Exact_Quality"),

    /**
     * The Phred posterior probability of an event with the same sign as the call that cover at least one target in
     * the segment.
     */
    SOME_QUALITY("Some_Quality"),

    /**
     * The Phred posterior probability of an event with any sign that covers al the target in the segment.
     */
    EVENT_QUALITY("Event_Quality"),

    /**
     * The Phred posterior probability that there is no event of any sign anywhere in the called segment.
     */
    NOEVENT_QUALITY("NoEvent_Quality"),

    /**
     * The Phred posterior probability of an event with the same sign as the call to start exactly at the first
     * target in the segment.
     */
    START_QUALITY("Start_Quality"),

    /**
     * The Phred posterior probability of an event with the same sign as the call to end exactly at the last target
     * in the segment (inclusive).
     */
    END_QUALITY("End_Quality")
    ;

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

    private static final SegmentTableColumns[] ACNV_MODELED_SEGMENT_ENUM_ARRAY =
            new SegmentTableColumns[]{SAMPLE, CONTIG, START, END, NUM_TARGETS, NUM_SNPS,
                    SEGMENT_MEAN_POSTERIOR_MODE, SEGMENT_MEAN_POSTERIOR_LOWER, SEGMENT_MEAN_POSTERIOR_UPPER,
                    SEGMENT_MEAN_POSTERIOR_10, SEGMENT_MEAN_POSTERIOR_20, SEGMENT_MEAN_POSTERIOR_30,
                    SEGMENT_MEAN_POSTERIOR_40, SEGMENT_MEAN_POSTERIOR_50, SEGMENT_MEAN_POSTERIOR_60,
                    SEGMENT_MEAN_POSTERIOR_70, SEGMENT_MEAN_POSTERIOR_80, SEGMENT_MEAN_POSTERIOR_90,
                    MINOR_ALLELE_FRACTION_POSTERIOR_MODE, MINOR_ALLELE_FRACTION_POSTERIOR_LOWER,
                    MINOR_ALLELE_FRACTION_POSTERIOR_UPPER,
                    MINOR_ALLELE_FRACTION_POSTERIOR_10, MINOR_ALLELE_FRACTION_POSTERIOR_20, MINOR_ALLELE_FRACTION_POSTERIOR_30,
                    MINOR_ALLELE_FRACTION_POSTERIOR_40, MINOR_ALLELE_FRACTION_POSTERIOR_50, MINOR_ALLELE_FRACTION_POSTERIOR_60,
                    MINOR_ALLELE_FRACTION_POSTERIOR_70, MINOR_ALLELE_FRACTION_POSTERIOR_80, MINOR_ALLELE_FRACTION_POSTERIOR_90};

    private static final SegmentTableColumns[] GERMLINE_CALL_OUTPUT_COLUMNS =
            new SegmentTableColumns[]{
                    SAMPLE, CONTIG, START, END, NUM_TARGETS, MEAN, SD, EXACT_QUALITY, SOME_QUALITY,
                    EVENT_QUALITY, NOEVENT_QUALITY, START_QUALITY, END_QUALITY
            };

    public static final String[] INTERVAL_COLUMN_NAME_ARRAY = toStringArray(INTERVAL_COLUMN_NAME_ENUM_ARRAY);

    public static final String[] MEAN_AND_CALL_COLUMN_NAME_ARRAY = toStringArray(MEAN_AND_CALL_COLUMN_NAME_ENUM_ARRAY);

    public static final String[] NO_CALL_COLUMN_NAME_ARRAY = toStringArray(NO_CALL_COLUMN_NAME_ENUM_ARRAY);

    public static final String[] NUM_TARGETS_AND_SNPS_COLUMN_NAME_ARRAY = toStringArray(NUM_TARGETS_AND_SNPS_ENUM_ARRAY);

    public static final String[] ACNV_MODELED_SEGMENT_COLUMN_NAME_ARRAY = toStringArray(ACNV_MODELED_SEGMENT_ENUM_ARRAY);

    private static String[] toStringArray(final SegmentTableColumns[] enumArray) {
        return Stream.of(enumArray).map(SegmentTableColumns::toString).toArray(String[]::new);
    }
}