package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import org.broadinstitute.hellbender.metrics.MultiLevelMetrics;

/** Metrics class for targeted pcr runs such as TSCA runs */
public final class TargetedPcrMetrics extends MultiLevelMetrics {

    /**  The name of the amplicon set used in this metrics collection run */
    public String CUSTOM_AMPLICON_SET;

    /** The number of bases in the reference genome used for alignment. */
    public long GENOME_SIZE;

    /** The number of unique bases covered by the intervals of all amplicons in the amplicon set */
    public long AMPLICON_TERRITORY;

    /** The number of unique bases covered by the intervals of all targets that should be covered */
    public long TARGET_TERRITORY;

    /** The total number of reads in the SAM/BAM file examine. */
    public long TOTAL_READS;

    /** The number of reads that pass the vendor's filter. */
    public long PF_READS;

    /** THe number of bases in the SAM/BAM file to be examined */
    public long PF_BASES;

    /** The number of PF reads that are not marked as duplicates. */
    public long PF_UNIQUE_READS;

    /** PF reads / total reads.  The percent of reads passing filter. */
    public double PCT_PF_READS;

    /** PF Unique Reads / Total Reads. */
    public double PCT_PF_UQ_READS;

    /** The number of PF unique reads that are aligned with mapping score > 0 to the reference genome. */
    public long PF_UQ_READS_ALIGNED;

    /** Tracks the number of read pairs that we see that are PF (used to calculate library size) */
    public long PF_SELECTED_PAIRS;

    /** Tracks the number of unique PF reads pairs we see (used to calc library size) */
    public long PF_SELECTED_UNIQUE_PAIRS;

    /** PF Reads Aligned / PF Reads. */
    public double PCT_PF_UQ_READS_ALIGNED;

    /** The number of PF unique bases that are aligned with mapping score > 0 to the reference genome. */
    public long PF_UQ_BASES_ALIGNED;

    /** The number of PF aligned amplified that mapped to an amplified region of the genome. */
    public long ON_AMPLICON_BASES;

    /** The number of PF aligned bases that mapped to within a fixed interval of an amplified region, but not on a baited region. */
    public long NEAR_AMPLICON_BASES;

    /** The number of PF aligned bases that mapped to neither on or near an amplicon. */
    public long OFF_AMPLICON_BASES;

    /** The number of PF aligned bases that mapped to a targeted region of the genome. */
    public long ON_TARGET_BASES;

    /** The number of PF aligned bases that are mapped in pair to a targeted region of the genome. */
    public long ON_TARGET_FROM_PAIR_BASES;

    /** On+Near Amplicon Bases / PF Bases Aligned. */
    public double PCT_AMPLIFIED_BASES;

    /** The percentage of aligned PF bases that mapped neither on or near an amplicon. */
    public double PCT_OFF_AMPLICON;

    /** The percentage of on+near amplicon bases that are on as opposed to near. */
    public double ON_AMPLICON_VS_SELECTED;

    /** The mean coverage of all amplicons in the experiment. */
    public double MEAN_AMPLICON_COVERAGE;

    /** The mean coverage of targets that recieved at least coverage depth = 2 at one base. */
    public double MEAN_TARGET_COVERAGE;

    /** The fold by which the amplicon region has been amplified above genomic background. */
    public double FOLD_ENRICHMENT;

    /** The number of targets that did not reach coverage=2 over any base. */
    public double ZERO_CVG_TARGETS_PCT;

    /**
     * The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to
     * the mean coverage level in those targets.
     */
    public double FOLD_80_BASE_PENALTY;

    /** The percentage of ALL target bases achieving 2X or greater coverage. */
    public double PCT_TARGET_BASES_2X;
    /** The percentage of ALL target bases achieving 10X or greater coverage. */
    public double PCT_TARGET_BASES_10X;
    /** The percentage of ALL target bases achieving 20X or greater coverage. */
    public double PCT_TARGET_BASES_20X;
	/** The percentage of ALL target bases achieving 30X or greater coverage. */
	public double PCT_TARGET_BASES_30X;

    /**
     * A measure of how undercovered <= 50% GC regions are relative to the mean. For each GC bin [0..50]
     * we calculate a = % of target territory, and b = % of aligned reads aligned to these targets.
     * AT DROPOUT is then abs(sum(a-b when a-b < 0)). E.g. if the value is 5% this implies that 5% of total
     * reads that should have mapped to GC<=50% regions mapped elsewhere.
     */
    public double AT_DROPOUT;

    /**
     * A measure of how undercovered >= 50% GC regions are relative to the mean. For each GC bin [50..100]
     * we calculate a = % of target territory, and b = % of aligned reads aligned to these targets.
     * GC DROPOUT is then abs(sum(a-b when a-b < 0)). E.g. if the value is 5% this implies that 5% of total
     * reads that should have mapped to GC>=50% regions mapped elsewhere.
     */
    public double GC_DROPOUT;
}
