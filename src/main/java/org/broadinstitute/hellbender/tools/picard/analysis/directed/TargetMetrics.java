package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import org.broadinstitute.hellbender.metrics.MultiLevelMetrics;

/**
 * For a sequencing run targeting specific regions of the genome this metric class holds metrics describing
 * how well those regions were targeted.
 */
public final class TargetMetrics extends MultiLevelMetrics {
    /**  The name of the PROBE_SET (BAIT SET, AMPLICON SET, ...) used in this metrics collection run */
    public String PROBE_SET;

    /** The number of unique bases covered by the intervals of all probes in the probe set */
    public long PROBE_TERRITORY;

    /** The number of unique bases covered by the intervals of all targets that should be covered */
    public long TARGET_TERRITORY;

    /** The number of bases in the reference genome used for alignment. */
    public long GENOME_SIZE;

    /** The total number of reads in the SAM/BAM file examined. */
    public long TOTAL_READS;

    /** The number of reads that pass the vendor's filter. */
    public long PF_READS;

    /** The number of bases in the SAM/BAM file to be examined */
    public long PF_BASES;

    /** The number of PF reads that are not marked as duplicates. */
    public long PF_UNIQUE_READS;

    // Tracks the number of read pairs that we see that are PF (used to calculate library size) */
    public long PF_SELECTED_PAIRS;

    // Tracks the number of unique PF reads pairs we see (used to calc library size)
    public long PF_SELECTED_UNIQUE_PAIRS;

    /** The number of PF unique reads that are aligned with mapping score > 0 to the reference genome. */
    public long PF_UQ_READS_ALIGNED;

    /** The number of PF unique bases that are aligned with mapping score > 0 to the reference genome. */
    public long PF_UQ_BASES_ALIGNED;

    /** The number of PF aligned probed that mapped to a baited region of the genome. */
    public long ON_PROBE_BASES;

    /** The number of PF aligned bases that mapped to within a fixed interval of a probed region, but not on a baited region. */
    public long NEAR_PROBE_BASES;

    /** The number of PF aligned bases that mapped to neither on or near a probe. */
    public long OFF_PROBE_BASES;

    /** The number of PF aligned bases that mapped to a targeted region of the genome. */
    public long ON_TARGET_BASES;

    /** The number of PF aligned bases that are mapped in pair to a targeted region of the genome. */
    public long ON_TARGET_FROM_PAIR_BASES;

    //metrics below here are derived after collection

    /** PF reads / total reads.  The percent of reads passing filter. */
    public double PCT_PF_READS;

    /** PF Unique Reads / Total Reads. */
    public double PCT_PF_UQ_READS;

    /** PF Reads Aligned / PF Reads. */
    public double PCT_PF_UQ_READS_ALIGNED;

    /** On+Near Bait Bases / PF Bases Aligned. */
    public double PCT_SELECTED_BASES;

    /** The percentage of aligned PF bases that mapped neither on or near a probe. */
    public double PCT_OFF_PROBE;

    /** The percentage of on+near probe bases that are on as opposed to near. */
    public double ON_PROBE_VS_SELECTED;

    /** The mean coverage of all probes in the experiment. */
    public double MEAN_PROBE_COVERAGE;

    /** The fold by which the probed region has been amplified above genomic background. */
    public double FOLD_ENRICHMENT;

    /** The mean coverage of targets that recieved at least coverage depth = 2 at one base. */
    public double MEAN_TARGET_COVERAGE;

    /** The number of targets that did not reach coverage=2 over any base. */
    public double ZERO_CVG_TARGETS_PCT;

    /**
     * The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to
     * the mean coverage level in those targets.
     */
    public double FOLD_80_BASE_PENALTY;

    /** The percentage of ALL target bases acheiving 2X or greater coverage. */
    public double PCT_TARGET_BASES_2X;
    /** The percentage of ALL target bases acheiving 10X or greater coverage. */
    public double PCT_TARGET_BASES_10X;
    /** The percentage of ALL target bases acheiving 20X or greater coverage. */
    public double PCT_TARGET_BASES_20X;
    /** The percentage of ALL target bases acheiving 30X or greater coverage. */
    public double PCT_TARGET_BASES_30X;
    /** The percentage of ALL target bases acheiving 40X or greater coverage. */
    public double PCT_TARGET_BASES_40X;
    /** The percentage of ALL target bases acheiving 50X or greater coverage. */
    public double PCT_TARGET_BASES_50X;
    /** The percentage of ALL target bases acheiving 100X or greater coverage. */
    public double PCT_TARGET_BASES_100X;

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