package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import org.broadinstitute.hellbender.metrics.MultiLevelMetrics;

/**
 * The set of metrics captured that are specific to a hybrid selection analysis.
 *
 * @author Tim Fennell
 */
public final class HsMetrics extends MultiLevelMetrics {
    /** The name of the bait set used in the hybrid selection. */
    public String BAIT_SET;

    /** The number of bases in the reference genome used for alignment. */
    public long GENOME_SIZE;

    /** The number of bases which have one or more baits on top of them. */
    public long BAIT_TERRITORY;

    /** The unique number of target bases in the experiment where target is usually targets etc. */
    public long TARGET_TERRITORY;

    /** Target terrirtoy / bait territory.  1 == perfectly efficient, 0.5 = half of baited bases are not target. */
    public double BAIT_DESIGN_EFFICIENCY;

    /** The total number of reads in the SAM/BAM file examine. */
    public long TOTAL_READS;

    /** The number of reads that pass the vendor's filter. */
    public long PF_READS;

    /** The number of PF reads that are not marked as duplicates. */
    public long PF_UNIQUE_READS;

    /** PF reads / total reads.  The percent of reads passing filter. */
    public double PCT_PF_READS;

    /** PF Unique Reads / Total Reads. */
    public double PCT_PF_UQ_READS;

    /** The number of PF unique reads that are aligned with mapping score > 0 to the reference genome. */
    public long PF_UQ_READS_ALIGNED;

    /** PF Reads Aligned / PF Reads. */
    public double PCT_PF_UQ_READS_ALIGNED;

    /** The number of bases in the PF aligned reads that are mapped to a reference base. Accounts for clipping and gaps. */
    public long PF_UQ_BASES_ALIGNED;

    /** The number of PF aligned bases that mapped to a baited region of the genome. */
    public long ON_BAIT_BASES;

    /** The number of PF aligned bases that mapped to within a fixed interval of a baited region, but not on a baited region. */
    public long NEAR_BAIT_BASES;

    /** The number of PF aligned bases that mapped to neither on or near a bait. */
    public long OFF_BAIT_BASES;

    /** The number of PF aligned bases that mapped to a targeted region of the genome. */
    public long ON_TARGET_BASES;

    /** On+Near Bait Bases / PF Bases Aligned. */
    public double PCT_SELECTED_BASES;

    /** The percentage of aligned PF bases that mapped neither on or near a bait. */
    public double PCT_OFF_BAIT;

    /** The percentage of on+near bait bases that are on as opposed to near. */
    public double ON_BAIT_VS_SELECTED;

    /** The mean coverage of all baits in the experiment. */
    public double MEAN_BAIT_COVERAGE;

    /** The mean coverage of targets that received at least coverage depth = 2 at one base. */
    public double MEAN_TARGET_COVERAGE;

    /** The number of aligned, de-duped, on-bait bases out of the PF bases available. */
    public double PCT_USABLE_BASES_ON_BAIT;

    /** The number of aligned, de-duped, on-target bases out of the PF bases available. */
    public double PCT_USABLE_BASES_ON_TARGET;

    /** The fold by which the baited region has been amplified above genomic background. */
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
	/** The percentage of ALL target bases achieving 40X or greater coverage. */
	public double PCT_TARGET_BASES_40X;
	/** The percentage of ALL target bases achieving 50X or greater coverage. */
	public double PCT_TARGET_BASES_50X;
	/** The percentage of ALL target bases achieving 100X or greater coverage. */
	public double PCT_TARGET_BASES_100X;

    /** The estimated number of unique molecules in the selected part of the library. */
    public Long HS_LIBRARY_SIZE;

    /**
     * The "hybrid selection penalty" incurred to get 80% of target bases to 10X. This metric
     * should be interpreted as: if I have a design with 10 megabases of target, and want to get
     * 10X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 10 * HS_PENALTY_10X.
     */
    public double HS_PENALTY_10X;

    /**
     * The "hybrid selection penalty" incurred to get 80% of target bases to 20X. This metric
     * should be interpreted as: if I have a design with 10 megabases of target, and want to get
     * 20X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 20 * HS_PENALTY_20X.
     */
    public double HS_PENALTY_20X;

	/**
	 * The "hybrid selection penalty" incurred to get 80% of target bases to 30X. This metric
	 * should be interpreted as: if I have a design with 10 megabases of target, and want to get
	 * 30X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 30 * HS_PENALTY_30X.
	 */
	public double HS_PENALTY_30X;

	/**
	 * The "hybrid selection penalty" incurred to get 80% of target bases to 40X. This metric
	 * should be interpreted as: if I have a design with 10 megabases of target, and want to get
	 * 40X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 40 * HS_PENALTY_40X.
	 */
	public double HS_PENALTY_40X;

	/**
	 * The "hybrid selection penalty" incurred to get 80% of target bases to 50X. This metric
	 * should be interpreted as: if I have a design with 10 megabases of target, and want to get
	 * 50X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 50 * HS_PENALTY_50X.
	 */
	public double HS_PENALTY_50X;

	/**
	 * The "hybrid selection penalty" incurred to get 80% of target bases to 100X. This metric
	 * should be interpreted as: if I have a design with 10 megabases of target, and want to get
	 * 100X coverage I need to sequence until PF_ALIGNED_BASES = 10^7 * 100 * HS_PENALTY_100X.
	 */
	public double HS_PENALTY_100X;

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
