package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import org.broadinstitute.hellbender.metrics.MultiLevelMetrics;

/**
 * Holds information about CpG sites encountered for RRBS processing QC
 * @author jgentry
 */
public final class RrbsCpgDetailMetrics extends MultiLevelMetrics {
	/** Sequence the CpG is seen in */
	public String SEQUENCE_NAME;
	/** Position within the sequence of the CpG site */
	public Integer POSITION;
	/** Number of times this CpG site was encountered */
	public Integer TOTAL_SITES;
	/** Number of times this CpG site was converted (TG for + strand, CA for - strand) */
	public Integer CONVERTED_SITES;
	/** TOTAL_BASES / CONVERTED_BASES */
	public Double PCT_CONVERTED;
}
