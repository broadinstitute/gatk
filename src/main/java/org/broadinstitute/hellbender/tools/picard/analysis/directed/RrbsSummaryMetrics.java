package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import org.broadinstitute.hellbender.metrics.MultiLevelMetrics;

/**
 * Holds summary statistics from RRBS processing QC
 *
 * @author jgentry
 */
public final class RrbsSummaryMetrics extends MultiLevelMetrics {
	/** Number of mapped reads processed */
	public Integer READS_ALIGNED;
	/** Number of times a non-CpG cytosine was encountered */
	public Integer NON_CPG_BASES;
	/** Number of times a non-CpG cytosine was converted (C->T for +, G->A for -) */
	public Integer NON_CPG_CONVERTED_BASES;
	/** NON_CPG_BASES / NON_CPG_CONVERTED_BASES */
	public Double PCT_NON_CPG_BASES_CONVERTED;
	/** Number of CpG sites encountered */
	public Integer CPG_BASES_SEEN;
	/** Number of CpG sites that were converted (TG for +, CA for -) */
	public Integer CPG_BASES_CONVERTED;
	/** CPG_BASES_SEEN / CPG_BASES_CONVERTED */
	public Double PCT_CPG_BASES_CONVERTED;
	/** Mean coverage of CpG sites */
	public Double MEAN_CPG_COVERAGE;
	/** Median coverage of CpG sites */
	public Integer MEDIAN_CPG_COVERAGE;
	/** Number of reads discarded for having no CpG sites */
	public Integer READS_WITH_NO_CPG;
	/** Number of reads discarded due to being too short */
	public Integer READS_IGNORED_SHORT;
	/** Number of reads discarded for exceeding the mismatch threshold */
	public Integer READS_IGNORED_MISMATCHES;
}
