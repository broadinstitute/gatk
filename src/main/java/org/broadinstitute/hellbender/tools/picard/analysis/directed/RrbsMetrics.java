package org.broadinstitute.hellbender.tools.picard.analysis.directed;

import htsjdk.samtools.metrics.MetricBase;

import java.util.*;

/**
 * Holds per-MetricAccumulationLevel metric information for the RRBS metrics. Required as the MultiLevelCollector
 * is designed around having a single metrics object and we have two being calculated so RrbsMetricsCollector builds
 * this object which can be teased apart downstream
 *
 * NB: This is purely for internal use, if used as a proper metric object it likely won't do what you want it to
 *
 * @author jgentry@broadinstitute.org
 */
final class RrbsMetrics extends MetricBase {
	private final RrbsSummaryMetrics summaryMetrics;
	private final List<RrbsCpgDetailMetrics> detailMetrics;

	public RrbsMetrics(final RrbsSummaryMetrics summaryMetrics, final List<RrbsCpgDetailMetrics> detailMetrics) {
		this.summaryMetrics = summaryMetrics;
		this.detailMetrics = detailMetrics;
	}

	public List<RrbsCpgDetailMetrics> getDetailMetrics() {
		return detailMetrics;
	}

	public RrbsSummaryMetrics getSummaryMetrics() {
		return summaryMetrics;
	}
}
