package org.broadinstitute.hellbender.tools.spark.pathseq.loggers;

import htsjdk.samtools.metrics.MetricBase;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.util.help.HelpConstants;

/**
 * Metrics that are calculated during the PathSeq scoring
 */
@SuppressWarnings("serial")
@DocumentedFeature(groupName = HelpConstants.DOC_CAT_METRICS, summary = HelpConstants.DOC_CAT_METRICS_SUMMARY)
public final class PSScoreMetrics extends MetricBase {

    /**
     * The number of non-host reads mapped to a pathogen
     */
    public Long MAPPED_READS;

    /**
     * The number of non-host reads that could not be mapped to a pathogen
     */
    public Long UNMAPPED_READS;

}
