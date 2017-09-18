package org.broadinstitute.hellbender.tools.spark.pathseq;

import htsjdk.samtools.metrics.MetricBase;

import java.io.Serializable;

/**
 * Metrics that are calculated during the PathSeq filter
 */
@SuppressWarnings("serial")
public final class PSScoreMetrics extends MetricBase implements Serializable {

    /**
     * The number of non-host reads mapped to a pathogen
     */
    public Long MAPPED_READS;

    /**
     * The number of non-host reads that could not be mapped to a pathogen
     */
    public Long UNMAPPED_READS;

}
