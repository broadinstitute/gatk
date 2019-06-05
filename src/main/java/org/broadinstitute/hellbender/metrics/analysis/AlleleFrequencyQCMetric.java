package org.broadinstitute.hellbender.metrics.analysis;

import htsjdk.samtools.metrics.MetricBase;

/**
 * Created by skwalker on 5/15/19.
 */
public class AlleleFrequencyQCMetric extends MetricBase {

    /** The name of the sample */
    public String SAMPLE;

    public String METRIC_TYPE;

    public Double METRIC_VALUE;

    public Double CHI_SQ_VALUE;

}
