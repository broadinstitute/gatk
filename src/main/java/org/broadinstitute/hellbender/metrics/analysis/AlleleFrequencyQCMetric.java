package org.broadinstitute.hellbender.metrics.analysis;

import htsjdk.samtools.metrics.MetricBase;

public class AlleleFrequencyQCMetric extends MetricBase {

    public String SAMPLE;

    public String METRIC_TYPE;

    public Double METRIC_VALUE;

    public Double CHI_SQ_VALUE;

}
