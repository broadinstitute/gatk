package org.broadinstitute.hellbender.metrics.analysis;

import htsjdk.samtools.metrics.MetricBase;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.util.help.HelpConstants;

@DocumentedFeature(groupName = HelpConstants.DOC_CAT_METRICS, summary = HelpConstants.DOC_CAT_METRICS_SUMMARY)
public class AlleleFrequencyQCMetric extends MetricBase {

    public String SAMPLE;

    public String METRIC_TYPE;

    public Double METRIC_VALUE;

    public Double CHI_SQ_VALUE;

}
