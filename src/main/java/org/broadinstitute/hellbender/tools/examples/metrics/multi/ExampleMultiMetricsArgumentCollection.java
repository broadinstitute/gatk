package org.broadinstitute.hellbender.tools.examples.metrics.multi;

import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.MetricAccumulationLevelArgumentCollection;
import org.broadinstitute.hellbender.metrics.MetricsArgumentCollection;

import java.io.Serializable;

/**
 * Example argument collection for multi-level metrics.
 */
public class ExampleMultiMetricsArgumentCollection extends MetricsArgumentCollection implements Serializable {

    public static final long serialVersionUID = 1L;

    //This example collector only uses the standard argument collection for selecting a collection level.
    @ArgumentCollection
    MetricAccumulationLevelArgumentCollection metricAccumulationLevel = new MetricAccumulationLevelArgumentCollection();
}
