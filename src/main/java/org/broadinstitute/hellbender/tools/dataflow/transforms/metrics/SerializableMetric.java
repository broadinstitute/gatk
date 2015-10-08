package org.broadinstitute.hellbender.tools.dataflow.transforms.metrics;

import htsjdk.samtools.metrics.MetricBase;

import java.io.Serializable;

/**
 * Serializable Metric base class
 */
public abstract class SerializableMetric extends MetricBase implements Serializable {
    private static final long serialVersionUID = 1l;
}
