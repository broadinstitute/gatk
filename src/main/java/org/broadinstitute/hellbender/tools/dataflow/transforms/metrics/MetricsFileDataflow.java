package org.broadinstitute.hellbender.tools.dataflow.transforms.metrics;

import htsjdk.samtools.metrics.MetricsFile;

import java.io.Serializable;
import java.io.StringWriter;

/**
 * Serializable version of a {@link MetricsFile} for use with dataflow
 */
public class MetricsFileDataflow<BEAN extends SerializableMetric, HKEY extends Comparable<HKEY> & Serializable> extends MetricsFile<BEAN, HKEY> implements Serializable {
    public static final long serialVersionUID = 1l;

    @Override
    public String toString() {
        StringWriter writer = new StringWriter();
        write(writer);
        return writer.toString();
    }
}
