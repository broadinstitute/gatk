package org.broadinstitute.hellbender.tools.dataflow.transforms.metrics;

import htsjdk.samtools.metrics.MetricsFile;

import java.io.Serializable;
import java.io.StringWriter;

/**
 * Serializable version of a {@link MetricsFile} for use with dataflow
 */
public class MetricsFileDataflow<BEAN extends SerializableMetric, HKEY extends Comparable<HKEY> & Serializable> extends MetricsFile<BEAN, HKEY> implements Serializable {
    public static final long serialVersionUID = 1L;

    @Override
    public String toString() {
        StringWriter writer = new StringWriter();  //Create a StringWriter which writes to an in memory buffer
        write(writer); //use MetricFile.write() to generate a string representation of the MetricFile and write it into our writer
        return writer.toString(); //return the string that was written
    }
}
