package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.metrics.MetricBase;

public abstract class MultiLevelMetrics extends MetricBase {
     /** The sample to which these metrics apply.  If null, it means they apply
     * to all reads in the file. */
    public String SAMPLE;

    /** The library to which these metrics apply.  If null, it means that the
     * metrics were accumulated at the sample level. */
    public String LIBRARY = null;

    /** The read group to which these metrics apply.  If null, it means that
     * the metrics were accumulated at the library or sample level.*/
    public String READ_GROUP = null;
}
