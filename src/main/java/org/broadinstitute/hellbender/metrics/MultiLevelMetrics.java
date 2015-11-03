package org.broadinstitute.hellbender.metrics;

import org.broadinstitute.hellbender.tools.dataflow.transforms.metrics.SerializableMetric;

import java.util.Comparator;

public abstract class MultiLevelMetrics extends SerializableMetric{
    private static final long serialVersionUID = 1l;

    /** The sample to which these metrics apply.  If null, it means they apply
     * to all reads in the file. */
    public String SAMPLE;

    /** The library to which these metrics apply.  If null, it means that the
     * metrics were accumulated at the sample level. */
    public String LIBRARY;

    /** The read group to which these metrics apply.  If null, it means that
     * the metrics were accumulated at the library or sample level.*/
    public String READ_GROUP;

    /**
     * Generates a comparator to sort by Sample, ReadGroup, and then Library.
     *
     * By specifying T explicitly, this can be used to generate a comparator for a sub class of MultiLevelMetrics which
     * can be chained to produce a comparator sorting first by this ordering and then by whatever additional criteria required.
     *
     * ex:
     *  Comparator<InsertSizeMetrics> comparator = MultiLevelMetrics.<InsertSizeMetrics>getMultiLevelMetricsComparator().thenComparing(insertSizeMetricsComparator);
     */
    public static <T extends MultiLevelMetrics> Comparator<T> getMultiLevelMetricsComparator(){
        return Comparator.comparing((T a) -> a.SAMPLE != null ? a.SAMPLE : "")
            .thenComparing(a -> a.READ_GROUP != null ? a.READ_GROUP : "")
            .thenComparing(a -> a.LIBRARY != null ? a.LIBRARY : "");
    }
}
