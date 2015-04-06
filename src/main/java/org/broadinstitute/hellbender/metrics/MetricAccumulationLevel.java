package org.broadinstitute.hellbender.metrics;

/**
 * For use with Picard metrics programs that may output metrics for multiple levels
 * of aggregation with an analysis.  Used to specify which metrics to output
 */
public enum MetricAccumulationLevel {
    ALL_READS, SAMPLE, LIBRARY, READ_GROUP
}
