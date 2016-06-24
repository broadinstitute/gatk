package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;

import java.io.Serializable;

/**
 *  PerRecordCollector - An interface for classes that collect data in order to generate one or more metrics.
 *      This process usually occurs in the following fashion:
 *      1. Loop through a data set (usually all records in a BAM file) and call collector.acceptRecord( data ),
 *         data in this step is usually added to metrics/Histogram objects
 *      2. Call collector.finish() - perform any final calculations necessary after ALL records have been accepted
 *      3. addMetricsToFile is then used to add any metric(s) or Histogram(s) to the given file
 *
 *      BEAN    - The Metric type we are generating
 *      HKEY    - The Key used in any Histograms, use a Wildcard(?) type if there are no Histograms
 *      ARGTYPE - Collectors are often used in groups of accumulation levels, in order to avoid recalculating
 *                any information needed by multiple collectors we allow different types of arguments that
 *                extend DefaultPerRecordCollectorArgs to accommodate any computed values
 */
public interface PerUnitMetricCollector<BEAN extends MetricBase, HKEY extends Comparable<HKEY>, ARGTYPE> extends Serializable {
    /**
     * Add a SAMRecord (with ReferenceSequence and Read Group info) to the metric(s) being calculated)
     * @param args Contains SAMRecord, SAMReadGroupRecord, ReferenceSequence of current record and any previously
     *             computed values that might be needed for this class
     */
    public void acceptRecord(final ARGTYPE args);

    /** When all records have been collected, compute any final values needed to finish constructing metrics/Histogram */
    public void finish();

    /**
     * Any metrics collected will be added to the metric file provided.
     * @param file MetricsFile to which all metrics created by this collector should be added
     */
    public void addMetricsToFile(final MetricsFile<BEAN, HKEY> file);
}

