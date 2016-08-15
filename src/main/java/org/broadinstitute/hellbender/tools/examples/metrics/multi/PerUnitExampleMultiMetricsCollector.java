package org.broadinstitute.hellbender.tools.examples.metrics.multi;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.hellbender.metrics.PerUnitMetricCollector;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.Serializable;
import java.util.Objects;

/**
 * A Collector for individual ExampleMultiMetrics for a given SAMPLE or SAMPLE/LIBRARY or
 * SAMPLE/LIBRARY/READ_GROUP (depending on aggregation levels)
 */
public final class PerUnitExampleMultiMetricsCollector
        implements PerUnitMetricCollector<
            ExampleMultiMetrics,
            Integer,
            PerUnitExampleMultiMetricsCollectorArgs>, Serializable
{
    private static final long serialVersionUID = 1L;

    // the metrics container for this level
    final private ExampleMultiMetrics metrics;

    public PerUnitExampleMultiMetricsCollector(
            final String sample,
            final String library,
            final String readGroup) {
        metrics= new ExampleMultiMetrics(sample, library, readGroup);
    }

    @Override
    public void acceptRecord(final PerUnitExampleMultiMetricsCollectorArgs args) {
        metrics.NUMREADS++;
    }

    @Override
    public void finish() { }

    @Override
    public void addMetricsToFile(final MetricsFile<ExampleMultiMetrics, Integer> file) {
        file.addMetric(metrics);
    }

    /**
     * Combine this PerUnitExampleMultiMetricsCollector with sourceCollector and return a combined
     * PerUnitExampleMultiMetricsCollector.
     * @param sourceCollector PerUnitExampleMultiMetricsCollector to combine in
     * @return PerUnitExampleMultiMetricsCollector representing the combination of the source collector
     * with this collector
     */
    public PerUnitExampleMultiMetricsCollector combine(PerUnitExampleMultiMetricsCollector sourceCollector) {
        Utils.nonNull(sourceCollector);
        final String validationMessage = "Internal error combining collectors";
        validateEquals(this.metrics.SAMPLE, sourceCollector.metrics.SAMPLE, validationMessage);
        validateEquals(this.metrics.LIBRARY, sourceCollector.metrics.LIBRARY, validationMessage);
        validateEquals(this.metrics.READ_GROUP, sourceCollector.metrics.READ_GROUP, validationMessage);

        final PerUnitExampleMultiMetricsCollector combinedCollector = new PerUnitExampleMultiMetricsCollector(
                this.metrics.SAMPLE,
                this.metrics.LIBRARY,
                this.metrics.READ_GROUP);
        combinedCollector.metrics.NUMREADS = this.metrics.NUMREADS + sourceCollector.metrics.NUMREADS;

        return combinedCollector;
    }

    // Safely validate that two strings are equal, even if one or both are null.
    private static void validateEquals(final String source, final String target, final String message) {
        if ( ! Objects.equals(source, target) ) {
            throw new IllegalArgumentException(
                    String.format("%s (%s : %s)",
                            message,
                            source == null ? "null" : source,
                            target == null ? "null" : target)
            );
        }
    }

}
