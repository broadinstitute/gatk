package org.broadinstitute.hellbender.tools.examples.metrics.multi;

import htsjdk.samtools.metrics.MetricBase;
import org.broadinstitute.hellbender.metrics.MultiLevelMetrics;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

/**
 * An example multi-level metrics collector that just counts the number of reads (per unit/level)
 */
public class ExampleMultiMetrics extends MultiLevelMetrics implements Serializable {

    private static final long serialVersionUID = 1;

    // The final metric collected for each level and for all reads. This variable is public, and
    // the name is all caps because it is accessed by the MetricsFile class via reflection, and the
    // variable name will appear in the final MetricsFile. See the javadoc for
    // <a href="http://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/metrics/MetricsFile.html"MetricsFile/a>
    // for more information.
    public int NUMREADS;

    public ExampleMultiMetrics(final String sample, final String library, final String readGroup) {
        this.SAMPLE = sample;
        this.LIBRARY = library;
        this.READ_GROUP = readGroup;
        NUMREADS = 0;
    }

    /**
     * Process a single read
     */
    public void addRead(final GATKRead read){
        NUMREADS++;
    }

    /**
     * Combines two ExampleMultiMetrics objects
     */
    public ExampleMultiMetrics combine(final ExampleMultiMetrics that) {
        Utils.nonNull(that);
        this.NUMREADS += that.NUMREADS;
        return this;
    }

    /**
     * Complete any calculations/processing. No-op for the Example collector.
     */
    public ExampleMultiMetrics finish() {
        return this;
    }

}
