package org.broadinstitute.hellbender.tools.examples.metrics.single;

import htsjdk.samtools.metrics.MetricBase;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;

public class ExampleSingleMetrics extends MetricBase implements Serializable {
    private static final long serialVersionUID = 1;

    // The final metric collected for each level and for all reads. This variable is public, and
    // the name is all caps because it is accessed by the MetricsFile class via reflection, and the
    // variable name will appear in the final MetricsFile. See the javadoc for
    // <a href="http://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/metrics/MetricsFile.html"MetricsFile/a>
    // for more information.
    public int NUMREADS = 0;

    //adds stats for the read
    public ExampleSingleMetrics addRead(final GATKRead read) {
        NUMREADS++;
        return this;
    }

    //combines two objects
    public ExampleSingleMetrics combine(final ExampleSingleMetrics that) {
        this.NUMREADS += that.NUMREADS;
        return this;
    }

    //completes the calculations - no-op
    public ExampleSingleMetrics finish() {
        return this;
    }

}
