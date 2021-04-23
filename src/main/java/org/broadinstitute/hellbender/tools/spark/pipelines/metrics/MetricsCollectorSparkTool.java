package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.metrics.Header;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.CommandLinePluginDescriptor;
import org.broadinstitute.barclay.argparser.CommandLinePluginProvider;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKReadFilterPluginDescriptor;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.metrics.MetricsArgumentCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.util.Collections;
import java.util.List;

/**
 * Base class for standalone Spark metrics collector tools. Subclasses should
 * adhere to the following pattern:
 *
 * - Declare an instance of a collector that implements the {@link MetricsCollectorSpark} interface
 * - Declare an instance of the input argument collection (of type T) for the collector
 * - Implement or override the {@link #getDefaultReadFilters}, {@link #initialize}, {@link #collectMetrics} and {@link #saveMetrics}
 *   methods by forwarding to the collector object
 *
 * The {link #runTool} method for this class will automatically put the collector through the
 * collection lifecycle.
 */
public abstract class MetricsCollectorSparkTool<T extends MetricsArgumentCollection>
        extends GATKSparkTool {

    private static final long serialVersionUID = 1l;

    /**
     * The following {@link MetricsCollectorSpark} methods must be implemented by subclasses
     * and should be forwarded to the embedded collector.
     */
    @Override
    abstract public List<ReadFilter> getDefaultReadFilters();
    abstract protected SortOrder getExpectedSortOrder();
    abstract protected void initialize(T inputArgs, SAMFileHeader samHeader, List<Header> defaultHeaders);
    abstract protected void collectMetrics(JavaRDD<GATKRead> filteredReads, SAMFileHeader samHeader);
    abstract protected void saveMetrics(String inputBaseName);

    /**
     * To be implemented by subclasses; return the fully initialized and populated
     * argument collection that will be passed to the collector
     */
    abstract protected T getInputArguments();

    @Override
    public final boolean requiresReads(){ return true; }

    /**
     * The runTool method used when metrics collector tools are run
     * "standalone". The filtered RDD is passed to the collectMetrics
     * method which does the bulk of the analysis.
     *
     * @param ctx our Spark context
     */
    @Override
    protected void runTool( JavaSparkContext ctx ) {
        if (requiresReference()) { //TODO: temporary until references are implemented
            throw new UnsupportedOperationException("Requires reference for collector not yet implemented");
        }
        ReadUtils.validateExpectedSortOrder(
                getHeaderForReads().getSortOrder(),
                getExpectedSortOrder(),
                false,
                getReadSourceName()
        );

        // Execute the collector lifecycle
        T collectorArgs = getInputArguments();
        if (null == collectorArgs) {
            // Indicates an incorrectly written collector. Metrics collector arg objects are derived
            // from MetricsArgumentCollection, so they always have at least an output argument.
            throw new IllegalStateException("A Spark metrics collector must return a non-a null argument object");
        }

        initialize(collectorArgs, getHeaderForReads(), getDefaultHeaders());
        final JavaRDD<GATKRead> filteredReads = getReads();
        collectMetrics(filteredReads, getHeaderForReads());
        saveMetrics(getReadSourceName());
    }

}

