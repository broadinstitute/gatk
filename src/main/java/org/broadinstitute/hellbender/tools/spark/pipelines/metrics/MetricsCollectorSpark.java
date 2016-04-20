package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.metrics.Header;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.metrics.MetricsArgumentCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.List;

/**
 * Interface implemented by Spark metrics collectors to allow them
 * to run as either a standalone Spark tool, or under the control of
 * CollectMultipleMetricsSpark, which reuses the same input RDD to run
 * multiple collectors.
 *
 * Spark collectors should be factored into the following implementation
 * classes:
 *
 * <li>
 * A class derived from MetricsArgumentCollection that defines the input
 * arguments for the collector in a form that is suitable as a command
 * line argument collection.
 * A class that implements this interface, parameterized with the
 * MetricsArgumentCollection-derived class.
 * </li>
 *
 * These two classes can then be used either in the standalone collector tool (see
 * CollectInsertSizeMetricsSpark as an example), or from CollectMultipleMetricsSpark.
 *
 * The general lifecycle of a Spark collector looks like this:
 *
 *     CollectorType collector = new CollectorType<CollectorArgType>()
 *     CollectorArgType args = // get metric-specific input arguments
 *
 *     // pass the input arguments back to the collector for initialization
 *     collector.initialize(args);
 *
 *     ReadFilter filter == collector.getReadFilter();
 *     collector.collectMetrics(
 *         getReads().filter(filter),
 *         getHeaderForReads(),
 *         getReadSourceName(),
 *         getAuthHolder()
 *     );
 *     collector.finishCollection();
 *
 */
public interface MetricsCollectorSpark<T extends MetricsArgumentCollection> extends Serializable
{
    public static final long serialVersionUID = 1L;

    /**
     * Return whether or not this collector requires a reference. Default implementation
     * returns false.
     * @return true if this collector requires that a reference be provided, otherwise false.
     */
    default boolean requiresReference() { return false;}

    /**
     * Return the sort order this collect requires. Collector consumers will validate
     * the expected sort order. If no specific sort order is required, return
     * <code>SortOrder.unsorted</code>. Default implementation returns
     * <code>SortOrder.coordinate</code>
     * @return SortOrder required for this collector or <code>SortOrder.unsorted</code>
     * to indicate any sort order is acceptable.
     */
    default SortOrder getExpectedSortOrder() { return SortOrder.coordinate; }

    /**
     * Give the collector's input arguments to the collector (if the collector
     * is being driven by a standalone tool. This method will always be called
     * before either getReadFilter (so the collector can use the arguments to
     * initialize the readFilter) or collectMetrics.
     * @param inputArgs an object that contains the argument values to be used for
     *                  the collector
     * @param defaultHeaders default metrics headers from the containing tool
     */
    void initialize(T inputArgs, List<Header> defaultHeaders);

    /**
     * Return the read filter required for this collector. The default implementation
     * returns a no-op read filter that allows all reads. Collectors that require
     * a more specific filter criteria should return an instance of the required filter.
     * @param samHeader the SAMFileHeader for the input BAM
     * @return the read filter required for this collector
     */
    default ReadFilter getReadFilter(SAMFileHeader samHeader) {
        return ReadFilterLibrary.ALLOW_ALL_READS;
    }

    /**
     * Do the actual metrics collection on the provided RDD.
     * @param filteredReads The reads to be analyzed for this collector. The reads will have already
     *                      been filtered by this collector's read filter.
     * @param samHeader The SAMFileHeader associated with the reads in the input RDD.
     * @param inputBaseName base name of the input file
     * @param authHolder authHolder
     */
    void collectMetrics(
            JavaRDD<GATKRead> filteredReads,
            SAMFileHeader samHeader,
            String inputBaseName,
            AuthHolder authHolder
    );

    /**
     * This method is called after collectMetrics has returned and
     * all Spark partitions have completed metrics collection. This
     * gives the collector to serialize the resulting metrics, usually
     * to a metrics file.
     */
    default void finishCollection() { }
}
