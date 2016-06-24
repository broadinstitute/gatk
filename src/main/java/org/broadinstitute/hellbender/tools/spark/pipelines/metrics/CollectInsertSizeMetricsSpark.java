package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.metrics.Header;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;

/**
 * Spark tool for collecting insert size metrics. Delegates to InsertSizeMetricsCollectorSpark.
 */
@CommandLineProgramProperties(
        summary        = "Program to collect insert size distribution information in SAM/BAM/CRAM file(s)",
        oneLineSummary = "Collect Insert Size Distribution on Spark",
        programGroup   = SparkProgramGroup.class)
public final class CollectInsertSizeMetricsSpark
        extends MetricsCollectorSparkTool<InsertSizeMetricsArgumentCollection> {

    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    InsertSizeMetricsArgumentCollection insertSizeCollectorArgs = new InsertSizeMetricsArgumentCollection();

    InsertSizeMetricsCollectorSpark insertSizeCollector = new InsertSizeMetricsCollectorSpark();

    @Override
    public SortOrder getExpectedSortOrder() { return insertSizeCollector.getExpectedSortOrder(); }

    @Override
    public InsertSizeMetricsArgumentCollection getInputArguments() {
        return insertSizeCollectorArgs;
    }

    /**
     * Expose the read filter required for this collector
     */
    @Override
    public ReadFilter getReadFilter(final SAMFileHeader samHeader) {
        return insertSizeCollector.getReadFilter(samHeader);
    }

    public void initialize(
            final InsertSizeMetricsArgumentCollection inputArgs,
            final List<Header> defaultHeaders) {
        insertSizeCollector.initialize(inputArgs, defaultHeaders);
    }

    @Override
    public void collectMetrics(
            final JavaRDD<GATKRead> filteredReads,
            final  SAMFileHeader samHeader,
            final String inputBaseName,
            final AuthHolder authHolder)
    {
        insertSizeCollector.collectMetrics(
                filteredReads,
                samHeader,
                inputBaseName,
                authHolder
        );
    }

    void finishCollection() { }

}
