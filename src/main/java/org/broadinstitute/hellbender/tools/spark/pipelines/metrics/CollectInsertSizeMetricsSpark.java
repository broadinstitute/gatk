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
import org.broadinstitute.hellbender.metrics.InsertSizeMetricsArgumentCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;

/**
 * Spark tool for collecting insert size metrics.
 */
@CommandLineProgramProperties(
        summary        = "Program to collect insert size distribution information in SAM/BAM/CRAM file(s)",
        oneLineSummary = "Collect Insert Size Distribution on Spark",
        programGroup   = SparkProgramGroup.class)
public final class CollectInsertSizeMetricsSpark
        extends MetricsCollectorSparkTool<InsertSizeMetricsArgumentCollection> {

    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    private InsertSizeMetricsArgumentCollection insertSizeArgs = new InsertSizeMetricsArgumentCollection();

    private InsertSizeMetricsCollectorSpark insertSizeCollector = new InsertSizeMetricsCollectorSpark();

    public InsertSizeMetricsArgumentCollection getInputArguments() {
        return insertSizeArgs;
    }

    protected SortOrder getExpectedSortOrder() { return insertSizeCollector.getExpectedSortOrder(); }

    protected void initialize(
            final InsertSizeMetricsArgumentCollection inputArgs,
            final SAMFileHeader samHeader,
            final List<Header> defaultHeaders)
    {
        insertSizeCollector.initialize(inputArgs, samHeader, defaultHeaders);
    }

    /**
     * Return the read filter required for this collector
     */
    @Override
    protected ReadFilter getReadFilter(final SAMFileHeader samHeader) {
        return insertSizeCollector.getReadFilter(samHeader);
    }

    @Override
    protected void collectMetrics(
            final JavaRDD<GATKRead> filteredReads,
            final SAMFileHeader samHeader)
    {
        insertSizeCollector.collectMetrics(filteredReads, samHeader);
    }

    @Override
    protected void finish(final String inputName, final AuthHolder authHolder) {
        insertSizeCollector.saveMetrics(inputName, authHolder);
    }
}
