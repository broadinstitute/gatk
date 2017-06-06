package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.metrics.Header;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
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
@DocumentedFeature
public final class CollectInsertSizeMetricsSpark
        extends MetricsCollectorSparkTool<InsertSizeMetricsArgumentCollection> {

    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    private InsertSizeMetricsArgumentCollection insertSizeArgs = new InsertSizeMetricsArgumentCollection();

    private InsertSizeMetricsCollectorSpark insertSizeCollector = new InsertSizeMetricsCollectorSpark();

    @Override
    public InsertSizeMetricsArgumentCollection getInputArguments() {
        return insertSizeArgs;
    }

    @Override
    protected SortOrder getExpectedSortOrder() { return insertSizeCollector.getExpectedSortOrder(); }

    @Override
    protected void initialize(
            final InsertSizeMetricsArgumentCollection inputArgs,
            final SAMFileHeader samHeader,
            final List<Header> defaultHeaders)
    {
        insertSizeCollector.initialize(inputArgs, samHeader, defaultHeaders);
    }

    /**
     * Return the read filters required for this collector
     */
    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return insertSizeCollector.getDefaultReadFilters();
    }

    @Override
    protected void collectMetrics(
            final JavaRDD<GATKRead> filteredReads,
            final SAMFileHeader samHeader)
    {
        insertSizeCollector.collectMetrics(filteredReads, samHeader);
    }

    @Override
    protected void saveMetrics(final String inputName) {
        insertSizeCollector.saveMetrics(inputName);
    }
}
