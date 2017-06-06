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
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.metrics.QualityYieldMetricsArgumentCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collections;
import java.util.List;

/**
 * Spark tool for collecting quality yield metrics. Delegates to QualityYieldMetricsCollectorSpark.
 */
@CommandLineProgramProperties(
        summary = "Collects quality yield metrics, a set of metrics that quantify the quality and yield of sequence data from a " +
                "SAM/BAM/CRAM input file.",
        oneLineSummary = "CollectQualityYieldMetrics on Spark",
        programGroup = SparkProgramGroup.class
)
@DocumentedFeature
public final class CollectQualityYieldMetricsSpark extends MetricsCollectorSparkTool<QualityYieldMetricsArgumentCollection> {

    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    private QualityYieldMetricsArgumentCollection qualityYieldArgs = new QualityYieldMetricsArgumentCollection();

    private QualityYieldMetricsCollectorSpark qualityYieldCollector = new QualityYieldMetricsCollectorSpark();

    @Override
    protected SortOrder getExpectedSortOrder() { return qualityYieldCollector.getExpectedSortOrder(); }

    @Override
    protected QualityYieldMetricsArgumentCollection getInputArguments() {
        return qualityYieldArgs;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.singletonList(ReadFilterLibrary.ALLOW_ALL_READS);
    }

    @Override
    protected void initialize(
            final QualityYieldMetricsArgumentCollection inputArgs,
            final SAMFileHeader samHeader,
            final List<Header> defaultHeaders) {
        qualityYieldCollector.initialize(inputArgs, samHeader, defaultHeaders);
    }

    @Override
    protected void collectMetrics(
            final JavaRDD<GATKRead> filteredReads,
            final SAMFileHeader samHeader)
    {
        qualityYieldCollector.collectMetrics(filteredReads, samHeader);
    }

    @Override
    protected void saveMetrics(final String inputName) {
        qualityYieldCollector.saveMetrics(inputName);
    }

}
