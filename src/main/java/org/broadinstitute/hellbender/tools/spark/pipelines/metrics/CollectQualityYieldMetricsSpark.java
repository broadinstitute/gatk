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
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.metrics.QualityYieldMetricsArgumentCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;

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
public final class CollectQualityYieldMetricsSpark extends MetricsCollectorSparkTool<QualityYieldMetricsArgumentCollection> {

    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    private QualityYieldMetricsArgumentCollection qualityYieldArgs = new QualityYieldMetricsArgumentCollection();

    private QualityYieldMetricsCollectorSpark qualityYieldCollector = new QualityYieldMetricsCollectorSpark();

    @Override
    public SortOrder getExpectedSortOrder() { return qualityYieldCollector.getExpectedSortOrder(); }

    @Override
    public QualityYieldMetricsArgumentCollection getInputArguments() {
        return qualityYieldArgs;
    }

    @Override
    public ReadFilter getReadFilter(SAMFileHeader samHeader) {
        return ReadFilterLibrary.ALLOW_ALL_READS;
    }

    @Override
    public void initialize(final QualityYieldMetricsArgumentCollection inputArgs, final List<Header> defaultHeaders) {
        qualityYieldCollector.initialize(inputArgs, defaultHeaders);
    }

    @Override
    public void collectMetrics(
            final JavaRDD<GATKRead> filteredReads,
            final SAMFileHeader samHeader,
            final String inputBaseName,
            final AuthHolder authHolder)
    {
        qualityYieldCollector.collectMetrics(
                filteredReads,
                samHeader,
                inputBaseName,
                authHolder
        );
    }

    void finishCollection() {}

}
