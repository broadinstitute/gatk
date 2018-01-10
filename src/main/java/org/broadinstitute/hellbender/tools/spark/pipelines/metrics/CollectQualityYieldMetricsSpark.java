package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.metrics.Header;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.metrics.QualityYieldMetricsArgumentCollection;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Collections;
import java.util.List;

/**
 * Collects quality yield metrics in SAM/BAM/CRAM file(s). The tool leverages the Spark framework for faster
 * operation.
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk CollectQualityYieldMetricsSpark \
 *   -R reference.fasta \
 *   -I input.bam \
 *   -O quality_yield_metrics.txt
 * </pre>
 *
 */@CommandLineProgramProperties(
        summary = "Collects quality yield metrics from SAM/BAM/CRAM file(s). The tool leverages the Spark " +
                "framework for faster operation.",
        oneLineSummary = "Collects quality yield metrics from SAM/BAM/CRAM file(s).",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
@BetaFeature
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
